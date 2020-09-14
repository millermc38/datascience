#' @export datascience
datascience<-function(data, response, covariates,rounding,family,parameter,intercept=T){

  ########################################################### DATA PREPARATION
  #In this section, we write a function that can detect character and factor variables, and turn them into dummy variables

  dummyize<-function(data){
    #Create a T/F vector testing if a vector is supposed to be a dummy
    factors_present<-sapply(data, class) == "character" | sapply(data, class) == "factor"
    factor_covariates<-(names(data)[factors_present])

    #It will be useful to store a list of non-factor covariates
    factors_not_present<-!factors_present
    non_factor_covariates<-(names(data)[factors_not_present])


    #We'll want information about the factors and their levels later to make the design matrix, so let's store them in a list
    factor_details<-rep(list(NA),length(factor_covariates))%>%
      set_names(factor_covariates)

    #If we detect factors, add a dummy variable for each level of each one. Then, drop the original variable.
    if(sum(factors_present>0)){ #If we detect factors...
      for(factor in factor_covariates){ #then for each factor...
        levels<-unique(data[,factor])%>%as.character #extract the unique levels
        factor_details[[factor]]<-levels #And store those levels for later

        #Now create the dummy columns through a loop
        for(i in 1:length(levels)){
          data<-data%>%mutate(!!levels[i]:=ifelse((!!sym(factor))==levels[i],1,0))
        }
      }

      #As a final step, you'll want to convert any true/false cols to 0,1.
      boolean_cols<-sapply(data, class)=="logical"
      data<-data%>%
        mutate_if(.predicate = is.logical,.funs = as.numeric)

      #Finally, we'll return the new, transformed dataset, all the info about the factors, and a list of non-factor covariates. This will be useful information later in the overall function.
      return(list(data%>%select(-!!(factor_covariates)),
                  factor_details,
                  non_factor_covariates)%>%set_names(c("dummy_data","factor_details","non_factor_covariates")))
    }else{
      return(data)
    }
  }

  #First, subset the data to only include the columns that the user asked for.
  data<-data%>%select(response,covariates)

  #We apply the dummyize function above if we detect any factor columns.
  if(class(dummyize(data))=="list"){  #Test if dataset has factor covariates. If it does:

    dum_data<-dummyize(data)$dummy_data #Get the new dummyized dataset

    #Now we'll sort the factor levels. Why? Because then we can easily create the reference level. This is somthing that I find to be a little troublesome about other packages: it is not super clear what the reference level is.
    sorted_factor_details<-lapply(X = dummyize(data)$factor_details,FUN = function(x){sort(x = x,decreasing = F)})

    #Make a list of the levels that will be collapsed into the intercept
    if(intercept==T){
    reference_vars<-sapply(X = sorted_factor_details,FUN = function(x){x[1]})
    }

    #We want to keep the newly created dummy variables, as well as the non-factor covariates selected by the user. There looks like a lot of objects here, but having these will make it easy to transform the data.
    if(intercept==T){
    vars_to_keep<-sapply(X = sorted_factor_details,FUN = function(x){x[-1]})%>%unlist%>%as.vector
    }else{
      vars_to_keep<-sapply(X = sorted_factor_details,FUN = function(x){x})%>%unlist%>%as.vector
    }
    non_factor_covs<-dummyize(data)$non_factor_covariates
    user_selected_non_factor_covs<-covariates[covariates %in% non_factor_covs]
    final_vars<-c(user_selected_non_factor_covs,vars_to_keep)

    #Finally, we use the previous blocks classification of what is what to produce the design matrix and response for analysis.
    if(intercept==T){
    X<-dum_data%>%select(c(user_selected_non_factor_covs,vars_to_keep))%>%as.matrix%>%cbind(1,.) #Add intercept
    }else{
      X<-dum_data%>%select(c(user_selected_non_factor_covs,vars_to_keep))%>%as.matrix
    }
    Y<-dum_data%>%select(response)%>%as.matrix

  }else{#If we didn't have to create dummy variables
    final_vars<-covariates
    if(intercept==T){
    X<-data%>%select(covariates)%>%as.matrix%>%cbind(1,.) #Add intercept
    }else{
      X<-data%>%select(covariates)%>%as.matrix
    }
    Y<-data%>%select(response)%>%as.matrix
  }

  ########################################################### LOESS (under development)

  if(family=="LOESS"){
    #Set up a place to store information for fitting the curve.

    Data<-cbind(Y,X)%>%data.frame
    local_constant_bandwith <- parameter
    curve_data<-seq(from=min(Data$x),
                    to=max(Data$x), by=.01)%>%
      data.frame()%>%
      setNames(c("Support"))%>%
      mutate(Lower=Support-(local_constant_bandwith),
             Upper=Support+(local_constant_bandwith))


    #Calcuting Density
    for(i in 1:nrow(curve_data)){
      points_in_window<-Data%>%filter(between(x,left = curve_data$Lower[i],
                                              right = curve_data$Upper[i]))

      gx_numerator_value=0
      gx_denominator_value=0
      for(j in 1:nrow(points_in_window)){
        u<-abs(points_in_window$x[j]-curve_data$Support[i])/local_constant_bandwith
        Epanechnikov_Kernel_j<-.75*(1-u^2)
        num_j <- Epanechnikov_Kernel_j*points_in_window$y[j]
        denom_j <- Epanechnikov_Kernel_j

        gx_denominator_value<-sum(gx_denominator_value,denom_j)
        gx_numerator_value<-sum(gx_numerator_value,num_j)
      }
      curve_data$gx_local_constant[i]<-gx_numerator_value/gx_denominator_value

    }

    return(ggplot()+
             geom_point(data = Data, aes(x=x,y=y))+
             geom_line(data = curve_data, aes(x=Support,y=gx_local_constant),color="red")+
             theme_light())

    report_key<-"LOESS"

  }

  ########################################################### RANDOM FOREST MODULE (under development)

  if(family=="rf"){

    rss_tracker<-data.frame(matrix(data = NA,nrow = nrow(data),ncol = (ncol(X)-1)))%>%
      set_names(final_vars)

    max_box_pop<-user_box_max+1 #Just artificially set box pop above users to keep while loop busy

    split_count<-1 #Keep track of how many splits you have

    while(max_box_pop>user_box_max){
      for(var in final_vars){
        for(i in 1:nrow(data)){


          split_name<-paste0("r_",split_count)

          data<-data%>%
            mutate((!!split_name):=ifelse((!!sym(var)) < data[i,var],1,2))

          r_means<-data%>%group_by(!!sym(split_name))%>%
            summarise(Y_mean=mean(!!sym(response)))

          #LEFT OF HERE
          r_mean_1<-ifelse(1 %in% as.matrix(r_means[,split_name]),r_means%>%filter((!!sym(split_name))==1)%>%select(Y_mean),0)[[1]]
          r_mean_2<-ifelse(2 %in% as.matrix(r_means[,split_name]),r_means%>%filter((!!sym(split_name))==2)%>%select(Y_mean),0)[[1]]

          r_1<-ifelse(r_mean_1==0,0,sum((data%>%filter((!!sym(split_name))==1)%>%select(!!sym(response))-r_mean_1)^2))
          r_2<-ifelse(r_mean_2==0,0,sum((data%>%filter((!!sym(split_name))==2)%>%select(!!sym(response))-r_mean_2)^2))

          rss_tracker[i,var]<-sum(r_1,r_2)

        }# End data
      }#End vars

      #Extract that best location, and then add a column to the df showing membership of a given split.
      best_split_loc<-which(rss_tracker == min(rss_tracker, na.rm = TRUE), arr.ind = TRUE)
      best_split_col<-names(rss_tracker)[best_split_loc[1,2]]
      best_split_row<-best_split_loc[1,1]


      #Now, save that as the split
      data<-data%>%
        mutate((!!split_name):=ifelse((!!sym(best_split_col)) < data[best_split_row,best_split_col],1,2))%>%
        unite(col = "box",starts_with(match = "r_"),remove = F) #This could be dangerous. Outside chance a dataset has this natively.

      boxes_summary<-data%>%
        group_by(box)%>%
        summarize(count=n())

      boxes<-boxes_summary$box

      max_box_pop<-max(boxes_summary$count)

      split_count<-split_count+1 #Iterate the split count up by one

    }#End boxes
  }#End random forest

  ########################################################### KNN MODULE
  if(family=="knn"){

    #Our goal is to make n x n matrix that shows the distance between all of the points.

    X<-X[,-1] #Remove intercept

    dist_mat<-matrix(data = NA,nrow = nrow(X),ncol = nrow(X))

    #Now, go through and calculate all the distances in the upper diagonal. I could probably vectorize this down the road.
    for(i in 1:nrow(dist_mat)){
      for(j in 1:ncol(dist_mat)){
        dist_mat[i,j]<-sqrt(sum((X[i,]-X[j,])^2)) #Euclidean distance.
      }
    }

    #Now, get the k-nearest neighbor positions
    dist_mat<-cbind(1:nrow(X),dist_mat)
    dist_list<-rep(list(NULL),nrow(X))
    for(i in 2:ncol(dist_mat-1)){
      k_nearest_points<-dist_mat[,c(1,i)]%>%data.frame%>%set_names(c("point","dist"))%>%
        arrange(dist)%>%slice(1:parameter)%>%select(point)
      dist_list[[i-1]]<-k_nearest_points
    }

    #For each data point, the estimate will take the mean of the response of those k-nearest neighbors
    fitted_vals<-rep(NA,nrow(X))
    for(i in 1:length(dist_list)){
      positions<-dist_list[[i]]%>%select(point)%>%unlist
      fitted_vals[i]<-mean(Y[positions,])
    }

    report_key<-"knn" #Administrative information for output
    return(data.frame(test_MSE=sqrt(mean((Y-fitted_vals)^2))))

  }#End KNN module

  ########################################################### GAUSSIAN REGRESSION MODULE (OLS)

  if(family=="gaussian"){

    #Get the beta estimates
    beta_hats<-(solve((t(X)%*%X)))%*%(t(X)%*%Y)

    #Get the fitted values and residuals
    preds<-X%*%beta_hats
    residuals<-Y-X%*%beta_hats

    #Get MSE
    MSE<-sum((residuals^2))/(nrow(X)-length(covariates)-1)

    #Get beta SEs
    beta_ses<-(MSE*solve(t(X) %*% X))%>%as.matrix%>%diag%>%as.numeric%>%sqrt

    #Get beta t values
    beta_ts<-beta_hats/beta_ses

    # Do statistical tests
    p_vals<-2*pt(q = abs(beta_ts),df =(nrow(X)-length(covariates)-1),lower.tail = FALSE)
    dist_of_test_stat<-"t"

    report_key<-"gaussian" #Administrative information for output
  }
  ########################################################### GLM REGRESSION (MLE)

  if(family %in% c("binomial","poisson")){

    #Make initial "guess" for parameter estimates
    beta_hats<-rep(0,ncol(X))

    #Threshhold. If beta estimates change by less than this, stop the algorithm.
    threshold<-.0000000001


    iteration_counter<-0 #Start at 1
    max_iterations<-1000 #Prevent an endless loop. Set the maximum number of iterations.

    #Make sure you allow while loop starting condition to be met
    improvement<-threshold+1 #Generically add one

    #Now do Newton-Raphson (a numerical method for maximum likelihood since MLE does not have closed form):
    while(iteration_counter <= max_iterations & improvement>=threshold){

      #Get weight matrix (depends on family)

      if(family=="binomial"){
        fits<-(exp(X%*%beta_hats)/(1+exp(X%*%beta_hats)))%>%as.vector
        W<-diag(fits*(1-fits))}

      if(family=="poisson"){
        fits<-exp(X%*%beta_hats)%>%as.vector
        W<-diag(fits)}

      #Update and print iterator
      print(paste0("Iteration ",iteration_counter))
      iteration_counter<-iteration_counter+1

      #Check if the hessian is invertible
      H<-(t(X)%*%W%*%X) #Hessian
      is_invertible<-tryCatch(solve(H),error=function(x)"error")%>%class

      if(is_invertible=="matrix"){#....if it is invertible, proceed:
        beta_hats_change<-solve(H)%*%t(X)%*%(Y - fits)
        improvement<-sum(abs(beta_hats_change))

        #Update the betas
        beta_hats<-beta_hats+beta_hats_change

      }else{#...when not invertible, stop:
        #Create
        improvement<-0}

      #Get beta SEs
      beta_ses<-solve(H)%>%as.matrix%>%diag%>%as.numeric%>%sqrt

      #Get beta t valyes
      beta_ts<-beta_hats/beta_ses

      # Do statistical tests
      p_vals<-2*pnorm(q = abs(beta_ts),mean=0,sd = 1,lower.tail = FALSE)

      #Specify distribution for output
      dist_of_test_stat<-"z"

    }#End of while loop

    report_key<-"glm" #Administrative information for output

  }#End of non-gaussian family cases

  ########################################################### MODEL RESULTS OUTPUT MODULE

  #Gather up any inversion errors
  if(report_key %in% c("gaussian","glm")){
    error_count<-if(family!="gaussian"){
      str_count(c(is_invertible),pattern = "character")%>%sum
    }else{
      0 #Just set to zero to bypass first condition in next step
    }


    if(error_count>0){#If any errors, then we'll print them out:
      if(is_invertible!="matrix"){print("Hessian is either computationally or exactly singular (not invertible)")}
    }else{#Else, print out the results...

      if(intercept==T){
        names=c("Intercept",final_vars)
      }else{names=c(final_vars)}
      #Main regression output:
      data.frame(parameter=names,
                 beta_hats%>%as.numeric%>%round(5),
                 beta_ses%>%as.numeric%>%round(5),
                 beta_ts%>%as.numeric%>%round(5),
                 p_vals%>%as.numeric)%>%
        set_names("parameter","Estimate","Std. Error",dist_of_test_stat,paste0("Pr(>|",dist_of_test_stat,"|)"))%>%print

      #Information about iterations
      if(family!="gaussian"){
        print(paste("Newton-Raphson iterations: ",iteration_counter-1,"(",max_iterations," allowed)" ))}

      #If we had dummy variables, show the reference levels
      if(class(dummyize(data))=="list"&intercept==T){
        print("Reference Levels:")
          print(reference_vars)
        }

    }
  }#End of regression reporting module
} #End of datascience function
