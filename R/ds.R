#' @export
datascience<-function(data, response, covariates,rounding,family,parameter){
  ########################################################### SETUP

  #First, we need to write a function that can handle dummy variables
  dummyize<-function(data){
    #Create a T/F vector testing if a vector is supposed to be a dummy
    factors_present<-sapply(data, class) == "character" | sapply(data, class) == "factor"
    factor_covariates<-(names(data)[factors_present])

    #It will be useful to have a list of non factor covariates as well in your output
    factors_not_present<-!factors_present
    non_factor_covariates<-(names(data)[factors_not_present])


    #We'll want information about the factors and their levels later to make the design matrix, so let's store them in a list
    factor_details<-rep(list(NA),length(factor_covariates))%>%
      set_names(factor_covariates)


    if(sum(factors_present>0)){
      #Remove it from the dataframe and create a dummy variable for each one
      for(factor in factor_covariates){
        levels<-unique(data[,factor])%>%as.character

        factor_details[[factor]]<-levels

        levels_len<-length(levels)
        for(i in 1:levels_len){
          data<-data%>%mutate(!!levels[i]:=ifelse((!!sym(factor))==levels[i],1,0))
        }
      }
      #As a final step, you'll want to convert any true false cols to 0,1.
      boolean_cols<-sapply(data, class)=="logical"
      data<-data%>%
        mutate_if(.predicate = is.logical,.funs = as.numeric)

      return(list(data%>%select(-!!(factor_covariates)),
                  factor_details,
                  non_factor_covariates)%>%set_names(c("dummy_data","factor_details","non_factor_covariates")))
    }else{
      return(data)
    }
  }

  #First, subset the data to only include inforamtion they asked for
  data<-data%>%select(response,covariates)




  #MUTATE THE DATAFRAME AS NEED!!!!!!!!!!!!!!!!!!!!!!!!!
  if(class(dummyize(data))=="list"){  #if we had to create dummy variables
    dum_data<-dummyize(data)$dummy_data
    sorted_factor_details<-lapply(X = dummyize(data)$factor_details,FUN = function(x){sort(x = x,decreasing = F)})

    #Make a list of the variables that will be collapsed into the intercept
    reference_vars<-sapply(X = sorted_factor_details,FUN = function(x){x[1]})

    #We want to keep the newlycreated dummy variables, as well as the non factor covariates selected by the user.
    vars_to_keep<-sapply(X = sorted_factor_details,FUN = function(x){x[-1]})%>%unlist%>%as.vector
    non_factor_covs<-dummyize(data)$non_factor_covariates
    user_selected_non_factor_covs<-covariates[covariates %in% non_factor_covs]
    final_vars<-c(user_selected_non_factor_covs,vars_to_keep)

    X<-dum_data%>%select(c(user_selected_non_factor_covs,vars_to_keep))%>%as.matrix%>%cbind(1,.) #Add intercept
    Y<-dum_data%>%select(response)%>%as.matrix
  }else{#if we didn't have to create dummy variables
    final_vars<-covariates
    X<-data%>%select(covariates)%>%as.matrix%>%cbind(1,.) #Add intercept
    Y<-data%>%select(response)%>%as.matrix
  }



  ########################################################### RF
  if(family=="rf"){


    rss_tracker<-data.frame(matrix(data = NA,nrow = nrow(data),ncol = (ncol(X)-1)))%>%
      set_names(final_vars)

    for(var in final_vars){
      for(i in 1:nrow(data)){

        split_name<-paste0("r_","r")

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

      }
    }

    #Extract that best location, and then add a column to the df showing membership of a given split.
    best_split_loc<-which(rss_tracker == min(rss_tracker, na.rm = TRUE), arr.ind = TRUE)
    best_split_col<-names(rss_tracker)[best_split_loc[1,2]]
    best_split_row<-best_split_loc[1,1]


    #Now, save that as the split
    data<-data%>%
      mutate((!!split_name):=ifelse((!!sym(best_split_col)) < data[best_split_row,best_split_col],1,2))


  }

  ########################################################### KNN
  if(family=="knn"){
    #First, make nxn matrix that shows the distance between all of the points

    X<-X[,-1] #No point in having an intercept in KNN

    dist_mat<-matrix(data = NA,nrow = nrow(X),ncol = nrow(X))

    #Now, go through and calculate all the distances in the upper diag. I could probably vectorize this
    for(i in 1:nrow(dist_mat)){
      for(j in 1:ncol(dist_mat)){
        dist_mat[i,j]<-sqrt(sum((X[i,]-X[j,])^2))
      }
    }

    #Don't want to consider distance from self
    # for(i in 1:nrow(dist_mat)){
    #   dist_mat[i,i]<-NA
    # }


    #now, get the k nearest neighbor positions
    dist_mat<-cbind(1:nrow(X),dist_mat)
    dist_list<-rep(list(NULL),nrow(X))
    for(i in 2:ncol(dist_mat-1)){

      k_nearest_points<-dist_mat[,c(1,i)]%>%data.frame%>%set_names(c("point","dist"))%>%
        arrange(dist)%>%slice(1:parameter)%>%select(point)

      dist_list[[i-1]]<-k_nearest_points
    }

    fitted_vals<-rep(NA,nrow(X))
    #For each one, take the mean of the response of those nearest neighbors
    for(i in 1:length(dist_list)){
      positions<-dist_list[[i]]%>%select(point)%>%unlist
      fitted_vals[i]<-mean(Y[positions,])
    }

    report_key<-"knn"
    return(
      data.frame(test_MSE=sqrt(mean((Y-fitted_vals)^2))))

  }#End module




  ########################################################### GAUSSIAN (OLS)
  if(family=="gaussian"){
    #Get the beta estimates
    beta_hats<-(solve((t(X)%*%X)))%*%(t(X)%*%Y) #WRITE ABOUT THEORY

    #Get the fitted values and residuals
    preds<-X%*%beta_hats
    residuals<-Y-X%*%beta_hats

    #Get MSE
    MSE<-sum((residuals^2))/(nrow(X)-length(covariates)-1)

    #Get beta ses
    beta_ses<-(MSE*solve(t(X) %*% X))%>%as.matrix%>%diag%>%as.numeric%>%sqrt

    #Get beta t vals
    beta_ts<-beta_hats/beta_ses

    # Do statistical tests
    p_vals<-2*pt(q = abs(beta_ts),df =(nrow(X)-length(covariates)-1),lower.tail = FALSE)
    dist_of_test_stat<-"t"

    report_key<-"gaussian"
  }
  ########################################################### GLM REGRESSION (MLE)
  if(family %in% c("binomial","poisson")){
    #Make initial "guess" for parameter estimates
    beta_hats<-rep(0,ncol(X))

    #Threshhold. If beta estimates change by less than this, stop the algorithm
    threshold<-.0000000001

    #Prevent and endless loop. Set the maximum number of iterations
    iteration_counter<-0 #Start at 1
    max_iterations<-1000

    #Make sure you allow while loop starting condition to be met
    improvement<-threshold+1 #Generically add one

    #Now do Newton-Raphson (a numerical method for maximum liklihood since not closed form):
    while(iteration_counter <= max_iterations & improvement>=threshold){


      #Get weight matrix (depends on family)

      if(family=="binomial"){
        fits<-(exp(X%*%beta_hats)/(1+exp(X%*%beta_hats)))%>%as.vector
        W<-diag(fits*(1-fits))}

      if(family=="poisson"){
        fits<-exp(X%*%beta_hats)%>%as.vector
        W<-diag(fits)}

      #Get change in betas from this
      print(paste0("Iteration ",iteration_counter))
      iteration_counter<-iteration_counter+1

      #Check if the hessian is invertible
      H<-(t(X)%*%W%*%X) #Hessian
      is_invertible<-tryCatch(solve(H),error=function(x)"error")%>%class

      if(is_invertible=="matrix"){#....if it is, proceed
        beta_hats_change<-solve(H)%*%t(X)%*%(Y - fits)
        improvement<-sum(abs(beta_hats_change))

        #Update the betas
        beta_hats<-beta_hats+beta_hats_change

      }else{#...when not invertible
        #Create
        improvement<-0}

      # beta_ses<-rep(NA,length(beta_hats))
      # beta_ts<-rep(NA,length(beta_hats))
      #p_vals<-rep(NA,length(beta_hats))

      #Get beta ses
      beta_ses<-solve(H)%>%as.matrix%>%diag%>%as.numeric%>%sqrt

      #Get beta t vals
      beta_ts<-beta_hats/beta_ses

      # Do statistical tests
      p_vals<-2*pnorm(q = abs(beta_ts),mean=0,sd = 1,lower.tail = FALSE)

      #Specify distribution
      dist_of_test_stat<-"z"
    }#End of while loop

    report_key<-"glm"
  }#End of non-gaussian family cases

  ###########################################################OUTPUT

  if(report_key %in% c("gaussian","glm")){
    error_count<-if(family!="gaussian"){
      str_count(c(is_invertible),pattern = "character")%>%sum
    }else{
      0 #Just set to zero to bypass first condition in next step
    }


    if(error_count>0){#....print errors
      if(is_invertible!="matrix"){print("Hessian is either computationally or exactly singular (not invertible)")}
    }else{#...print results as usual
      data.frame(parameter=c("Intercept",final_vars),
                 beta_hats%>%as.numeric%>%round(5),
                 beta_ses%>%as.numeric%>%round(5),
                 beta_ts%>%as.numeric%>%round(5),
                 p_vals%>%as.numeric)%>%
        set_names("parameter","Estimate","Std. Error",dist_of_test_stat,paste0("Pr(>|",dist_of_test_stat,"|)"))%>%print
      if(family!="gaussian"){
        print(paste("Newton-Raphson iterations: ",iteration_counter-1,"(",max_iterations," allowed)" ))}

      if(class(dummyize(data))=="list"){
        print("Reference Levels:")
        print(reference_vars)}

    }
    #What is F-test in regression?
  }#End of regression reporting module
}
