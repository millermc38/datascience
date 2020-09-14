#' @export cart

cart<-function(data, response, covariates,family,min_pre_partition_size,min_terminal_node_size,print_progress=F){
if(family=="tree"){ #test

  #Keep track of best splits
  best_split_tracker<-data.frame()
  report_key<-"tree"

  #First, we need to assign all observations to a "terminal_node". They are all part of the same terminal_node at the start.
  rf_data<-data%>%
    mutate(point_id=row_number())%>%
    mutate(terminal_node=1)#They all start in the same terminal_node

  #Just artificially set terminal_node pop above users to keep while loop busy
  max_box_pop<-min_terminal_node_size+1

  #Keep track of how many splits you have
  split_count<-1

  while(max_box_pop>min_terminal_node_size){

    #We need to track the RSS for any combination of var, terminal_node, data point split for any given r_i, or split i
    rss_tracker<-expand.grid(unique(rf_data$terminal_node),covariates,seq(1:nrow(rf_data)))%>%
      mutate(rss=NA)%>% #NA does not seem to work, but if I don't do this, then some invalid terminal_node/cov/point combos will be left at zero. Cannot have that!
      set_names(c("terminal_node","covariate","point","rss"))

    for(box_i in unique(rf_data$terminal_node)){

      #First, let's just filter down to the terminal_node/var/point combo of interest. We do this in seperate steps, because this eliminates searching that subsequent steps have to do.
      rf_data_filtered_step1<-rf_data%>%
        filter(terminal_node==box_i)

      if(nrow(rf_data_filtered_step1)<min_pre_partition_size){
        next
      }

      for(var in final_vars){
        for(i in rf_data_filtered_step1$point_id){

          split_name<-paste0("r_",split_count)

          if(print_progress==T){
            print(paste("Split",split_name,"terminal_node",box_i,"Covariate",var,"Obs",i,sep = ": "))}

          #Assign regions to a new column: r_i
          rf_data_filtered_step2<-rf_data_filtered_step1%>%
            mutate((!!split_name):=ifelse((!!sym(var)) < rf_data_filtered_step1[rf_data_filtered_step1$point_id==i,var],1,2))

          #Now, we need to at least have an r1 and and r2, so if there are not at least two points in each region, we need to move on to the next point in the loop. We also need to calculate the mean for each group in case we want to get the rss

          r_means<-rf_data_filtered_step2%>%
            group_by(!!sym(split_name))%>%
            summarise(Y_mean=mean(!!sym(response)),r_counts=n())%>%
            mutate(rss=NA)%>% #Store RSS for each r.
            data.frame


          if((min(r_means$r_counts)<min_terminal_node_size)==T|nrow(r_means)<2){
            next
          }else{
            #Calculate rss for each region
            for(r in 1:nrow(r_means)){
              #If the r_i appears the r_means, filter the frame down to that mean, otherwise, select 0.
              r_mean<-ifelse(r %in% as.matrix(r_means[,split_name]), #Check if r is even a valid region (if there is just one region, for example)
                             r_means%>%filter((!!sym(split_name))==r)%>%select(Y_mean)#...then, get that mean
                             ,0)[[1]]#...else, just set rss to zero

              #Now, find the rss if it was not zero!
              r_means[r,"rss"]<-ifelse(r_mean==0,0,sum((rf_data_filtered_step2%>%
                                                          filter((!!sym(split_name))==r)%>%
                                                          select(!!sym(response))-r_mean)^2))
              #For some reason this loop overwrites r_means[1,2] over with a zero...why?
            }

            rss_tracker[rss_tracker$terminal_node==box_i&
                          rss_tracker$point==i&
                          rss_tracker$covariate==var,"rss"]<-sum(r_means$rss)
          }#End if else condition for min amount in a region
        }# End data
      }#End vars
    }#End terminal_node


    #It is possible that all of the boxes returned NA even if boxes have a handful of observations. This occurs when the terminal_node can be split, but doing so would cause a terminal_node to have fewer observations than the user specified.

    if(mean(is.na(rss_tracker$rss))==1){

      print("while loop ending: no more splits")
      rf_fitted<-left_join(x = rf_data,y = boxes_summary,by="terminal_node")%>%
        select(!!sym(response),mean,everything())
      return(list(rf_fitted, best_split_tracker))
      #left off here...why is it not breaking here and heading to "return"?
      max_box_pop<-0 #break while loop

    }else{

      #print(rss_tracker)
      #Extract that best location
      best_split<-rss_tracker%>%filter(rss==min(rss,na.rm = T))%>%
        slice(1) #Might be multiple best splits, just pick the top one

      best_split_tracker<-rbind(best_split_tracker,best_split)


      #Mark true for r1
      rf_condition_for_r1<-rf_data$terminal_node==best_split$terminal_node&
        rf_data[,as.character(best_split$covariate)]<rf_data[best_split$point,as.character(best_split$covariate)]

      #Mark true for r2
      rf_condition_for_r2<-rf_data$terminal_node==best_split$terminal_node&
        rf_data[,as.character(best_split$covariate)]>=rf_data[best_split$point,as.character(best_split$covariate)]

      #Mark true for neither
      rf_condition_for_none<-rf_data$terminal_node!=best_split$terminal_node

      #Now throw that into the dataframe
      rf_data<-rf_data%>%
        mutate((!!split_name):=ifelse(rf_condition_for_r1==T,1,ifelse(rf_condition_for_r2==T,2,NA)))%>%
        unite(col = "terminal_node",starts_with(match = "r_"),remove = F) #This could be dangerous. Outside chance a dataset has this natively.

      boxes_summary<-rf_data%>%
        group_by(terminal_node)%>%
        summarize(count=n(),mean=mean(!!sym(response)))

      max_box_pop<-max(boxes_summary$count)

      split_count<-split_count+1 #Iterate the split count up by one

    }#End condition where we check if a split is possible

  }#End boxes


  #Now merge the fitted values on the data
  rf_fitted<-left_join(x = rf_data,y = boxes_summary,by="terminal_node")%>%
    select(!!sym(response),mean,everything())

  #Merge the cutpoints onto the best_split_tracker
  best_split_tracker<-left_join(x = best_split_tracker,y = rf_data,by=c("point"="point_id"))

  return(list(rf_fitted, best_split_tracker))

  #Think it would be bestter to view it more like an actual tree. Start with a main dataset, filter down the the relevant half whenever it forks, save the rule, and keep splitting. It can split in parallel that way.
}#End random forest
}#End functino

