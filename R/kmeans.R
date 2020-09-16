#' @export kmeanscluster
#Manual settings
# data<-iris
# clusters<-3
# max_iterations=100
# print_iterations=T
# standardize=T
kmeanscluster<-function(data,clusters,max_iterations,print_iterations,standardize){

  data_original<-data
  ind<-sapply(X = data,is.numeric)

if(standardize==T){
data[ind]<-lapply(data[ind], scale)
}

# Note we have to cover the scenario where the number of clusters specified by the user does not cleanly divide the number of observations.
if(nrow(data)%%clusters==0){tags=c(rep(1:clusters,floor(nrow(data)/clusters)))
}else{tags=c(rep(1:clusters,floor(nrow(data)/clusters)),c(1:3)[1:(nrow(data)%%3)])}

#Handle the data label ids and cluster ids outside of the dataset, since we are doing operations on the main data and these columns might otherwise accidently get run through.
cluster_tracker<-data.frame(cluster=as.character(sample(x=tags,
                                                        size = nrow(data),
                                                        replace = F)),
                            id=as.character(1:nrow(data)))

#Now, for each group, we will calculate the centroid
all_centroids<-cbind(cluster_tracker,data)%>% #Temporarily read the ids/data together for centroid calcs
  group_by(cluster)%>%
  summarise_if(.predicate = is.numeric,.funs = mean)%>%
  arrange(cluster)%>%#Key, or else the for loop in the next block will get things wrong
  select(-cluster)

old_cluster_assignments<-as.character(cluster_tracker$cluster)
continue_while_loop<-T
reclassification_iteration<-0

while(continue_while_loop==T & (max_iterations>reclassification_iteration)==T){

# We should find the distance of each point from each centroid
for(i in 1:clusters){

  #The sweep function is basically like broadcasting in python, with the default operations being subtraction. Not sure how you would do other operations, but subtraction is all we need!
  data_minus_centroid<-as.matrix(sweep(x = data[ind],
                                       MARGIN = 2,
                                       STATS = as.matrix(all_centroids[i,])))

  cluster_tracker[,paste0("dist_to_clust",i,"_cent")]<-sqrt(rowSums(data_minus_centroid^2))

}

#Now for more tidyverse, unfortunately! Find the nearest centroid, then update the cluster of a given point
cluster_tracker_updated<-cluster_tracker%>%
  pivot_longer(cols = starts_with("dist_to_clust"),names_to = "centroid",values_to = "dist")%>%
  group_by(id,cluster)%>%
  filter(dist==min(dist))%>%
  ungroup%>%
  mutate(cluster=substring(text = centroid,first = 14,last = 14))%>%
  mutate(id=as.character(id))

#Update the iteration
reclassification_iteration<-reclassification_iteration+1

#Let the user know they have reached the maximum iterations. While loops only break off at the start of each loop, so doing it after updating reclassification_iteration is no problem.


#Print status update on iterations
if(print_iterations==T){
  print(paste0("Reclassification iteraction: ",reclassification_iteration))}

if((reclassification_iteration)==max_iterations){
  message(paste0("The maximum number of iterations you allowed, ",max_iterations," , has been met."))
}

#Now, lets test to see if the cluster assignment has changed or not
if(all.equal(cluster_tracker_updated$cluster,old_cluster_assignments)==T){
  continue_while_loop=F
}else{
  old_cluster_assignments<-cluster_tracker_updated$cluster
}

}#End of while loop

clustered_data<-inner_join(x = data_original%>%mutate(id=as.character(1:nrow(data))),
                           y = cluster_tracker_updated,
                           by="id")%>%
  select(-id,-dist,-centroid)
return(clustered_data)

}#End of function
