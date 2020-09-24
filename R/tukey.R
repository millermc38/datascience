#' @export tukey
tukey<-function(model_results,conf_level){
  all_contrasts<-expand.grid(results$parameter_estimates$Parameter,
                             results$parameter_estimates$Parameter)%>%
    filter(Var1!=Var2)%>%
    left_join(x = .,y=results$parameter_estimates[,c("Parameter","Estimate")],by=c("Var1"="Parameter"))%>%
    left_join(x = .,y=results$parameter_estimates[,c("Parameter","Estimate")],by=c("Var2"="Parameter"))%>%
    set_names("Var1","Var2","Var1_val","Var2_val")%>%
    mutate(difference=Var1_val-Var2_val)%>%
    as_tibble%>%
    #WARNING: Lots of temporary coding here that I will update in the future. Use of distinct is bound to fail at some point, and samples size (50) is not set dynamically and does not generalize from iris
    distinct(abs(difference),.keep_all = T)%>%
    select(Var1,Var2,difference)%>%
    mutate(lower=difference-(qtukey(p = conf_level,df = 150-3,nmeans =3 )/sqrt(2))*sqrt(2*(results$performance$MSE/50)),
           upper=difference+(qtukey(p = conf_level,df = 150-3,nmeans =3 )/sqrt(2))*sqrt(2*(results$performance$MSE/50)))%>%
    mutate(Contrast=paste0(Var1,"-",Var2))%>%
    select(Contrast,difference,lower,upper,-Var1,-Var2)
  return(all_contrasts)
}
