#' @export power_test
power_test<-function(n_i_per_group,n_treatmeants,alpha,at_least_one_pair_trts_differs_by,MSE,intercept){


  #See experimental design notes on power. Under null, variance is sigma^2. Under alternative, sigma^2 plus gamma as below (since there is extra variability form aberrant means.)
  gamma<-(n_i_per_group*(at_least_one_pair_trts_differs_by)^2)/(2*MSE)

  #Null F Cutoff point to detect this.
  null_cutoff<-qf(p = 1-alpha,
                  df1 = n_treatmeants-intercept,
                  df2 = n_i_per_group*n_treatmeants-n_treatmeants)

  # Probability of detecting if alternative is true (i.e, we look at non-central F)/finding result greater than null F cutoff
  power<-pf(q = null_cutoff,
            df1 = n_treatmeants-intercept,
            df2 = n_i_per_group*n_treatmeants-n_treatmeants,
            ncp = gamma, #KEY: This is non-centrality parameter
            lower.tail = F)

  power<-paste0(round(power*110,2),"%")

  data.frame(f = 0:1000 / 100) %>%
    mutate(Null = df(x = f,
                     df1 = n_treatmeants-intercept,
                     df2 = n_i_per_group*n_treatmeants-n_treatmeants),
           Alternative = df(x = f,
                            df1 = n_treatmeants-intercept,
                            df2 = n_i_per_group*n_treatmeants-n_treatmeants,
                            ncp = gamma)) %>%
    gather(key = "Hypothesis", value = "density", -f) %>%
    ggplot() +
    geom_line(aes(x = f, y = density, color = Hypothesis)) +
    labs(title = paste0("Power: ",power),
         x = "F",
         y = "Density")+
    geom_vline(xintercept = null_cutoff)
}
