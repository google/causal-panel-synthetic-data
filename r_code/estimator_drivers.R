
#TODO(alexdkellogg): check whether I should rename id_var to be "entry" in here,
#    so I can just call entry in the main, rather than as.name(id_var)
estimate_gsynth_tot<- function(data_full, placebo_adjust=F, 
                                   bias_stat="mean", bias_horizon=5,
                                   id_var = "entry", 
                                   time_var = "period", 
                                   treat_indicator = "treatperiod_0", 
                                   outcome_var = "target", x_in = NULL,
                                   counterfac_var = "counter_factual", ...){
  
  gsynth_series=estimate_gsynth_series(data_full,
                               id_var = "entry", 
                               time_var = "period", 
                               treat_indicator = "treatperiod_0", 
                               outcome_var = "target", x_in = NULL,
                               counterfac_var = "counter_factual", ...)
  
  tot_series=compute_tot_se_jackknife(gsynth_series, stat_in=bias_stat,
                                      post_treat_only=F)
  
  
  if(placebo_adjust){
    #create a placebo using NN matching of treated to donors
    placebo_data=create_placebo_df(data_full)
    #estimate the method using this placebo set, where we know true TOT
    placebo_estimates=estimate_gsynth_series(placebo_data,
                                           id_var , 
                                           time_var, 
                                           treat_indicator, 
                                           outcome_var, x_in,
                                           counterfac_var,
                                           ...)
    #Compute the TOT per period, from which we recover the bias per period 
    placebo_tot=compute_tot_se_jackknife(placebo_estimates, stat_in=bias_stat,
                                         post_treat_only = F)
    debias_tib=placebo_tot %>%
      dplyr::mutate(abs_placebo_bias=jackknife_mean_abs_tot-mean_abs_cf_tot,
                    pct_placebo_bias=jackknife_mean_pct_tot-mean_pct_cf_tot) %>%
      dplyr::select(post_period_t,abs_placebo_bias,pct_placebo_bias )
    #Subtract the bias from the point effect for the relevant time periods.
    #Also, compute bias adjusted unit-level treatment effects
    #TODO(alexdkellogg): attach the original data ("est"), think about names
    #    for the output that's been debiased.
    tot_series=tot_series %>%
      dplyr::left_join(debias_tib, by="post_period_t") %>%
      dplyr::mutate(dplyr::across(.cols = tidyr::matches("jackknife.*abs|observed.*abs"),
                                  .fns = function(x){
                                    return(x-abs_placebo_bias)
                                  },
                                  #.names = "debiased_{col}"
                                  .names = "{col}"),
                    dplyr::across(.cols = tidyr::matches("jackknife.*pct|observed.*pct"),
                                  .fns = function(x){
                                    return(x-pct_placebo_bias)
                                  },
                                  .names = "{col}"))
  }
  return(tot_series)
  
} 