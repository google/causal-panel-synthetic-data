pacman::p_load(dplyr, furrr, stringr, tidyr, tibble)
compute_cf_tot <- function(tib_inp, stat_in_inp,t_var_inp,
                           treat_t_inp, outcome_inp, cf_inp){
  #Helper function to our bootstrap and jackknife procedures, used to compute
  #the relevant counterfactual effects
  
  #Args
  #tib_inp: original data
  #stat_in_inp: string for the desired statistic -- mean, median, aggregate.
  #t_var_inp: string name of the time column.
  #treat_t_inp: string name of the column indicating treatment start period.
  #outcome_inp: string name of the outcome variable of interest.
  #cf_inp: string name of the counterfactual outcome
  
  #output: tibble with the counterfactual effect by time since treatment
  
  #column names for output
  cf_abs_col=paste(stat_in_inp, "_abs_cf_tot", sep="")
  cf_pct_col=paste(stat_in_inp, "_pct_cf_tot", sep="")
  
  #Computing the relevant effect and storing as tibble
  cf_tot_df <- tib_inp %>%
    dplyr::mutate(post_period_t = !!as.name(t_var_inp) -
                    !!as.name(treat_t_inp)) %>%
    dplyr::filter(post_period_t >= 0) %>%
    dplyr::group_by(post_period_t) %>%
    dplyr::summarise(
      !!as.name(cf_abs_col):=dplyr::case_when(
        stat_in_inp=="aggregate"~sum(!!as.name(outcome_inp) -
                                   !!as.name(cf_inp)),
        stat_in_inp=="mean"~mean(!!as.name(outcome_inp)- 
                               !!as.name(cf_inp)),
        stat_in_inp=="median"~median(!!as.name(outcome_inp) - 
                                   !!as.name(cf_inp))),
      !!as.name(cf_pct_col):=dplyr::case_when(
        stat_in_inp=="aggregate"~(sum(!!as.name(outcome_inp)) / 
                                sum(!!as.name(cf_inp))) - 1,
        stat_in_inp=="mean"~mean((!!as.name(outcome_inp) / 
                                !!as.name(cf_inp)) - 1),
        stat_in_inp=="median"~median((!!as.name(outcome_inp) / 
                                    !!as.name(cf_inp)) - 1))) %>%
    dplyr::ungroup()
  
  return(cf_tot_df)
}



compute_tot_se_jackknife <- function(estimated_series_df, time_var = "period", treat_period_var = "Treatment_Period",
                                     pred_var = "point.pred", outcome_var = "response",
                                     stat_in = "mean", alpha_ci = 0.95,
                                     compute_cf_eff = T, counterfac_var = "counter_factual") {
  # Computes jackknife estimates (from resample package) of the treatment effect on the treated, by period
  
  # Args
  # estimated_series_df: long-form dataframe with the estimated effects, aka output from one of the methods above.
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_period_var:column name of the period for which the unit's treatment starts (assumes it ends at T max)
  # stat_in: string indicating the statistic which we hope to get a bootstrap mean and se for, either "mean", "aggregate", or "median" typical
  # alpha_ci: number between 0 and 1, indicating the confidence interval desires
  # compute_cf_eff: boolean flag for whether a column of counterfactual effects should be appended to output
  # counterfac_eff: string variable name of the counterfactual effect to estimate, typically "cf_point.effect or "cf_pct.effect" from above
  
  # Output
  # Tibble containing, by post_treat_period, the sample stat and bootstrapped mean of boot_var, as well as upper and lower bounds
  # that make up the alpha_ci*100% confidence interval
  
  #Generate a post treatment time indicator -- time since treatment.
  effect_series_postonly <- estimated_series_df %>%
    mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
    filter(post_period_t >= 0)
  
  #Compute the relevant effect using the counterfactual data
  if (compute_cf_eff) {
    cf_tot_df<-compute_cf_tot(tib_inp=estimated_series_df,
                              stat_in_inp=stat_in,t_var_inp=time_var,
                              treat_t_inp=treat_period_var,
                              outcome_inp=outcome_var, cf_inp=counterfac_var)
  }
  
  #Generate a list of observed effects by time-since-treatment
  effects_to_jk<-effect_series_postonly %>%
    dplyr::mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var),
                  pct_eff = (!!as.name(outcome_var) / !!as.name(pred_var)) - 1) %>%
    split(.[["post_period_t"]])
  
  jk_comp_abs<-effects_to_jk %>%
    furrr::future_map(~ .[["abs_eff"]])
  #Jackknife the relevant statistic for both absolute and percentage terms
  jk_comp_abs<- switch (stat_in,
                        "mean" =  jk_comp_abs %>% 
                          furrr::future_map(~ jackknife(., mean)),
                        "median"= jk_comp_abs %>% 
                          furrr::future_map(~ jackknife(., median)),
                        "aggregate"=jk_comp_abs %>% 
                          furrr::future_map(~ jackknife(., sum))
  )
  jk_comp_pct<-effects_to_jk %>%
    furrr::future_map(~ .[["pct_eff"]])
  #Jackknife the relevant statistic for both absolute and percentage terms
  jk_comp_pct<- switch (stat_in,
                        "mean" =  jk_comp_pct %>% 
                          furrr::future_map(~ jackknife(., mean)),
                        "median"= jk_comp_pct %>% 
                          furrr::future_map(~ jackknife(., median)),
                        "aggregate"=jk_comp_pct %>% 
                          furrr::future_map(~ jackknife(., function(x) {
                            (sum(x[[outcome_var]]) / sum(x[[pred_var]])) - 1
                          }))
  )
  #Gather the bounds on the confidence intervals and format as tibble
  jk_comp_ci <-format_jackknife(
    jk_est=jk_comp_abs,ci_in=alpha_ci,col_name= 
      paste(stat_in,"_abs_tot", sep="")) %>%
    dplyr::inner_join(
      format_jackknife(jk_est=jk_comp_pct,ci_in=alpha_ci, col_name= 
                         paste(stat_in,"_pct_tot",sep="")), by="post_period_t")
  #Add the number of treated units per period
  jk_ci_byT <- effect_series_postonly %>%
    dplyr::group_by(post_period_t) %>%
    dplyr::summarise(treated_n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(jk_comp_ci, by="post_period_t")
  #Add counterfactual if relevant
  if(compute_cf_eff){
    jk_ci_byT=jk_ci_byT %>% dplyr::left_join(cf_tot_df,by = "post_period_t")
  }
  return(jk_ci_byT)
}





compute_ci_bounds_bootstrap <- function(estimated_series_df, time_var = "period", treat_period_var = "Treatment_Period",
                                        pred_var = "point.pred", outcome_var = "response",
                                        stat_in = "mean", alpha_ci = 0.95,
                                        compute_cf_eff = T, counterfac_var = "counter_factual",
                                        nboots = 10000) {
  # Computes bootstrap estimates (from resample package) of the treatment effect on the treated, by period
  
  # Args
  # estimated_series_df: long-form dataframe with the estimated effects, aka output from one of the methods above.
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_period_var:column name of the period for which the unit's treatment starts (assumes it ends at T max)
  # boot_var: string indicating the effect we wish to bootstrap, either "point.effect" or "pct.effect" works
  # stat_in: string indicating the statistic which we hope to get a bootstrap mean and se for, either "mean", "aggregate", "median" typical
  # alpha_ci: number between 0 and 1, indicating the confidence interval desires
  # compute_cf_eff: boolean flag for whether a column of counterfactual effects should be appended to output
  # counterfac_eff: string variable name of the counterfactual effect to estimate, typically "cf_point.effect or "cf_pct.effect" from above
  # nboots, number of bootstrap samples
  
  # Output
  # Tibble containing, by post_treat_period, the sample stat and bootstrapped mean of boot_var, as well as upper and lower bounds
  # that make up the alpha_ci*100% confidence interval
  
  #Generate a post treatment time indicator -- time since treatment.
  effect_series_postonly <- estimated_series_df %>%
    mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
    filter(post_period_t >= 0)
  
  #Compute the relevant effect using the counterfactual data
  if (compute_cf_eff) {
    cf_tot_df<-compute_cf_tot(tib_inp=estimated_series_df,
                              stat_in_inp=stat_in,t_var_inp=time_var,
                              treat_t_inp=treat_period_var,
                              outcome_inp=outcome_var, cf_inp=counterfac_var)
  }
  
  #Generate a list of observed effects by time-since-treatment
  effects_to_boot<-effect_series_postonly %>%
    dplyr::mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var),
                  pct_eff = (!!as.name(outcome_var) / !!as.name(pred_var)) - 1) %>%
    split(.[["post_period_t"]])
  
  boot_comp_abs<-effects_to_boot %>%
    furrr::future_map(~ .[["abs_eff"]])
  #Jackknife the relevant statistic for both absolute and percentage terms
  boot_comp_abs<- switch (stat_in,
                        "mean" =  boot_comp_abs %>% 
                          furrr::future_map(~ bootstrap(., mean, R = nboots)),
                        "median"= boot_comp_abs %>% 
                          furrr::future_map(~ bootstrap(., median, R = nboots)),
                        "aggregate"=boot_comp_abs %>% 
                          furrr::future_map(~ bootstrap(., sum, R = nboots))
  )
  boot_comp_pct<-effects_to_boot %>%
    furrr::future_map(~ .[["pct_eff"]])
  #Jackknife the relevant statistic for both absolute and percentage terms
  boot_comp_pct<- switch (stat_in,
                        "mean" =  boot_comp_pct %>% 
                          furrr::future_map(~ bootstrap(., mean, R = nboots)),
                        "median"= boot_comp_pct %>% 
                          furrr::future_map(~ bootstrap(., median, R = nboots)),
                        "aggregate"=boot_comp_pct %>% 
                          furrr::future_map(~ bootstrap(., function(x) {
                            (sum(x[[outcome_var]]) / sum(x[[pred_var]])) - 1
                          }, R = nboots))
  )
  #Gather the bounds on the confidence intervals and format as tibble
  boot_comp_ci <-format_jackknife(
    jk_est=boot_comp_abs,ci_in=alpha_ci,col_name= 
      paste(stat_in,"_abs_tot", sep="")) %>%
    dplyr::inner_join(
      format_jackknife(jk_est=boot_comp_pct,ci_in=alpha_ci, col_name= 
                         paste(stat_in,"_pct_tot",sep="")), 
      by="post_period_t") %>%
    rename_with( ~stringr::str_replace(.x,"jackknife", "bootstrap"),
                 tidyr::starts_with("jackknife"))
    
  #Add the number of treated units per period
  boot_ci_byT <- effect_series_postonly %>%
    dplyr::group_by(post_period_t) %>%
    dplyr::summarise(treated_n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(boot_comp_ci, by="post_period_t")
  #Add counterfactual if relevant
  if(compute_cf_eff){
    boot_ci_byT=boot_ci_byT %>% dplyr::inner_join(cf_tot_df,by = "post_period_t")
  }
  return(boot_ci_byT)
}
  

format_jackknife<-function(jk_est, ci_in, col_name){
  #Helper function to compute the confidence bands for the jackknife
  #Also reformats the output for convenient printing
  
  #Args
  
  
  statname<- paste("jackknife_", col_name, sep = "")
  obs_statname <- paste("observed_", col_name, sep = "")
  
  jk_est_ci<- 
    furrr::future_pmap(list(jk_est,
                            list(probs = 
                                   c(0.5 - ci_in / 2,0.5 + ci_in / 2))),
                       CI.t) %>%
    do.call(rbind, .) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(post_period_t = seq_len(nrow(.))-1,
                  !!as.name(statname) := jk_est %>% 
                    furrr::future_map(~ .[["stats"]]) %>% 
                    furrr::future_map(~ .[["Mean"]]) %>%
                    unlist(),
                  !!as.name(obs_statname) := jk_est %>% 
                    furrr::future_map(~ .[["stats"]]) %>% 
                    furrr::future_map(~ .[["Observed"]]) %>%
                    unlist()
    )
  
  names(jk_est_ci)[1:2] <- c(paste("jackknife_lb_",col_name, sep=""),
                             paste("jackknife_ub_", col_name, sep = ""))
  
  return(jk_est_ci)
}
