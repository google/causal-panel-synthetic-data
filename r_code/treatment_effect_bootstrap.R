library(CausalImpact)
library(Matrix)
library(tsfeatures)
library(tsibble)
library(ggfortify)
library(gsynth)
library(augsynth)
library(tidyr)
library(panelView)
library(synthdid)
library(resample)
library(mvtnorm)
library(janitor)
library("qpcR")
library(dplyr)
library(furrr)
library(gridExtra)
library(quadprog)
library(rstatix)
library(ForecastComb)


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


#TODO(alexdkellogg): keep cleaning this code to cut it down -- use format jackknife
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
                        "median"= jkcomp_abs %>% 
                          furrr::future_map(~ jackknife(., median)),
                        "aggregate"=jkcomp_abs %>% 
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
  
  
  statname_abs <- paste("jackknife_", stat_in, "_abs_tot", sep = "")
  statname_pct <- paste("jackknife_", stat_in, "_pct_tot", sep = "")
  obs_statname_abs <- paste("observed_", stat_in, "_abs_tot", sep = "")
  obs_statname_pct <- paste("observed_", stat_in, "_pct_tot", sep = "")
  jackknife_comp_ci_abs <- furrr::future_pmap(list(jk_comp_abs, list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = effect_series_postonly %>% 
             distinct(post_period_t) %>%
             pull()) %>%
    mutate(
      !!as.name(statname_abs) := jk_comp_abs %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist(),
      !!as.name(obs_statname_abs) := jk_comp_abs %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Observed"]]) %>%
        unlist()
    )
  names(jackknife_comp_ci_abs)[1:2] <- c("jackknife_lb_abs_tot", "jackknife_ub_abs_tot")
  
  jackknife_comp_ci_pct <- furrr::future_pmap(list(jk_comp_pct,
                                                   list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = effect_series_postonly %>%
             distinct(post_period_t) %>%
             pull()) %>%
    mutate(
      !!as.name(statname_pct) := jk_comp_pct %>%
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist(),
      !!as.name(obs_statname_pct) := jk_comp_pct %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Observed"]]) %>%
        unlist()
    )
  
  names(jackknife_comp_ci_pct)[1:2] <- c("jackknife_lb_pct_tot", "jackknife_ub_pct_tot")
  
  jackknife_comp_ci <- jackknife_comp_ci_abs %>% inner_join(jackknife_comp_ci_pct, by = "post_period_t")
  # Manual version of the above code
  # jackknife_mean_bias= jackknife_comp%>% furrr::future_map(~ .$"stats") %>% furrr::future_map(~ .$"Mean") %>%
  #   unlist() %>% tibble("de_biased_mean"=.,  "post_period_t"=effect_series_postonly %>% distinct(post_period_t) %>% pull())
  #
  # jackknife_se=jackknife_comp%>% furrr::future_map(~ .$"stats") %>% furrr::future_map(~ .$"SE") %>%
  #   unlist() %>% tibble("se"=.,"post_period_t"=effect_series_postonly %>% distinct(post_period_t) %>% pull())
  
  
  stat_tot_sample <- effect_series_postonly %>%
    group_by(post_period_t) %>%
    summarise(treated_n = n()) %>%
    ungroup()
  
  
  if (compute_cf_eff) {
    jackknife_ci_by_postperiod <- stat_tot_sample %>%
      left_join(
        jackknife_comp_ci,
        by = "post_period_t"
      ) %>%
      left_join(
        cf_tot_df,
        by = "post_period_t"
      )
  }
  if (!compute_cf_eff) {
    jackknife_ci_by_postperiod <- stat_tot_sample %>% 
      left_join(
        jackknife_comp_ci,
        by = "post_period_t"
      )
  }
  
  
  return(jackknife_ci_by_postperiod)
}




#TODO(alexdkellogg): keep cleaning this code to cut it down; use the CF helper.
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
  
  effect_series_postonly <- estimated_series_df %>%
    mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
    filter(post_period_t >= 0)
  
  if (compute_cf_eff) {
    cf_tot_df <- estimated_series_df %>%
      mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
      filter(post_period_t >= 0)
    if (tolower(stat_in) == "aggregate") {
      cf_tot_df <- cf_tot_df %>%
        group_by(post_period_t) %>%
        summarise(
          agg_abs_cf_tot = sum(!!as.name(outcome_var) - !!as.name(counterfac_var)),
          agg_pct_cf_tot = (sum(!!as.name(outcome_var)) / sum(!!as.name(counterfac_var))) - 1
        ) %>%
        ungroup()
    }
    if (tolower(stat_in) == "mean") {
      cf_tot_df <- cf_tot_df %>%
        group_by(post_period_t) %>%
        summarise(
          mean_cf_tot = mean(!!as.name(counterfac_var)),
          mean_pct_cf_tot = mean((!!as.name(outcome_var) / !!as.name(counterfac_var)) - 1)
        ) %>%
        ungroup()
    }
    if (tolower(stat_in) == "median") {
      cf_tot_df <- cf_tot_df %>%
        group_by(post_period_t) %>%
        summarise(
          median_abs_cf_tot = median(!!as.name(outcome_var) - !!as.name(counterfac_var)),
          median_pct_cf_tot = median((!!as.name(outcome_var) / !!as.name(counterfac_var)) - 1)
        ) %>%
        ungroup()
    }
  }
  
  # is the double future-map here the slow thing? How does the speed compare to lapply?
  if (tolower(stat_in) == "aggregate") {
    # Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ bootstrap(., mean, R = nboots))
    
    bootstrap_comp_pct <- effect_series_postonly %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ bootstrap(., function(x) {
        (sum(x[[outcome_var]]) / sum(x[[pred_var]])) - 1
      }, R = nboots))
  }
  if (tolower(stat_in) == "mean") {
    # Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ bootstrap(., mean, R = nboots))
    
    # Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp_pct <- effect_series_postonly %>%
      mutate(pct_eff = (!!as.name(outcome_var) / !!as.name(pred_var)) - 1) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["pct_eff"]]) %>%
      furrr::future_map(~ bootstrap(., mean, R = nboots))
  }
  
  if (tolower(stat_in) == "median") {
    # Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ bootstrap(., median, R = nboots))
    
    # Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp_pct <- effect_series_postonly %>%
      mutate(pct_eff = (!!as.name(outcome_var) / !!as.name(pred_var)) - 1) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["pct_eff"]]) %>%
      furrr::future_map(~ bootstrap(., median, R = nboots))
  }
  
  statname_abs <- paste("bootstrap_", stat_in, "_abs_tot", sep = "")
  statname_pct <- paste("bootstrap_", stat_in, "_pct_tot", sep = "")
  obs_statname_abs <- paste("observed_", stat_in, "_abs_tot", sep = "")
  obs_statname_pct <- paste("observed_", stat_in, "_pct_tot", sep = "")
  bootstrap_comp_ci_abs <- furrr::future_pmap(list(bootstrap_comp_abs, list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = effect_series_postonly %>%
             distinct(post_period_t) %>% 
             pull()) %>%
    mutate(
      !!as.name(statname_abs) := bootstrap_comp_abs %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist(),
      !!as.name(obs_statname_abs) := bootstrap_comp_abs %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Observed"]]) %>%
        unlist()
    )
  names(bootstrap_comp_ci_abs)[1:2] <- c("bootstrap_lb_abs_tot", "bootstrap_ub_abs_tot")
  
  bootstrap_comp_ci_pct <- furrr::future_pmap(list(bootstrap_comp_pct, list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = effect_series_postonly %>%
             distinct(post_period_t) %>%
             pull()) %>%
    mutate(
      !!as.name(statname_pct) := bootstrap_comp_pct %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist(),
      !!as.name(obs_statname_pct) := bootstrap_comp_pct %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Observed"]]) %>%
        unlist()
    )
  names(bootstrap_comp_ci_pct)[1:2] <- c("bootstrap_lb_pct_tot", "bootstrap_ub_pct_tot")
  
  bootstrap_comp_ci <- bootstrap_comp_ci_abs %>% inner_join(bootstrap_comp_ci_pct, by = "post_period_t")
  
  
  stat_tot_sample <- effect_series_postonly %>%
    group_by(post_period_t) %>%
    summarise(treated_n = n()) %>%
    ungroup()
  
  
  
  if (compute_cf_eff) {
    bootstrap_ci_by_postperiod <- stat_tot_sample %>%
      left_join(
        bootstrap_comp_ci,
        by = "post_period_t"
      ) %>%
      left_join(
        cf_tot_df,
        by = "post_period_t"
      )
  }
  if (!compute_cf_eff) {
    bootstrap_ci_by_postperiod <- stat_tot_sample %>% left_join(
      bootstrap_comp_ci,
      by = "post_period_t"
    )
  }
  
  
  
  return(bootstrap_ci_by_postperiod)
}


#TODO(alexdkellogg): Should there be weights based on num of treated entries?
#  in one sample, could be 750 treated at post=15, vs 600 in the other
#  Should we jackknife the bias directly?
compute_jackknife_bias<-function(tot_list, alpha_ci = 0.95){
  tot_tib=tot_list %>% dplyr::bind_rows()

  effect_var_abs<-names(tot_tib %>%
                          dplyr::select(tidyselect::matches("observed.*abs")))
  effect_var_pct<-names(tot_tib %>%
                          dplyr::select(tidyselect::matches("observed.*pct")))

  cfvar_pct<-names(tot_tib %>%
                     dplyr::select(tidyselect::contains("pct_cf")))
  
  cfvar_abs<-setdiff(names(tot_tib %>%
                             dplyr::select(tidyselect::contains("cf"))),
                     cfvar_pct)
    



  jackknife_comp_abs <- tot_tib %>%
    dplyr::mutate(bias=!!as.name(effect_var_abs)-!!as.name(cfvar_abs)) %>%
    split(.[["post_period_t"]]) %>%
    furrr::future_map(~ .[["bias"]]) %>%
    furrr::future_map(~ jackknife(., mean))
  
  #Weighted average, but no variation per group...
  # dplyr::group_by(post_period_t) %>%
  #   dplyr::mutate(weighted_avg = sum(treated_n*!!as.name(effect_var_abs))/sum(treated_n) )%>%

  jackknife_comp_pct <- tot_tib %>%
    dplyr::mutate(bias=!!as.name(effect_var_pct)-!!as.name(cfvar_pct)) %>%
    split(.[["post_period_t"]]) %>%
    furrr::future_map(~ .[["bias"]]) %>%
    furrr::future_map(~ jackknife(., mean))
  

  #get the 95% CI for the jackknife over post treat time and dataset
  statname_abs <- paste("jackknife_abs_bias", sep = "")
  statname_pct <- paste("jackknife_pct_bias", sep = "")
  obs_statname_abs <- paste("observed_abs_bias", sep = "")
  obs_statname_pct <- paste("observed_pct_bias", sep = "")
  jackknife_comp_ci_abs <- furrr::future_pmap(list(jackknife_comp_abs, list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = tot_tib %>% 
             distinct(post_period_t) %>%
             pull()) %>%
    mutate(
      !!as.name(statname_abs) := jackknife_comp_abs %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist(),
      !!as.name(obs_statname_abs) := jackknife_comp_abs %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Observed"]]) %>%
        unlist()
    )
  
  names(jackknife_comp_ci_abs)[1:2] <- c("jackknife_lb_abs_bias", "jackknife_ub_abs_bias")
  
  jackknife_comp_ci_pct <- furrr::future_pmap(list(jackknife_comp_pct,
                                                   list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = tot_tib %>%
             distinct(post_period_t) %>%
             pull()) %>%
    mutate(
      !!as.name(statname_pct) := jackknife_comp_pct %>%
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist(),
      !!as.name(obs_statname_pct) := jackknife_comp_pct %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Observed"]]) %>%
        unlist()
    )
  
  names(jackknife_comp_ci_pct)[1:2] <- c("jackknife_lb_pct_bias", "jackknife_ub_pct_bias")
  
  jackknife_comp_ci <- jackknife_comp_ci_abs %>% inner_join(jackknife_comp_ci_pct, by = "post_period_t") %>%
    dplyr::select(post_period_t, dplyr::everything()) %>%
    dplyr::mutate(cf_bias_pct=0, cf_bias_abs=0)
  
  

  
  return(jackknife_comp_ci)
  
}