library(CausalImpact)
library(Matrix)
library(tsfeatures)
library(tsibble)
library(ggfortify)
library(gsynth)
library(augsynth)
library(tidyverse)
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

# New fixes: recreating "point effect" and "pct.effect" variables within this code, wnat to allow that as an option
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
          mean_cf_tot = mean(!!as.name(outcome_var)- !!as.name(counterfac_var)),
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
  
  # store the bootstrap metric by post treatment period
  if (tolower(stat_in) == "aggregate") {
    # Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ jackknife(., sum))
    
    jackknife_comp_pct <- effect_series_postonly %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ jackknife(., function(x) {
        (sum(x[[outcome_var]]) / sum(x[[pred_var]])) - 1
      }))
  }
  
  if (tolower(stat_in) == "mean") {
    # Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ jackknife(., mean))
    
    # Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp_pct <- effect_series_postonly %>%
      mutate(pct_eff = (!!as.name(outcome_var) / !!as.name(pred_var)) - 1) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["pct_eff"]]) %>%
      furrr::future_map(~ jackknife(., mean))
  }
  
  if (tolower(stat_in) == "median") {
    # Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ jackknife(., median))
    
    # Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp_pct <- effect_series_postonly %>%
      mutate(pct_eff = (!!as.name(outcome_var) / !!as.name(pred_var)) - 1) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["pct_eff"]]) %>%
      furrr::future_map(~ jackknife(., median))
  }
  
  browser()
  statname_abs <- paste("jackknife_", stat_in, "_abs_tot", sep = "")
  statname_pct <- paste("jackknife_", stat_in, "_pct_tot", sep = "")
  obs_statname_abs <- paste("observed_", stat_in, "_abs_tot", sep = "")
  obs_statname_pct <- paste("observed_", stat_in, "_pct_tot", sep = "")
  jackknife_comp_ci_abs <- furrr::future_pmap(list(jackknife_comp_abs, list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = effect_series_postonly %>% 
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
  names(jackknife_comp_ci_abs)[1:2] <- c("jackknife_lb_abs_tot", "jackknife_ub_abs_tot")
  
  jackknife_comp_ci_pct <- furrr::future_pmap(list(jackknife_comp_pct,
                                                   list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = effect_series_postonly %>%
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




# Uncertainty Bootstrap
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
