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



#################################################################
# Metrics
#################################################################

compute_avg_metric_per_t <- function(estimated_series_df, time_var = "period", outcome_var = "response", prediction_var = "point.pred",
                                     counterfac_var = "counter_factual", treat_period_var = "Treatment_Period",
                                     metric_str = "mae", pct_eff_flag = F) {
  
  # Create the gap plot for each of the post treatment periods
  
  # Args
  # estimated_series_df: long-form dataframe with predicted counterfactual, actual counterfactual, period, and Treatment Period by ID.
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # outcome_var: the true y var for the time series
  # prediction_var: predicted counterfactual by time and id
  # counterfac_var: true counterfactual by time and id
  # treat_period_var: str name of the column which indicates the treatment period for each id_var (at each time_var)
  # metric_str: string indicating desired metric. can be mae, mse, or both
  # pct_eff_flag: binary flag for whether the percent error should be computed
  
  # Output
  # Dataframe containing number of rows equal to the longest post-treat period
  # and for each post treat period (from 0 -- time of treat, to max), an average of the metric for all observations
  # that experienced that particular post treat period (eg treated in t=5, total T=15 means 10 post treat periods)
  
  # absolute error computation
  if (pct_eff_flag == FALSE) {
    if (tolower(metric_str) == "mse") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_error_sq = (!!as.name(prediction_var) - !!as.name(counterfac_var))^2) %>%
        group_by(post_period_t) %>%
        summarise(MSE_TE_byT = mean(indiv_error_sq))
    }
    
    if (tolower(metric_str) == "mae") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error = abs(!!as.name(prediction_var) - !!as.name(counterfac_var))) %>%
        group_by(post_period_t) %>%
        summarise(MAE_TE_byT = mean(indiv_t_error))
    }
    
    
    if (tolower(metric_str) == "both") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(
          indiv_t_error = abs(!!as.name(prediction_var) - !!as.name(counterfac_var)),
          indiv_error_sq = (!!as.name(prediction_var) - !!as.name(counterfac_var))^2
        ) %>%
        group_by(post_period_t) %>%
        summarise(
          MAE_TE_byT = mean(indiv_t_error),
          MSE_TE_byT = mean(indiv_error_sq)
        )
    }
  }
  
  # percent error computation
  if (pct_eff_flag) {
    if (tolower(metric_str) == "mse") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueT- EstTE=(Target/counterFac)-1-(Target/Pred-1)
        mutate(indiv_error_sq = (!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))^2) %>%
        group_by(post_period_t) %>%
        summarise(MSE_TE_byT = mean(indiv_error_sq))
    }
    
    if (tolower(metric_str) == "mae") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error = abs(!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))) %>%
        group_by(post_period_t) %>%
        summarise(MAE_TE_byT = mean(indiv_t_error))
    }
    
    
    if (tolower(metric_str) == "both") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(
          indiv_t_error = abs(!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var))),
          indiv_error_sq = (!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))^2
        ) %>%
        group_by(post_period_t) %>%
        summarise(
          MAE_TE_byT = mean(indiv_t_error),
          MSE_TE_byT = mean(indiv_error_sq)
        )
    }
  }
  
  
  
  
  return(avg_metric_series)
}



compute_avg_metric <- function(estimated_series_df, time_var = "period", outcome_var = "response", prediction_var = "point.pred",
                               counterfac_var = "counter_factual", treat_period_var = "Treatment_Period",
                               metric_str = "mae", pct_eff_flag = F) {
  # Estimates the relevant metric for the predicted Treatment Effect by comparing the predicted and true counterfactuals
  
  # Args
  # estimated_series_df: long-form dataframe with predicted counterfactual, actual counterfactual, period, and Treatment Period by ID.
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # outcome_var: the true y var for the time series
  # prediction_var: predicted counterfactual by time and id
  # counterfac_var: true counterfactual by time and id
  # treat_period_var: str name of the column which indicates the treatment period for each id_var (at each time_var)
  # metric_str: string indicating desired metric. can be mae, mse, or both
  # pct_eff_flag: binary flag for whether the percent error should be computed
  
  # Output
  # tibble containing an average of the metric for all observations, over all T
  
  # absolute error computation
  if (pct_eff_flag == FALSE) {
    if (tolower(metric_str) == "mse") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_error_sq = (!!as.name(prediction_var) - !!as.name(counterfac_var))^2) %>%
        summarise(MSE_TE = mean(indiv_error_sq))
    }
    
    if (tolower(metric_str) == "mae") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error = abs(!!as.name(prediction_var) - !!as.name(counterfac_var))) %>%
        summarise(MAE_TE = mean(indiv_t_error))
    }
    
    
    if (tolower(metric_str) == "both") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(
          indiv_t_error = abs(!!as.name(prediction_var) - !!as.name(counterfac_var)),
          indiv_error_sq = (!!as.name(prediction_var) - !!as.name(counterfac_var))^2
        ) %>%
        summarise(
          MAE_TE = mean(indiv_t_error),
          MSE_TE = mean(indiv_error_sq)
        )
    }
  }
  
  # percent error computation
  if (pct_eff_flag) {
    if (tolower(metric_str) == "mse") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueT- EstTE=(Target/counterFac)-1-(Target/Pred-1)
        mutate(indiv_error_sq = (!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))^2) %>%
        summarise(MSE_TE = mean(indiv_error_sq))
    }
    
    if (tolower(metric_str) == "mae") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error = abs(!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))) %>%
        summarise(MAE_TE = mean(indiv_t_error))
    }
    
    
    if (tolower(metric_str) == "both") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(
          indiv_t_error = abs(!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var))),
          indiv_error_sq = (!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))^2
        ) %>%
        summarise(
          MAE_TE = mean(indiv_t_error),
          MSE_TE = mean(indiv_error_sq)
        )
    }
  }
  
  
  
  
  return(avg_metric)
}


compute_tot_coverage <- function(bootstrapped_series_df, time_var = "post_period_t", ub_var = "jackknife_ub_pct_tot",
                                 lb_var = "jackknife_lb_pct_tot", counterfac_var = "mean_pct_cf_tot") {
  coverage_out <- bootstrapped_series_df %>%
    mutate(
      contains_truth = as.numeric(!!as.name(counterfac_var) <= !!as.name(ub_var) &
                                    !!as.name(counterfac_var) >= !!as.name(lb_var))
    ) %>%
    pull(contains_truth) %>%
    mean()
  
  return(coverage_out)
}