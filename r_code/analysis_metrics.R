pacman::p_load(dplyr)
#################################################################
# Metrics
#################################################################


metric_gen_helper<-function(input_tib, t_var_inp, treat_t_inp,
                            pred_inp, cf_inp, outcome_inp,
                            pct_flag_inp, data_type=c("raw", "tot")){
  #Helper function to our metric computation procedures, handles individual
  #level computation
  
  #Args
  #input_tib: original data
  #t_var_inp: string name of the time column.
  #treat_t_inp: string name of the column indicating treatment start period.
  #outcome_inp: string name of the outcome variable of interest.
  #cf_inp: string name of the counterfactual outcome
  #pred_inp: string name for the predicted outcome column
  
  #output: tibble with the metrics by entry and time since treatment
  
  #Create tibble with metrics (MAE and MSE) for each entry-time combination
  if(data_type=="raw"){
    input_tib=input_tib %>%
      dplyr::mutate(post_period_t = !!as.name(t_var_inp) - 
                      !!as.name(treat_t_inp))
  }
  avg_metric <- input_tib %>%
    dplyr::filter(post_period_t >= 0) %>% 
    dplyr::mutate(indiv_t_error = dplyr::case_when(
      !pct_flag_inp~abs(!!as.name(pred_inp) - !!as.name(cf_inp)),
      pct_flag_inp~abs(!!as.name(outcome_inp) / !!as.name(cf_inp) -
                         (!!as.name(outcome_inp) / !!as.name(pred_inp)))
    ),
    indiv_error_sq=dplyr::case_when(
      !pct_flag_inp~(!!as.name(pred_inp) - !!as.name(cf_inp))^2,
      pct_flag_inp~(!!as.name(outcome_inp) / !!as.name(cf_inp) - 
                      (!!as.name(outcome_inp) / !!as.name(pred_inp)))^2
    ))
  
  return(avg_metric)
}



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
  # Tibble containing number of rows equal to the longest post-treat period
  # and for each post treat period (from 0 -- time of treat, to max), an average of the metric for all observations
  # that experienced that particular post treat period (eg treated in t=5, total T=15 means 10 post treat periods)
  
  #Create tibble with metrics (MAE and MSE) for each entry-time combination
  avg_metric_series <- 
    metric_gen_helper(input_tib=estimated_series_df, t_var_inp=time_var, 
                      treat_t_inp=treat_period_var,pred_inp=prediction_var,
                      cf_inp=counterfac_var, outcome_inp=outcome_var,
                      pct_flag_inp=pct_eff_flag) %>%
    dplyr::group_by(post_period_t)
  
  #Compute the average metrics over te entries for each time period
  metric_series_output= switch(tolower(metric_str),
                        "mse" = avg_metric_series %>%
                          dplyr::summarise(MSE_TE_byT = mean(indiv_error_sq)),
                        "mae" = avg_metric_series %>% 
                          dplyr::summarise(MAE_TE_byT = mean(indiv_t_error)),
                        "both" = avg_metric_series %>% 
                          dplyr::summarise(MAE_TE_byT = mean(indiv_t_error),
                                           MSE_TE_byT = mean(indiv_error_sq))
  )
  
  
  return(metric_series_output)
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
  

  #Create tibble with metrics (MAE and MSE) for each entry-time combination
  avg_metric <-
    metric_gen_helper(input_tib=estimated_series_df, t_var_inp=time_var, 
                      treat_t_inp=treat_period_var,pred_inp=prediction_var,
                      cf_inp=counterfac_var, outcome_inp=outcome_var,
                      pct_flag_inp=pct_eff_flag)
  
  metric_output= switch(tolower(metric_str),
    "mse" = avg_metric %>%
      dplyr::summarise(MSE_TE = mean(indiv_error_sq)),
    "mae" = avg_metric %>% 
      dplyr::summarise(MAE_TE = mean(indiv_t_error)),
    "both" = avg_metric %>% 
      dplyr::summarise(MAE_TE = mean(indiv_t_error),
                       MSE_TE = mean(indiv_error_sq))
  )
  
  
  
  return(metric_output)

  
}
