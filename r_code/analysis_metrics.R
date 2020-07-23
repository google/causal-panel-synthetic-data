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



#################################################################
# Metrics
#################################################################

metric_gen_helper<-function(input_tib, t_var_inp, treat_t_inp,
                            pred_inp, cf_inp, outcome_inp,
                            pct_flag_inp){
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
  avg_metric <- input_tib %>%
    dplyr::mutate(post_period_t = !!as.name(t_var_inp) - 
                    !!as.name(treat_t_inp)) %>%
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

#TODO(alexdkellogg): update for over seed use.
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




compute_jackknife_metrics<-function(estimated_tib_list, time_var = "period", 
                            outcome_var = "response",
                            prediction_var = "point.pred",
                            counterfac_var = "counter_factual",
                            treat_period_var = "Treatment_Period",
                            pct_eff_flag = T,
                            alpha_ci=0.95,
                            horizon=NULL){
  
  #horizon: if null, average over all t. Otherwise, compute grouped avg up to horizon
  #bind our data over each replicate sample
  estimated_tib=estimated_tib_list %>% dplyr::bind_rows()
  #store the metrics for each unit time combination
  avg_metric_series <- 
    metric_gen_helper(input_tib=estimated_tib, t_var_inp=time_var, 
                      treat_t_inp=treat_period_var,pred_inp=prediction_var,
                      cf_inp=counterfac_var, outcome_inp=outcome_var,
                      pct_flag_inp=pct_eff_flag)
  #Split by time if a horizon is specified
  if(is.null(horizon)){
    effects_to_jk<-avg_metric_series %>%
      dplyr::select(post_period_t, indiv_t_error, indiv_error_sq)
    
    jk_comp_mae=effects_to_jk %>% 
      dplyr::pull(indiv_t_error) %>%
      resample::jackknife(., mean, args.stat = list(na.rm=T))
    
    jk_comp_mse=effects_to_jk %>% 
      dplyr::pull(indiv_error_sq) %>%
      resample::jackknife(mean, args.stat = list(na.rm=T))

    jackknife_comp_ci=tibble::tibble(
      jackknife_lb_mae=jk_comp_mae$stats$Mean-1.96*jk_comp_mae$stats$SE,
      jackknife_ub_mae=jk_comp_mae$stats$Mean+1.96*jk_comp_mae$stats$SE,
      jackknife_mae=jk_comp_mae$stats$Mean,
      observed_mae=jk_comp_mae$stats$Observed,
      jackknife_lb_mse=jk_comp_mse$stats$Mean-1.96*jk_comp_mse$stats$SE,
      jackknife_ub_mse=jk_comp_mse$stats$Mean+1.96*jk_comp_mse$stats$SE,
      jackknife_mse=jk_comp_mse$stats$Mean,
      observed_mse=jk_comp_mse$stats$Observed)
  }
  else{
    #Take advantage of parallel processing if possible
    effects_to_jk=avg_metric_series %>%
      dplyr::select(post_period_t, indiv_t_error, indiv_error_sq) %>%
      dplyr::filter(post_period_t<horizon) %>%
      split(.[["post_period_t"]])
  
  jk_comp_mae= effects_to_jk %>%
    furrr::future_map(~ .[["indiv_t_error"]]) %>% 
    furrr::future_map(~ resample::jackknife(
      ., mean, args.stat = list(na.rm=T)) )
  
  jk_comp_mse=effects_to_jk %>%
    furrr::future_map(~ .[["indiv_error_sq"]]) %>% 
    furrr::future_map(~ resample::jackknife(
      ., mean,args.stat = list(na.rm=T)))
  
  jk_comp_ci_mae= format_jackknife(jk_est=jk_comp_mae,ci_in=alpha_ci, 
                                   col_name= "mae")
  
  jk_comp_ci_mse= format_jackknife(jk_est=jk_comp_mse,ci_in=alpha_ci, 
                                   col_name= "mse")
  jackknife_comp_ci <- jk_comp_ci_mae %>% 
    dplyr::inner_join(jk_comp_ci_mse, by = "post_period_t")
  }
  return(jackknife_comp_ci)
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

