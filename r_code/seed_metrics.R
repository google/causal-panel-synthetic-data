source(here("r_code/analysis_metrics.R"))
pacman::p_load(dplyr, furrr, tidyr, resample, tibble)

compute_jackknife_bias<-function(tot_list, alpha_ci = 0.95, horizon=5){
  tot_tib=tot_list %>% 
    dplyr::bind_rows() %>%
    dplyr::filter(post_period_t<horizon)
  
  effect_var_abs<-names(tot_tib %>%
                          dplyr::select(tidyselect::matches("observed.*abs")))
  effect_var_pct<-names(tot_tib %>%
                          dplyr::select(tidyselect::matches("observed.*pct")))
  
  cfvar_pct<-names(tot_tib %>% dplyr::select(tidyselect::contains("pct_cf")))
  
  cfvar_abs<-setdiff(names(tot_tib %>%
                             dplyr::select(tidyselect::contains("cf"))),
                     cfvar_pct)
  
  jk_comp_abs <- tot_tib %>%
    dplyr::mutate(bias=!!as.name(effect_var_abs)-!!as.name(cfvar_abs)) %>%
    split(.[["post_period_t"]]) %>%
    furrr::future_map(~ .[["bias"]]) %>%
    furrr::future_map(~ jackknife(., mean))
  
  #Weighted average, but no variation per group...
  # dplyr::group_by(post_period_t) %>%
  #   dplyr::mutate(weighted_avg = sum(treated_n*!!as.name(effect_var_abs))/sum(treated_n) )%>%
  
  jk_comp_pct <- tot_tib %>%
    dplyr::mutate(bias=!!as.name(effect_var_pct)-!!as.name(cfvar_pct)) %>%
    split(.[["post_period_t"]]) %>%
    furrr::future_map(~ .[["bias"]]) %>%
    furrr::future_map(~ jackknife(., mean))
  
  post_range=c(min(tot_tib$post_period_t),max(tot_tib$post_period_t) )
  #Gather the bounds on the confidence intervals and format as tibble
  jk_comp_ci <-format_jackknife(
    jk_est=jk_comp_abs,ci_in=alpha_ci,col_name="abs_bias",range=post_range) %>%
    dplyr::inner_join(
      format_jackknife( jk_est=jk_comp_pct,ci_in=alpha_ci, 
                        col_name="pct_bias", range=post_range), 
      by="post_period_t") %>%
    dplyr::select(post_period_t, dplyr::everything()) %>%
    dplyr::mutate(cf_bias_pct=0, cf_bias_abs=0)
  
  return(jk_comp_ci)
}


compute_tot_variance<-function(tot_list, horizon=5){
  tot_tib=tot_list %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(post_period_t<horizon)
  
  effect_var_abs<-names(tot_tib %>%
                          dplyr::select(tidyselect::matches("observed.*abs")))
  effect_var_pct<-names(tot_tib %>%
                          dplyr::select(tidyselect::matches("observed.*pct")))
  
  cfvar_pct<-names(tot_tib %>% dplyr::select(tidyselect::contains("pct_cf")))
  
  cfvar_abs<-setdiff(names(tot_tib %>%
                             dplyr::select(tidyselect::contains("cf"))),
                     cfvar_pct)
  
  jk_comp_abs <- tot_tib %>%
    split(.[["post_period_t"]]) %>%
    furrr::future_map(~ .[[effect_var_abs]]) %>%
    furrr::future_map(~ jackknife(., var))
  
  #Weighted average, but no variation per group...
  # dplyr::group_by(post_period_t) %>%
  #   dplyr::mutate(weighted_avg = sum(treated_n*!!as.name(effect_var_abs))/sum(treated_n) )%>%
  
  jk_comp_pct <- tot_tib %>%
    split(.[["post_period_t"]]) %>%
    furrr::future_map(~ .[[effect_var_pct]]) %>%
    furrr::future_map(~ jackknife(., var))
  
  #Gather the bounds on the confidence intervals and format as tibble
  post_range=c(min(tot_tib$post_period_t),max(tot_tib$post_period_t) )
  jk_comp_ci <-format_jackknife(
    jk_est=jk_comp_abs,ci_in=0.95,col_name="abs_tot_var", range=post_range) %>%
    dplyr::inner_join(
      format_jackknife(jk_est=jk_comp_pct,ci_in=0.95, 
                       col_name="pct_tot_var", range=post_range), 
      by="post_period_t") %>%
    dplyr::select(post_period_t, jackknife_abs_tot_var, observed_abs_tot_var,
                  jackknife_pct_tot_var,observed_pct_tot_var)
  
  return(jk_comp_ci)
}


#TODO(alexdkellogg): update for over seed use.
compute_tot_coverage <- function(tot_list, horizon=5) {
  
  tot_tib=tot_list %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(post_period_t<horizon)
  
  ub_var_abs<-names(tot_tib %>%
                          dplyr::select(tidyselect::matches("ub.*abs")))
  lb_var_abs<-names(tot_tib %>%
                      dplyr::select(tidyselect::matches("lb.*abs")))
  lb_var_pct<-names(tot_tib %>%
                          dplyr::select(tidyselect::matches("lb.*pct")))
  ub_var_pct<-names(tot_tib %>%
                          dplyr::select(tidyselect::matches("ub.*pct")))
  
  cfvar_pct<-names(tot_tib %>% dplyr::select(tidyselect::contains("pct_cf")))
  
  cfvar_abs<-setdiff(names(tot_tib %>%
                             dplyr::select(tidyselect::contains("cf"))),
                     cfvar_pct)
  
  coverage_out <- tot_tib %>%
    dplyr::group_by(post_period_t) %>%
    dplyr::mutate(
      contains_truth_abs = as.numeric(!!as.name(cfvar_abs) <= 
                                    !!as.name(ub_var_abs) &
                                    !!as.name(cfvar_abs) >= 
                                    !!as.name(lb_var_abs)),
      contains_truth_pct = as.numeric(!!as.name(cfvar_pct) <= 
                                        !!as.name(ub_var_pct) &
                                        !!as.name(cfvar_pct) >= 
                                        !!as.name(lb_var_pct)),
    ) %>%
    dplyr::summarise(coverage_abs=mean(contains_truth_abs),
                     coverage_pct=mean(contains_truth_pct))
  
  return(coverage_out)
}



#TODO(alexdkellogg): clean this up -- was written for c(outcome_var = "response",
#    prediction_var = "point.pred",counterfac_var = "counter_factual"); hacked 
#    to convert it to RMSE of ATT rather than Y_hat
compute_jackknife_metrics<-function(estimated_tib_list, 
                                    data_type=c("tot", "raw"),
                                    time_var = "period", 
                                    prediction_var = "observed_mean_abs_tot", 
                                    counterfac_var = "mean_abs_cf_tot", 
                                    outcome_var = "mean_abs_cf_tot",
                                    treat_period_var = "Treatment_Period",
                                    pct_eff_flag = F,
                                    alpha_ci=0.95,
                                    horizon=5){
  # if (missing(data_type)) {
  #   stop("Please input argument data_type as `raw` or `tot`.")
  # } else{
  #   data_type <- match.arg(data_type)
  # }
  
  #horizon: if null, average over all t. Otherwise, compute grouped avg up to horizon
  #bind our data over each replicate sample
  estimated_tib=estimated_tib_list %>% dplyr::bind_rows()
  #store the metrics for each unit time combination
  avg_metric_series <- 
    metric_gen_helper(input_tib=estimated_tib, t_var_inp=time_var, 
                      treat_t_inp=treat_period_var,pred_inp=prediction_var,
                      cf_inp=counterfac_var, outcome_inp=outcome_var,
                      pct_flag_inp=pct_eff_flag,data_type=data_type)
  #Split by time if a horizon is specified
  if(is.null(horizon)){
    effects_to_jk<-avg_metric_series %>%
      dplyr::select(post_period_t, indiv_t_error, indiv_error_sq)
    
    # jk_comp_mae=effects_to_jk %>% 
    #   dplyr::pull(indiv_t_error) %>%
    #   resample::jackknife(., mean, args.stat = list(na.rm=T))
    
    jk_comp_mse=effects_to_jk %>%
      dplyr::pull(indiv_error_sq) %>%
      resample::jackknife(mean, args.stat = list(na.rm=T))
    
    jackknife_comp_ci=tibble::tibble(
      # jackknife_mae=jk_comp_mae$stats$Mean,
      # observed_mae=jk_comp_mae$stats$Observed,
      jackknife_rmse=sqrt(jk_comp_mse$stats$Mean),
      observed_rmse=sqrt(jk_comp_mse$stats$Observed))
  }
  else{
    #Take advantage of parallel processing if possible
    effects_to_jk=avg_metric_series %>%
      dplyr::select(post_period_t, indiv_t_error, indiv_error_sq) %>%
      dplyr::filter(post_period_t<horizon) %>%
      split(.[["post_period_t"]])
    
    # jk_comp_mae= effects_to_jk %>%
    #   furrr::future_map(~ .[["indiv_t_error"]]) %>% 
    #   furrr::future_map(~ resample::jackknife(
    #     ., mean, args.stat = list(na.rm=T)) )
    
    jk_comp_mse=effects_to_jk %>%
      furrr::future_map(~ .[["indiv_error_sq"]]) %>% 
      furrr::future_map(~ resample::jackknife(
        ., mean,args.stat = list(na.rm=T)))
    
    # jk_comp_ci_mae= format_jackknife(jk_est=jk_comp_mae,ci_in=alpha_ci, 
    #                                  col_name= "mae")
    # jackknife_comp_ci <- jk_comp_ci_mae %>% 
    #   dplyr::inner_join(jk_comp_ci_mse, by = "post_period_t")
    post_range=c(min(avg_metric_series$post_period_t),horizon-1)
    jk_comp_ci_mse= format_jackknife(jk_est=jk_comp_mse,ci_in=alpha_ci, 
                                     col_name= "mse", range=post_range) %>%
      dplyr::mutate(jackknife_rmse=sqrt(jackknife_mse),
                    observed_rmse=sqrt(observed_mse)) %>% 
      dplyr::select(post_period_t, jackknife_rmse, observed_rmse) 
    
    
  }
  return(jk_comp_ci_mse)
}




