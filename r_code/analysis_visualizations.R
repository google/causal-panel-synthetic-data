pacman::p_load(dplyr, ggplot2)


#################################################################
# Plots
#################################################################



# print PDF Likely does not work with future_map???
create_gap_ci_plot <- function(bootstrapped_effects_df, time_var = "post_period_t",
                               cf_plot = T, pct_flag=T,
                               print_to_pdf = NULL, plot_title = NULL, plot_x_lab = NULL, plot_y_lab = NULL) {
  
  # Computes bootstrap estimates (from resample package) of the treatment effect on the treated, by period
  
  # Args
  # bootstrapped_effects_df: data on the bootstrapped (and potentially counterfactual) effects for each post treat period,
  # from compute_ci_bounds_bootstrap or jackknife version
  # time_var:column name of numeric period number indicating the post treatment time period, in increasing order (eg 0 is the first time)
  # effect_var: string indicating the effect we wish to plot, either "jackknife_mean_tot" or "jackknife_median_tot" works, depending on input
  # upper_ci: string indicating thename of the upper CI variable in bootstrapped_effects_df
  # lower_ci: string indicating thename of the lower CI variable in bootstrapped_effects_df
  # cf_plot: boolean flag for whether a counterfactual effect should be plotted
  # cf_var: string variable name of the counterfactual effect to plot, typically "mean_cf_tot or "median_cf_tot" from bootstrapped_effects_df
  # print_to_pdf: string, file path where the pdfs should be printed.
  # plot_title: title of the plot to be printed
  # plot_x_lab: title of the x axis label
  # plot_y_lab: title of the y axis label
  
  # Output
  # gap plot, potentially with counterfactual predictions
  
  #Automatically get the names of the variables to plot
  effect_type=ifelse(pct_flag,"pct", "abs")
  bootstrapped_effects_df=bootstrapped_effects_df %>%
    dplyr::select(tidyselect::all_of(time_var)) %>%
    dplyr::bind_cols(
      bootstrapped_effects_df %>%
        dplyr::select(tidyselect::matches("tot|bias")) %>%
        dplyr::select(tidyselect::contains(effect_type))
    )
    
  
  
  lower_ci=names(bootstrapped_effects_df %>%
                     dplyr::select(tidyselect::contains("lb")))
  upper_ci=names(bootstrapped_effects_df %>%
                   dplyr::select(tidyselect::contains("ub")))
  effect_var=setdiff(names(bootstrapped_effects_df %>%
                             dplyr::select(
                               tidyselect::matches("jackknife|bootstrap"))), 
                     c(lower_ci, upper_ci))
  
  if (!is.null(print_to_pdf)) {
    pdf(print_to_pdf)
  }

  if (cf_plot) {
    cf_var=names(bootstrapped_effects_df %>%
                   dplyr::select(tidyselect::contains("cf")))
    plot_out <- bootstrapped_effects_df %>%
      ggplot2::ggplot(aes(x = !!as.name(time_var), y = !!as.name(effect_var), color = "Estimate")) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(aes(ymin = !!as.name(lower_ci), ymax = !!as.name(upper_ci)), alpha = 0.3, color = NA) +
      ggplot2::geom_line(aes(y = !!as.name(cf_var), color = "True")) +
      ggplot2::ggtitle(ifelse(is.null(plot_title), paste("Gap Plot of", effect_var), plot_title)) +
      ggplot2::labs(
        x = ifelse(is.null(plot_x_lab), time_var, plot_x_lab),
        y = ifelse(is.null(plot_y_lab), effect_var, plot_y_lab)
      ) +
      ggplot2::scale_colour_manual(
        name = "",
        values = c("Estimate" = "black", "True" = "red"), labels = c("Estimate", "True")
      )
  }
  
  if (!cf_plot) {
    plot_out <- bootstrapped_effects_df %>%
      ggplot2::ggplot(aes(x = !!as.name(time_var), y = !!as.name(effect_var))) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(aes(ymin = !!as.name(lower_ci), ymax = !!as.name(upper_ci)), alpha = 0.3) +
      ggplot2::ggtitle(ifelse(is.null(plot_title), paste("Gap Plot of", effect_var), plot_title)) +
      ggplot2::labs(
        x = ifelse(is.null(plot_x_lab), time_var, plot_x_lab),
        y = ifelse(is.null(plot_y_lab), effect_var, plot_y_lab)
      )
  }
  
  
  if (!is.null(print_to_pdf)) {
    print(plot_out)
    dev.off()
  }
  return(plot_out)
}


plot_individual_cf<-function(est_data_inp, time_var="period", id_var="entry",
                             outcome_var="response",
                             prediction_var="point.pred",
                             treat_period_var="Treatment_Period"){
  
  big_ids=est_data_inp %>%
    dplyr::filter(!!as.name(time_var)==1) %>%
    dplyr::filter(!!as.name(outcome_var)>
                    quantile(!!as.name(outcome_var),0.95)) %>%
    dplyr::select(tidyselect::all_of(c(id_var, treat_period_var)))
  
  random_ids=est_data_inp %>%
    dplyr::filter(!!as.name(time_var)==1) %>%
    dplyr::filter(!!as.name(id_var) %in% 
                    setdiff(!!as.name(id_var), big_ids %>% 
                              dplyr::pull(!!as.name(id_var)))) %>%
    dplyr::sample_n(3) %>%
    dplyr::select(tidyselect::all_of(c(id_var, treat_period_var)))
  
  indiv_ids=big_ids %>% 
    dplyr::bind_rows(random_ids) %>% 
    dplyr::pull(!!as.name(id_var))
  
  indiv_treat=big_ids %>% 
    dplyr::bind_rows(random_ids) %>% 
    dplyr::pull(!!as.name(treat_period_var))
  
  plot_out <- furrr::future_map2(.x=indiv_ids, .y=indiv_treat,
                                 .f=individual_cf_helper, data_inp=est_data_inp)
  
  return(plot_out)
  
}


individual_cf_helper<-function(data_inp, id_inp, treat_inp,
                             id_var="entry",time_var="period" ,
                             outcome_var="response",
                             prediction_var="point.pred"){
  data_inp %>%
    dplyr::filter(!!as.name(id_var)==id_inp) %>%
    ggplot2::ggplot(aes(x = !!as.name(time_var), y = !!as.name(outcome_var))) +
    ggplot2::geom_line(aes(color="Outcome")) +
    ggplot2::geom_line(aes(y=!!as.name(prediction_var), color="Predicted" ))+
    ggplot2::geom_vline(xintercept = treat_inp, color="red") +
    ggplot2::scale_color_manual(name="Type", values = c("black", "darkgray"))+
    ggplot2::theme_bw()+
    ggplot2::labs(x="Time", y="Outcome", 
                  title = "Counterfactual vs Outcome Series",
                  subtitle = paste("ID=",id_inp))
}





plot_full_gap<-function(est_data_inp,  time_var="period", id_var="entry",
                        outcome_var="response", prediction_var="point.pred",
                        treat_period_var="Treatment_Period", start_plot=-20,
                        pct_flag=F){
  single_dataset=F
  #Turn a single dataframe into a list, but note the type of data
  if(inherits(est_data_inp, "data.frame")){
    single_dataset=T
    est_data_inp=list(est_data_inp)
  }
  #Step 1: make a call to the jackknife function 
  gap_full=furrr::future_map(est_data_inp, compute_tot_se_jackknife, 
                             stat_in="mean", post_treat_only=F)
  gap_full_forplot=furrr::future_map(gap_full, 
                                     function(x){
                                       return(
                                         x %>% 
                                           dplyr::filter(post_period_t>=start_plot)
                                         )})
  #Step 2: make a call to the plot with updated title
  #If we are simulating (or have multiple data with draws from eps), 
  #compute the bias first. If we are only dealing with one dataset, print as is.
  if(single_dataset==T){
    plot_out<- create_gap_ci_plot(gap_full_forplot[[1]], 
                       plot_title="Gap Plot", 
                       plot_x_lab="Time Period",
                       plot_y_lab="ToT", pct_flag = pct_flag)+
      geom_vline(xintercept = 0, color="blue")
  }
  else{
    bias_forplot=compute_jackknife_bias(gap_full_forplot,horizon = max(gap_full_forplot[[1]]$post_period_t) )
    plot_out <- create_gap_ci_plot(bias_forplot, 
                                   plot_title="Gap Plot", 
                                   plot_x_lab="Time Period",
                                   plot_y_lab="ToT", pct_flag = pct_flag)+
      geom_vline(xintercept = 0, color="blue")
  }
  
  return(plot_out)
  
}




#TODO(alexdkellogg): not sure if we need this anymore
plot_tot_bias_per_t <- function(estimated_series_df, time_var = "period", id_var = "entry", outcome_var = "response", prediction_var = "point.pred",
                                counterfac_var = "counter_factual", treat_period_var = "Treatment_Period",
                                max_post_t = 12, pct_eff_flag = F, plot_title = NULL) {
  
  # Create the gap plot for each of the post treatment periods
  
  # Args
  # estimated_series_df: long-form dataframe with predicted counterfactual, actual counterfactual, period, and Treatment Period by ID.
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # id_var: column name of the unit id
  # outcome_var: the true y var for the time series
  # prediction_var: predicted counterfactual by time and id
  # counterfac_var: true counterfactual by time and id
  # treat_period_var: str name of the column which indicates the treatment period for each id_var (at each time_var)
  # max_post_t: numeric variable, the number of time periods over which we are averaging per unit
  # pct_eff_flag: binary flag for whether the percent error should be computed
  # plot_title: string title of the plot
  
  # Output
  # Dataframe containing number of rows equal to the longest post-treat period
  # and for each post treat period (from 0 -- time of treat, to max), an average of the metric for all observations
  # that experienced that particular post treat period (eg treated in t=5, total T=15 means 10 post treat periods)
  
  # absolute error computation
  if (pct_eff_flag == FALSE) {
    avg_bias_series <- estimated_series_df %>%
      mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
      filter(post_period_t >= 0, post_period_t <= max_post_t) %>%
      group_by(!!as.name(id_var)) %>%
      summarise(bias_per_id = mean(!!as.name(prediction_var) - !!as.name(counterfac_var))) %>%
      ggplot2::ggplot(aes(x = bias_per_id)) +
      ggplot2::geom_density(fill = "blue", alpha = 0.4) +
      theme_bw() +
      theme(
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
      ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      scale_x_continuous(name = "ATT by ID") +
      ggplot2::ggtitle(ifelse(is.null(plot_title), paste("ToT Bias by Unit, averaged over", max_post_t, "post-treat Periods"), plot_title))
  }
  
  # percent error computation
  if (pct_eff_flag) {
    avg_bias_series <- estimated_series_df %>%
      mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
      filter(post_period_t >= 0) %>%
      group_by(!!as.name(id_var)) %>%
      summarise(bias_per_id = mean(!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))) %>%
      ggplot2::ggplot(aes(x = bias_per_id)) +
      ggplot2::geom_density(fill = "blue", alpha = 0.4) +
      theme_bw() +
      theme(
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
      ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      scale_x_continuous(name = "ATT by ID") +
      ggplot2::ggtitle(ifelse(is.null(plot_title), paste("ToT Pct Bias by Unit, averaged over", max_post_t, "post-treat Periods"), plot_title))
  }
  
  
  
  
  return(avg_bias_series)
}
