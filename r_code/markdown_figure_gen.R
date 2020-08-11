source(here("r_code/analysis_visualizations.R"))
source(here("r_code/factor_dgp.R"))
source(here("r_code/seed_metrics.R"))

gt_set_up <- function(dataset_name, method_vec){
  load(here::here("Data", "Variations",dataset_name ))
  data_name=stringr::str_remove(tolower(dataset_name), ".rdata")
  
  metric_tibble_helper<-function(metric_name, method_name){
    tib_out=eval(as.name(paste(method_name,metric_name, sep="_"))) %>% 
      dplyr::mutate( method=method_name)
  }
  
  #set up the coverage tibble
  coverage_tib=
    lapply(method_vec, metric_tibble_helper, metric_name="coverage") %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(metric="coverage", metric_value=coverage_abs)%>%
    dplyr::select(post_period_t,metric, method, metric_value)
  
  #set up the rmse tibble
  rmse_tib=
    lapply(method_vec, metric_tibble_helper, metric_name="overall_metrics") %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(metric="rmse", metric_value=jackknife_rmse)%>%
    dplyr::select(post_period_t,metric, method, metric_value)
  
  #set up the bias tibble
  bias_tib=
    lapply(method_vec, metric_tibble_helper, metric_name="bias") %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(metric="bias", metric_value=jackknife_abs_bias) %>%
    dplyr::select(post_period_t, metric,metric_value, method)
  
  metrics_tib=coverage_tib %>%
    dplyr::bind_rows(rmse_tib) %>%
    dplyr::bind_rows(bias_tib) %>%
    tidyr::pivot_wider(names_from = c( method ), values_from=metric_value)
  
  return(metrics_tib)
}



markdown_estimation_output<-function(tib_to_gt_data, datalist, plot_indiv=T,
                                     method_vec){
  load(here::here("Data", "Variations",datalist ))
  
  data_name=stringr::str_remove(tolower(datalist), ".rdata")
  print(data_name)
  #print(all_dgp_params[[match(data_name,names(all_dgp_params))]])

  args.list <- c(lapply(method_vec, function(x){
    eval(as.name(paste(x,"bias", "plot", sep="_")))
  }) ,list(nrow=4,ncol=2))
  do.call(gridExtra::grid.arrange,args.list )

  
  indiv_cf_plotter<-function(method_name, rand_int, num_ex){
    lapply(
      plot_individual_cf(
        eval(as.name(paste(method_name,"est", sep="_")))[[rand_int]])[seq_len(num_ex)],
      function(x){
        x+ggplot2::labs(caption=method_name)
      })
  }
  
  
  random_data=sample(1:n_seeds,1)
  if(plot_indiv){
    all_cf_plots= lapply(method_vec, indiv_cf_plotter,
                         rand_int=random_data, num_ex=2)
    
    indiv_cf_combined=c(all_cf_plots)
    
    lapply(indiv_cf_combined,print)
  }
  
  gridExtra::grid.arrange(tsfeature_pc_by_treatment_plot(formatted_data[[random_data]])+
                            ggplot2::labs(caption = data_name), nrow=1)
  
  print(tsfeature_by_treat_df(formatted_data[[random_data]]) %>% 
          dplyr::select(-c(".y.", group1, group2)))
  
  
  return(tib_to_gt_data %>%
           gt::gt(rowname_col = "post_period_t", groupname_col = "metric") %>%
           gt::tab_header(
             title = gt::md("Metrics by Method"),
             subtitle = gt::md(data_name)
           ) %>%
           gt::tab_source_note(gt::md("Notes:")) %>%
           gt::tab_stubhead(label = "Method") %>%
           gt::fmt_number(columns=dplyr::everything(), decimals = 3))
  
  
}


