pacman::p_load(dplyr, furrr,here,tictoc)

source(here::here("r_code/factor_dgp.R"))
source(here::here("r_code/panel_estimation.R"))
source(here::here("r_code/estimator_drivers.R"))
source(here::here("r_code/analysis_metrics.R"))
source(here::here("r_code/analysis_visualizations.R"))
source(here::here("r_code/treatment_effect_bootstrap.R"))
source(here::here("r_code/seed_metrics.R"))
source(here::here("r_code/ts_feature_analysis.R"))
source(here::here("r_code/gfoo_estimator.R"))
source(here::here("r_code/ensemble_estimators.R"))
source(here::here("r_code/placebo_creation.R"))
plan(multiprocess, workers=availableCores()-3)
set.seed(1982)
n_seeds <- 75
seeds <- sample(1000:9999, size = n_seeds)

options(future.globals.maxSize= 891289600) #850mb=850*1024^2


aa_dgp_params<-list(
  "aa_high_acf"=list(date_start="2010-01-01",
       first_treat="2017-07-01",
       date_end="2020-01-01",
       num_entries=50,
       prop_treated=0.25,
       treat_impact_sd = 0, 
       treat_impact_mean = 0, 
       rho=0.9,
       num_factors=4,
       rescale_y_mean = 2.5e3,
       cov_overlap_scale = 0,
       seed=42),
  "aa_high_acf_loading_shift"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0, 
    rho=0.9,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    cov_overlap_scale = 0,
    intercept_scale = 0,
    loading_scale = 0.75,
    seed=42
  ),
  "aa_low_acf"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0, 
    rho=0.1,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    cov_overlap_scale = 0,
    seed=42
  ),
  "aa_low_acf_int_shift"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0, 
    rho=0.1,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    intercept_scale = 2.15,
    cov_overlap_scale = 0,
    seed=42
  ),
  "aa_low_acf_load_int_shifts"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0, 
    rho=0.1,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    cov_overlap_scale = 0,
    intercept_scale = 1.9,
    loading_scale = 0.75,
    seed=42
  ),
  "aa_noisy_factors"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0, 
    rho=0.9,
    num_factors=6,
    rescale_y_mean = 2.5e3,
    cov_overlap_scale = 0,
    loading_scale = 0,
    seed=42
  ),
  
  "aa_noisy_factors_load_shift"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0, 
    rho=0.9,
    num_factors=6,
    rescale_y_mean = 2.5e3,
    loading_scale = 0.75,
    seed=42
  ),
  
  "aa_noisy_factors_lowacf"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0, 
    rho=0.2,
    num_factors=6,
    rescale_y_mean = 2.5e3,
    cov_overlap_scale = 0,
    loading_scale = 0,
    seed=42
  ),
  
  "aa_noisy_factors_load_shift_lowacf"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0, 
    rho=0.2,
    num_factors=6,
    rescale_y_mean = 2.5e3,
    loading_scale = 0.75,
    seed=42
  )
)

ab_dgp_params<-list(
  "ab_no_het"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0.1,
    treat_decay_mean=0.7,
    treat_decay_sd=0,
    rho=0.75,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    cov_overlap_scale = 0,
    seed=42 ),
  "ab_no_het_loading_shift"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0.1,
    treat_decay_mean=0.7,
    treat_decay_sd=0,
    rho=0.75,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    loading_scale = 0.75,
    seed=42 ),
  "ab_impact_het_loading_shift"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0.1, 
    treat_impact_mean = 0.1,
    treat_decay_mean=0.7,
    treat_decay_sd=0,
    rho=0.75,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    loading_scale = 0.75,
    seed=42 
  ),
  "ab_decay_het_loading_shift"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0.1,
    treat_decay_mean=0.7,
    treat_decay_sd=0.1,
    rho=0.75,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    loading_scale = 0.75,
    seed=42 
  ),
  "ab_decay_impact_het_loading_shift"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0.1,
    treat_decay_mean=0.7,
    treat_decay_sd=0.1,
    rho=0.75,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    loading_scale = 0.75,
    seed=42 
  ),
  "ab_impact_het"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0.1, 
    treat_impact_mean = 0.1,
    treat_decay_mean=0.7,
    treat_decay_sd=0,
    rho=0.75,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    loading_scale = 0,
    seed=42 
  ),
  "ab_decay_het"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0.1,
    treat_decay_mean=0.7,
    treat_decay_sd=0.1,
    rho=0.75,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    loading_scale = 0,
    seed=42 
  ),
  "ab_decay_impact_het"=list(
    date_start="2010-01-01",
    first_treat="2017-07-01",
    date_end="2020-01-01",
    num_entries=200,
    prop_treated=0.25,
    treat_impact_sd = 0, 
    treat_impact_mean = 0.1,
    treat_decay_mean=0.7,
    treat_decay_sd=0.1,
    rho=0.75,
    num_factors=4,
    rescale_y_mean = 2.5e3,
    loading_scale = 0,
    seed=42 
  )
)


use_ensemble=T
#TODO(alexdkellogg): doing double the work (calling tot on post and all periods)
list_of_dgps=c(aa_dgp_params,ab_dgp_params)
for(i in seq_len(length(list_of_dgps))){
  tic("Starting DGP")
  data_requested=do.call(factor_synthetic_dgp,list_of_dgps[[i]])
  noised_data=furrr::future_map(.x=seeds, 
                                .f=~noisify_draw(
                                  data_inp=data_requested,
                                  seed=.x,
                                  sig_y=0.2))
  formatted_data=furrr::future_map(noised_data, format_for_est)
  toc()
  
  
  tic("Estimating Gsynth")
  gsynth_est=furrr::future_map(formatted_data, estimate_gsynth_series, se=F)
  gsynth_tot=furrr::future_map(gsynth_est, compute_tot_se_jackknife, stat_in="mean")
  
  gsynth_bias=compute_jackknife_bias(gsynth_tot)
  # gsynth_bias_plot=create_gap_ci_plot(gsynth_bias, 
  #                                     plot_title="Gsynth Bias", 
  #                                     plot_x_lab="Post-Treat Time",
  #                                     plot_y_lab="ATT Bias", pct_flag = F)
  gsynth_bias_plot=plot_full_gap(gsynth_est, data_type = "raw")+
    ggplot2::labs(y="ATT Bias", title="Gsynth Bias")
  gsynth_overall_metrics=compute_jackknife_metrics(gsynth_tot)
  gsynth_tot_var=compute_tot_variance(gsynth_tot)
  gsynth_coverage=compute_tot_coverage(gsynth_tot)
  toc()
  
  tic("Estimating Gsynth, Placebo Debiased")
  gsynth_debiased_tot=furrr::future_map(formatted_data,  estimate_gsynth_tot,
                                        placebo_adjust = T, se=F)
  
  gsynth_debiased_bias=compute_jackknife_bias(gsynth_debiased_tot) %>%
    dplyr::filter(post_period_t>=0)
  gsynth_debiased_bias_plot=plot_full_gap(gsynth_debiased_tot, data_type = "tot")+
    ggplot2::labs(y="ATT Bias", title="Gsynth Debiased Bias")
  gsynth_debiased_overall_metrics=compute_jackknife_metrics(gsynth_debiased_tot)
  gsynth_debiased_tot_var=compute_tot_variance(gsynth_debiased_tot)
  gsynth_debiased_coverage=compute_tot_coverage(gsynth_debiased_tot)
  toc()

  tic("Estimating SCDID")
  scdid_est=furrr::future_map(formatted_data, estimate_scdid_series)
  scdid_tot=furrr::future_map(scdid_est, compute_tot_se_jackknife, stat_in="mean")

  scdid_bias=compute_jackknife_bias(scdid_tot)
  # scdid_bias_plot=create_gap_ci_plot(scdid_bias,
  #                                       plot_title="SCDID Bias",
  #                                       plot_x_lab="Post-Treat Time",
  #                                       plot_y_lab="ATT Bias", pct_flag = F)
  scdid_bias_plot=plot_full_gap(scdid_est, data_type = "raw")+
    ggplot2::labs(y="ATT Bias", title="SCDID Bias")
  scdid_overall_metrics=compute_jackknife_metrics(scdid_tot)
  scdid_tot_var=compute_tot_variance(scdid_tot)
  scdid_coverage=compute_tot_coverage(scdid_tot)
  toc()
  
  tic("Estimating SCDID unconstrained")
  scdid_uncon_est=furrr::future_map(formatted_data, estimate_scdid_series, constrained=F)
  scdid_uncon_tot=furrr::future_map(scdid_uncon_est, compute_tot_se_jackknife, stat_in="mean")
  
  scdid_uncon_bias=compute_jackknife_bias(scdid_uncon_tot)
  # scdid_bias_plot=create_gap_ci_plot(scdid_bias,
  #                                       plot_title="SCDID Bias",
  #                                       plot_x_lab="Post-Treat Time",
  #                                       plot_y_lab="ATT Bias", pct_flag = F)
  scdid_uncon_bias_plot=plot_full_gap(scdid_uncon_est,data_type = "raw")+
    ggplot2::labs(y="ATT Bias", title="SCDID Unconstrained Bias")
  scdid_uncon_overall_metrics=compute_jackknife_metrics(scdid_uncon_tot)
  scdid_uncon_tot_var=compute_tot_variance(scdid_uncon_tot)
  scdid_uncon_coverage=compute_tot_coverage(scdid_uncon_tot)
  toc()
  
  tic("Estimating SCM unconstrained")
  scm_est=furrr::future_map(formatted_data, estimate_scm_series)
  scm_tot=furrr::future_map(scm_est, compute_tot_se_jackknife, stat_in="mean")
  
  scm_bias=compute_jackknife_bias(scm_tot)
  # scdid_bias_plot=create_gap_ci_plot(scdid_bias,
  #                                       plot_title="SCDID Bias",
  #                                       plot_x_lab="Post-Treat Time",
  #                                       plot_y_lab="ATT Bias", pct_flag = F)
  scm_bias_plot=plot_full_gap(scm_est,data_type = "raw")+
    ggplot2::labs(y="ATT Bias", title="SCM Bias")
  scm_overall_metrics=compute_jackknife_metrics(scm_tot)
  scm_tot_var=compute_tot_variance(scm_tot)
  scm_coverage=compute_tot_coverage(scm_tot)
  toc()
  
  tic("Estimating Gfoo")
  gfoo_est=furrr::future_map(formatted_data, estimate_gfoo_series, season_freq=4)
  gfoo_tot=furrr::future_map(gfoo_est, compute_tot_se_jackknife, stat_in="mean")
  
  gfoo_bias=compute_jackknife_bias(gfoo_tot)
  # scdid_bias_plot=create_gap_ci_plot(scdid_bias,
  #                                       plot_title="SCDID Bias",
  #                                       plot_x_lab="Post-Treat Time",
  #                                       plot_y_lab="ATT Bias", pct_flag = F)
  gfoo_bias_plot=plot_full_gap(gfoo_est,data_type = "raw")+
    ggplot2::labs(y="ATT Bias", title="Gfoo Bias")
  gfoo_overall_metrics=compute_jackknife_metrics(gfoo_tot)
  gfoo_tot_var=compute_tot_variance(gfoo_tot)
  gfoo_coverage=compute_tot_coverage(gfoo_tot)
  toc()
  

  tic("Estimating MC")
  mc_est=furrr::future_map(formatted_data, estimate_gsynth_series, se=F,
                   estimator="mc")
  mc_tot=furrr::future_map(mc_est, compute_tot_se_jackknife, stat_in="mean")

  mc_bias=compute_jackknife_bias(mc_tot)
  # mc_bias_plot=create_gap_ci_plot(mc_bias,
  #                                    plot_title="MC Bias",
  #                                    plot_x_lab="Post-Treat Time",
  #                                    plot_y_lab="ATT Bias",  pct_flag = F)
  mc_bias_plot=plot_full_gap(mc_est,data_type = "raw")+
    ggplot2::labs(y="ATT Bias", title="MC Bias")
  mc_overall_metrics=compute_jackknife_metrics(mc_tot)
  mc_tot_var=compute_tot_variance(mc_tot)
  mc_coverage=compute_tot_coverage(mc_tot)
  toc()

  tic("Estimating Causal Impact")
  causalimpact_est=furrr::future_map(formatted_data, estimate_causalimpact_series)
  causalimpact_tot=furrr::future_map(causalimpact_est, compute_tot_se_jackknife, stat_in="mean")

  causalimpact_bias=compute_jackknife_bias(causalimpact_tot)
  # causalimpact_bias_plot=create_gap_ci_plot(causalimpact_bias,
  #                                              plot_title="Causal Impact Bias",
  #                                              plot_x_lab="Post-Treat Time",
  #                                              plot_y_lab="ATT Bias",  pct_flag = F)
  causalimpact_bias_plot=plot_full_gap(causalimpact_est,data_type = "raw")+
    ggplot2::labs(y="ATT Bias", title="Causal Impact Bias")
  causalimpact_overall_metrics=compute_jackknife_metrics(causalimpact_tot)
  causalimpact_tot_var=compute_tot_variance(causalimpact_tot)
  causalimpact_coverage=compute_tot_coverage(causalimpact_tot)
  toc()
  
  tic("Estimating BART")
  bart_est=furrr::future_map(formatted_data, estimate_bart_series)
  bart_tot=furrr::future_map(bart_est, compute_tot_se_jackknife, stat_in="mean")
  
  bart_bias=compute_jackknife_bias(bart_tot)

  bart_bias_plot=plot_full_gap(bart_est,data_type = "raw")+
    ggplot2::labs(y="ATT Bias", title="BART Bias")
  bart_overall_metrics=compute_jackknife_metrics(bart_tot)
  bart_tot_var=compute_tot_variance(bart_tot)
  bart_coverage=compute_tot_coverage(bart_tot)
  toc()
  
  if(use_ensemble){
    #For ensemble -- create an AA version of the data, estimate weights, apply
    tic("Estimating Ensemble")
    
    
    ensemble_est=furrr::future_map2(.x=formatted_data, 
                       .y=preprocess_ensemble_input(list(gsynth_est, gfoo_est,
                                                         scdid_uncon_est)),
                       .f=~estimate_ensemble(
                         method_names=c("gsynth", "gfoo", "scdid_uncon"),
                         true_data=.x, pred_list=.y))
    
    ensemble_tot=furrr::future_map(ensemble_est, compute_tot_se_jackknife, stat_in="mean")
    ensemble_bias=compute_jackknife_bias(ensemble_tot)
    # ensemble_bias_plot=create_gap_ci_plot(ensemble_bias,
    #                                           plot_title="Ensemble Bias",
    #                                           plot_x_lab="Post-Treat Time",
    #                                           plot_y_lab="ATT Bias",  pct_flag = F)
    ensemble_bias_plot=plot_full_gap(ensemble_est)+
      ggplot2::labs(y="ATT Bias", title="Ensemble Bias")
    ensemble_overall_metrics=compute_jackknife_metrics(ensemble_tot)
    ensemble_tot_var=compute_tot_variance(ensemble_tot)
    ensemble_coverage=compute_tot_coverage(ensemble_tot)
    toc()  
  }
  #save.image(here::here("Data", "Variations", paste(names(list_of_dgps)[i],".RData",sep = "")))
}





