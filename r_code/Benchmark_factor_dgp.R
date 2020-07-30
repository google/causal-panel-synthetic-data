pacman::p_load(dplyr, furrr,here,tictoc)

source(here::here("r_code/factor_dgp.R"))
source(here::here("r_code/panel_estimation.R"))
source(here::here("r_code/analysis_metrics.R"))
source(here::here("r_code/analysis_visualizations.R"))
source(here::here("r_code/treatment_effect_bootstrap.R"))
source(here::here("r_code/seed_metrics.R"))
plan(multiprocess, workers=availableCores()-1)
set.seed(1982)
n_seeds <- 50
seeds <- sample(1000:9999, size = n_seeds)

options(future.globals.maxSize= 891289600) #850mb=850*1024^2


aa_dgp_params<-list(
  "aa_high_acf"=list(date_start="2010-01-01",
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
  "aa_low_acf_sel_covariate_shift"=list(
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
    selection="observables",
    cov_overlap_scale = 0.75,
    loading_scale = 0,
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


list_of_dgps=ab_dgp_params
for(i in seq_len(length(list_of_dgps))){
  tic("Starting DGP")
  data_requested=do.call(factor_synthetic_dgp,list_of_dgps[[i]])
  noised_data=furrr::future_map(.x=seeds, 
                                .f=~noisify_draw(
                                  data_inp=data_requested,
                                  seed=.x))
  formatted_data=future_map(noised_data, format_for_est)
  toc()
  
  
  tic("Estimating Gsynth")
  gsynth_est=future_map(formatted_data, estimate_gsynth_series, se=F)
  gsynth_tot=future_map(gsynth_est, compute_tot_se_jackknife, stat_in="mean")
  
  gsynth_bias=compute_jackknife_bias(gsynth_tot)
  gsynth_bias_plot=create_gap_ci_plot(gsynth_bias, 
                                      plot_title="Gsynth Bias", 
                                      plot_x_lab="Post-Treat Time",
                                      plot_y_lab="ATT Bias", pct_flag = F)
  gsynth_overall_metrics=compute_jackknife_metrics(gsynth_est)
  gsynth_tot_var=compute_tot_variance(gsynth_tot)
  gsynth_coverage=compute_tot_coverage(gsynth_tot)
  toc()

  tic("Estimating SCDID")
  scdid_est=future_map(formatted_data, estimate_scdid_series)
  scdid_tot=future_map(scdid_est, compute_tot_se_jackknife, stat_in="mean")

  scdid_bias=compute_jackknife_bias(scdid_tot)
  scdid_bias_plot=create_gap_ci_plot(scdid_bias,
                                        plot_title="SCDID Bias",
                                        plot_x_lab="Post-Treat Time",
                                        plot_y_lab="ATT Bias", pct_flag = F)
  scdid_overall_metrics=compute_jackknife_metrics(scdid_est)
  scdid_tot_var=compute_tot_variance(scdid_tot)
  scdid_coverage=compute_tot_coverage(scdid_tot)
  toc()

  tic("Estimating MC")
  mc_est=future_map(formatted_data, estimate_gsynth_series, se=F,
                   estimator="mc")
  mc_tot=future_map(mc_est, compute_tot_se_jackknife, stat_in="mean")

  mc_bias=compute_jackknife_bias(mc_tot)
  mc_bias_plot=create_gap_ci_plot(mc_bias,
                                     plot_title="MC Bias",
                                     plot_x_lab="Post-Treat Time",
                                     plot_y_lab="ATT Bias",  pct_flag = F)
  mc_overall_metrics=compute_jackknife_metrics(mc_est)
  mc_tot_var=compute_tot_variance(mc_tot)
  mc_coverage=compute_tot_coverage(mc_tot)
  toc()

  tic("Estimating Causal Impact")
  causalimpact_est=future_map(formatted_data, estimate_causalimpact_series)
  causalimpact_tot=future_map(causalimpact_est, compute_tot_se_jackknife, stat_in="mean")

  causalimpact_bias=compute_jackknife_bias(causalimpact_tot)
  causalimpact_bias_plot=create_gap_ci_plot(causalimpact_bias,
                                               plot_title="Causal Impact Bias",
                                               plot_x_lab="Post-Treat Time",
                                               plot_y_lab="ATT Bias",  pct_flag = F)
  causalimpact_overall_metrics=compute_jackknife_metrics(causalimpact_est)
  causalimpact_tot_var=compute_tot_variance(causalimpact_tot)
  causalimpact_coverage=compute_tot_coverage(causalimpact_tot)
  toc()
  
  save.image(here::here(paste("Data/",names(list_of_dgps)[i],".RData",sep = "")))
  #save.image(here::here(paste("Data/",glue::glue("Data{i}.RData"),sep = "")))
}





