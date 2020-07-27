pacman::p_load(dplyr, furrr,here,tictoc)

source(here("r_code/factor_dgp.R"))
source(here("r_code/panel_estimation.R"))
source(here("r_code/analysis_metrics.R"))
source(here("r_code/analysis_visualizations.R"))
source(here("r_code/treatment_effect_bootstrap.R"))
source(here("r_code/seed_metrics.R"))
plan(multiprocess, workers=availableCores()-1)
set.seed(1982)
n_seeds <- 2
seeds <- sample(1000:9999, size = n_seeds)

options(future.globals.maxSize= 891289600) #850mb=850*1024^2

tic("Starting DGP")
AA_data_no_sel_base=factor_synthetic_dgp(date_start="2010-01-01",
                                         first_treat="2017-07-01",
                                         date_end="2020-01-01",
                                         treat_impact_sd = 0, 
                                         treat_impact_mean = 0, 
                                         rho=0.9,
                                         rescale_y_mean = 2.5e3,
                                         cov_overlap_scale = 0,
                                         seed=42)

AA_data_no_sel_unformatted=furrr::future_map(.x=seeds, 
                                             .f=~noisify_draw(data_inp=
                                                           AA_data_no_sel_base,
                                                         seed=.x))

AA_data_no_sel=future_map(AA_data_no_sel_unformatted, format_for_est)
toc()

tic("Starting DGP, Selection")
AA_data_sel_base=factor_synthetic_dgp(date_start="2010-01-01",
                                             first_treat="2017-07-01",
                                             date_end="2020-01-01",
                                             treat_impact_sd = 0, 
                                             treat_impact_mean = 0, 
                                             rho=0.9,
                                             rescale_y_mean = 2.5e3,
                                             cov_overlap_scale = 0,
                                             intercept_scale = 0.75,
                                             seed=42)
                                                                      
AA_data_sel_unformatted=furrr::future_map(.x=seeds, 
                                             .f=~noisify_draw(data_inp=
                                                                AA_data_sel_base,
                                                              seed=.x))

AA_data_sel=future_map(AA_data_sel_unformatted, format_for_est)
toc()



tic("Estimating Gsynth, Placebo, No Selection")
gsynth_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se_est=F)
gsynth_tot_AA=future_map(gsynth_AA, compute_tot_se_jackknife, stat_in="mean")

gsynth_AA_bias=compute_jackknife_bias(gsynth_tot_AA)
gsynth_bias_plot=create_gap_ci_plot(gsynth_AA_bias, 
                                    plot_title="Gsynth Jackknife bias Post Treat Period", 
                                    plot_x_lab="Post-Treat Time",
                                    plot_y_lab="ATT Bias", pct_flag = F)
gsynth_overall_metrics=compute_jackknife_metrics(gsynth_AA, pct_eff_flag = F)
gsynth_AA_tot_var=compute_tot_variance(gsynth_tot_AA)
gsynth_AA_cov=compute_tot_coverage(gsynth_tot_AA)
toc()



tic("Estimating Gsynth, Placebo, Selection")
gsynth_AA_sel=future_map(AA_data_sel, estimate_gsynth_series, se_est=F)
gsynth_tot_AA_sel=future_map(gsynth_AA_sel, compute_tot_se_jackknife,
                             stat_in="mean")

gsynth_AA_sel_bias=compute_jackknife_bias(gsynth_tot_AA_sel)
gsynth_sel_bias_plot=create_gap_ci_plot(gsynth_AA_sel_bias, 
                                    plot_title="Gsynth Jackknife bias Post Treat Period", 
                                    plot_x_lab="Post-Treat Time",
                                    plot_y_lab="ATT Bias", pct_flag = F)
gsynth_sel_overall_metrics=compute_jackknife_metrics(gsynth_AA_sel, pct_eff_flag = F)
gsynth_sel_AA_tot_var=compute_tot_variance(gsynth_tot_AA_sel)
gsynth_AA_sel_cov=compute_tot_coverage(gsynth_tot_AA_sel )
toc()




tic("Estimating SCDID")
scdid_AA=future_map(AA_data_no_sel, estimate_scdid_series)
scdid_tot_AA=future_map(scdid_AA, compute_tot_se_jackknife, stat_in="mean")

scdid_AA_bias=compute_jackknife_bias(scdid_tot_AA)
scdid_AA_bias_plot=create_gap_ci_plot(scdid_AA_bias, 
                                        plot_title="Gsynth Jackknife bias Post Treat Period", 
                                        plot_x_lab="Post-Treat Time",
                                        plot_y_lab="ATT Bias", pct_flag = F)
scdid_AA_overall_metrics=compute_jackknife_metrics(scdid_AA, pct_flag = F)
scdid_AA_tot_var=compute_tot_variance(scdid_tot_AA)
scdid_AA_cov=compute_tot_coverage(scdid_tot_AA)
toc()

tic("Estimating SCDID Sel")
scdid_AA_sel=future_map(AA_data_sel, estimate_scdid_series)
scdid_tot_AA_sel=future_map(scdid_AA_sel, compute_tot_se_jackknife, stat_in="mean")

scdid_AA_sel_bias=compute_jackknife_bias(scdid_tot_AA_sel)
scdid_AA_sel_bias_plot=create_gap_ci_plot(scdid_AA_sel_bias, 
                                      plot_title="Gsynth Jackknife bias Post Treat Period", 
                                      plot_x_lab="Post-Treat Time",
                                      plot_y_lab="ATT Bias")
scdid_AA_sel_overall_metrics=compute_jackknife_metrics(scdid_AA_sel, pct_flag = F)
scdid_AA_sel_tot_var=compute_tot_variance(scdid_tot_AA_sel)
scdid_AA_sel_cov=compute_tot_coverage(scdid_tot_AA_sel)
toc()

tic("Estimating MC")

mc_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se_est=F,
                 estimator_type="mc")
mc_tot_AA=future_map(mc_AA, compute_tot_se_jackknife, stat_in="mean")
mc_AA_plots=future_map(mc_tot_AA,
                       create_gap_ci_plot,plot_title="mc Mean ToT by 
                           Post Treat Period, Jackknife CI", plot_x_lab="Post
                           Treat Period", plot_y_lab="Mean ToT")
mc_AA_metrics_byT=future_map(mc_AA,compute_avg_metric_per_t,
                             metric_str="mae", pct_eff_flag=T)
mc_AA_metrics=future_map(mc_AA,compute_avg_metric,
                         metric_str="both", pct_eff_flag=T)

mc_AA_bias=compute_jackknife_bias(mc_tot_AA)
mc_AA_bias_plot=create_gap_ci_plot(mc_AA_bias, 
                                      plot_title="Gsynth Jackknife bias Post Treat Period", 
                                      plot_x_lab="Post-Treat Time",
                                      plot_y_lab="ATT Bias")
mc_AA_overall_metrics=compute_jackknife_metrics(mc_AA, pct_flag = F)
mc_AA_tot_var=compute_tot_variance(mc_tot_AA)
mc_AA_cov=compute_tot_coverage(mc_tot_AA)
toc()




tic("Estimating MC sel")

mc_AA_sel=future_map(AA_data_sel, estimate_gsynth_series, se_est=F,
                 estimator_type="mc")
mc_tot_AA_sel=future_map(mc_AA_sel, compute_tot_se_jackknife, stat_in="mean")
mc_AA_sel_plots=future_map(mc_tot_AA_sel,
                       create_gap_ci_plot,plot_title="mc Mean ToT by 
                           Post Treat Period, Jackknife CI", plot_x_lab="Post
                           Treat Period", plot_y_lab="Mean ToT")
mc_AA_sel_metrics_byT=future_map(mc_AA_sel,compute_avg_metric_per_t,
                             metric_str="mae", pct_eff_flag=T)
mc_AA_sel_metrics=future_map(mc_AA_sel,compute_avg_metric,
                         metric_str="both", pct_eff_flag=T)

mc_AA_sel_bias=compute_jackknife_bias(mc_tot_AA_sel)
mc_AA_sel_bias_plot=create_gap_ci_plot(mc_AA_sel_bias, 
                                          plot_title="Gsynth Jackknife bias Post Treat Period", 
                                          plot_x_lab="Post-Treat Time",
                                          plot_y_lab="ATT Bias")
mc_AA_sel_overall_metrics=compute_jackknife_metrics(mc_AA_sel, pct_flag = F)
mc_AA_sel_tot_var=compute_tot_variance(mc_tot_AA_sel)
mc_AA_sel_cov=compute_tot_coverage(mc_tot_AA_sel)
toc()

