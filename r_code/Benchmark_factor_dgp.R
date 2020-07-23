library(dplyr)
library(furrr)
library(here)
library(tictoc)

source(here("r_code/panel_estimation.R"))
source(here("r_code/analysis_metrics.R"))
source(here("r_code/analysis_visualizations.R"))
source(here("r_code/treatment_effect_bootstrap.R"))
source(here("r_code/factor_DGP.R"))
cluster_v1=makeClusterPSOCK(5)
plan(cluster, workers=cluster_v1)
set.seed(1982)
n_seeds <- 50
seeds <- sample(1000:9999, size = n_seeds)

options(future.globals.maxSize= 891289600) #850mb=850*1024^2
tic("Starting DGP")
AA_data_no_sel_unformatted=future_map(.x=seeds, .f=~factor_synthetic_dgp(date_start="2010-01-01",
                                                             first_treat="2017-07-01",
                                                             date_end="2020-01-01",
                                                             treat_impact_sd = 0, 
                                                             treat_impact_mean = 0, 
                                                             rho=0.9,
                                                             rescale_y_mean = 2.5e3,
                                                             cov_overlap_scale = 0,
                                                             seed=.x))

AA_data_no_sel=future_map(AA_data_no_sel_unformatted, format_for_est)
toc()




tic("Estimating Gsynth")

gsynth_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se_est=F,
                     normalize_flag=T)
gsynth_tot_AA=future_map(gsynth_AA, compute_tot_se_jackknife, stat_in="mean")
gsynth_AA_plots=future_map(gsynth_tot_AA,
                           create_gap_ci_plot,plot_title="Gsynth Mean ToT by 
                           Post Treat Period, Jackknife CI", plot_x_lab="Post
                           Treat Period", plot_y_lab="Mean ToT")
gsynth_AA_metrics_byT=future_map(gsynth_AA,compute_avg_metric_per_t,
                              metric_str="mae", pct_eff_flag=T)
gsynth_AA_metrics=future_map(gsynth_AA,compute_avg_metric,
                              metric_str="both", pct_eff_flag=T)

gsynth_AA_coverage=future_map(gsynth_tot_AA, 
                              compute_tot_coverage)

gsynth_AA_bias=compute_jackknife_bias(gsynth_tot_AA)
gsynth_bias_plot=create_gap_ci_plot(gsynth_AA_bias, 
                                    plot_title="Gsynth Jackknife bias Post Treat Period", 
                                    plot_x_lab="Post-Treat Time",
                                    plot_y_lab="ATT Bias (%)")
gsynth_overall_metrics=compute_jackknife_metrics(gsynth_AA, horizon = 5)
toc()


tic("Estimating SCDID")
scdid_AA=future_map(AA_data_no_sel, estimate_scdid_series)
scdid_tot_AA=future_map(scdid_AA, compute_tot_se_jackknife, stat_in="mean")
scdid_AA_plots=future_map(scdid_tot_AA,
                           create_gap_ci_plot,plot_title="scdid Mean ToT by 
                           Post Treat Period, Jackknife CI", plot_x_lab="Post
                           Treat Period", plot_y_lab="Mean ToT")
scdid_AA_metrics_byT=future_map(scdid_AA,compute_avg_metric_per_t,
                                 metric_str="mae", pct_eff_flag=T)
scdid_AA_metrics=future_map(scdid_AA,compute_avg_metric,
                             metric_str="both", pct_eff_flag=T)

scdid_AA_coverage=future_map(scdid_tot_AA, 
                              compute_tot_coverage)
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

mc_AA_coverage=future_map(mc_tot_AA, 
                              compute_tot_coverage)
toc()