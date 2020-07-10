####################################################
# Start Example:

#First, need to load in the "causal_panel_benchmark_functions.R" file and functions
####################################################
source("causal_panel_benchmark_functions_v2.R")


set.seed(42)

library(tictoc)
library(here)
#First step: ingest data into a list of DF.
#
# data_list <-list.files(here("Data/"), pattern = "^synthetic.*csv$") %>% paste(here(),"/Data/",., sep="") %>%
#   lapply(.,read_csv) %>% map(., ~ (.x %>% dplyr::select(-X1))) %>% 
#   map(., ~ (.x %>% mutate(period=case_when(
#     min(period)==0~period+1,
#     TRUE~period
#   ))  )) %>% map(~clean_names(.))
# #Rename the list so that each element (df) is labelled properly
# names(data_list)=list.files(pattern = "^synthetic.*csv$")


data_list <-list.files(here("Data/"), pattern = "^synthetic_data_small.*csv$") %>% paste(here(),"/Data/",., sep="") %>%
  lapply(.,read_csv) %>% map(., ~ (.x %>% dplyr::select(-X1))) %>% 
  map(., ~ (.x %>% mutate(period=case_when(
    min(period)==0~period+1,
    TRUE~period
  ))  )) %>% map(~clean_names(.))
#Rename the list so that each element (df) is labelled properly
names(data_list)=list.files(pattern = "^synthetic_data_small.*csv$")


################################################
#Compute the AB results
################################################
plan(multiprocess, workers=availableCores()-1)
#GSYNTH IFE PROCESS
#estimate the series of impacts
gsynth_ife_AB=future_pmap(list(data_list, se_est=F), estimate_gsynth_series)
#get bootstrapped CIs, counterfactual effects


gsynth_ife_bootstrapped_pct_att_AB=future_pmap(list(gsynth_ife_AB,  boot_var="pct.effect", counterfac_eff="cf_pct.effect" ),compute_tot_se_jackknife)
#store gap plots
gsynth_ife_bootstrapped_pct_att_AB_plots=future_pmap(list(gsynth_ife_bootstrapped_pct_att_AB,
                                                          plot_title="Gsynth Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
gsynth_ife_pct_metric_byT_AB=future_pmap(list(gsynth_ife_AB, metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
gsynth_ife_pct_metric_AB=future_pmap(list(gsynth_ife_AB, metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
gsynth_ife_pct_coverage_AB=future_map(gsynth_ife_bootstrapped_pct_att_AB, compute_tot_coverage)



#SCDID PROCESS
scdid_results_AB=future_pmap(data_list, estimate_scdid_series,
                             id_var="entry",time_var="period", (treat_indicator="treatperiod_0"),  (outcome_var="target"),
                             (counterfac_var="counter_factual"),(pre_SDID = 0), (post_SDID = NULL),
                            (nn_SDID = NULL),  (scale_SDID = 100),  (period_SDID = 30)) 
  


#get bootstrapped CIs, counterfactual effects
scdid_bootstrapped_pct_att_AB=future_pmap(list(scdid_results_AB, boot_var="pct.effect", counterfac_eff="cf_pct.effect"),compute_tot_se_jackknife)
#store gap plots
scdid_bootstrapped_pct_att_AB_plots=future_pmap(list(scdid_bootstrapped_pct_att_AB,
                                                     plot_title="SCDID Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
scdid_pct_metric_byT_AB=future_pmap(list(scdid_results_AB, metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
scdid_pct_metric_AB=future_pmap(list(scdid_results_AB, metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
scdid_pct_coverage_AB=future_map(scdid_bootstrapped_pct_att_AB, compute_tot_coverage)



#Gsynth IFE EM PROCESS
gsynth_ife_EM_AB=future_pmap(list(data_list, se_est=F, EM_flag=T), estimate_gsynth_series)
#get bootstrapped CIs, counterfactual effects
gsynth_ife_EM_bootstrapped_pct_att_AB=future_pmap(list(gsynth_ife_EM_AB,boot_var="pct.effect", counterfac_eff="cf_pct.effect" ),compute_tot_se_jackknife)
#store gap plots
gsynth_ife_EM_bootstrapped_pct_att_AB_plots=future_pmap(list(gsynth_ife_EM_bootstrapped_pct_att_AB),create_gap_ci_plot)
#compute MSE and MAE bu post period
gsynth_ife_EM_pct_metric_byT_AB=future_pmap(list(gsynth_ife_EM_AB,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
gsynth_ife_EM_pct_metric_AB=future_pmap(list(gsynth_ife_EM_AB,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
gsynth_ife_EM_pct_coverage_AB=future_map(gsynth_ife_EM_bootstrapped_pct_att_AB, compute_tot_coverage)



#Gsynth MC PROCESS
gsynth_mc_AB=future_pmap(list(data_list, se_est=F, estimator_type="mc"), estimate_gsynth_series)

#get bootstrapped CIs, counterfactual effects
gsynth_mc_bootstrapped_pct_att_AB=future_pmap(list(gsynth_mc_AB,boot_var="pct.effect", counterfac_eff="cf_pct.effect" ),compute_tot_se_jackknife)
#store gap plots
gsynth_mc_bootstrapped_pct_att_AB_plots=future_pmap(list(gsynth_mc_bootstrapped_pct_att_AB),create_gap_ci_plot)
#compute MSE and MAE bu post period
gsynth_mc_pct_metric_byT_AB=future_pmap(list(gsynth_mc_AB,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
gsynth_mc_pct_metric_AB=future_pmap(list(gsynth_mc_AB,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
gsynth_mc_pct_coverage_AB=future_map(gsynth_mc_bootstrapped_pct_att_AB, compute_tot_coverage)



#Causal Impact PRocess
#Gsynth MC PROCESS
causalimpact_results_AB=future_pmap(list(data_list), estimate_causalimpact_series)

#get bootstrapped CIs, counterfactual effects
causalimpact_bootstrapped_pct_att_AB=future_pmap(list(causalimpact_results_AB, boot_var="pct.effect", counterfac_eff="cf_pct.effect" ),compute_tot_se_jackknife)
#store gap plots
causalimpact_bootstrapped_pct_att_AB_plots=future_pmap(list(causalimpact_bootstrapped_pct_att_AB,
                                                            plot_title="Causal Impact Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
causalimpact_pct_metric_byT_AB=future_pmap(list(causalimpact_results_AB,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
causalimpact_pct_metric_AB=future_pmap(list(causalimpact_results_AB,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
causalimpact_pct_coverage_AB=future_map(causalimpact_bootstrapped_pct_att_AB, compute_tot_coverage)






################################################
#Compute the AA results
################################################

placebo_data_list=future_pmap(list(data_list), create_placebo_df)

#GSYNTH IFE PROCESS
#estimate the series of impacts
gsynth_ife_AA=future_pmap(list(placebo_data_list, se_est=F), estimate_gsynth_series)
#get bootstrapped CIs, counterfactual effects
gsynth_ife_bootstrapped_pct_att_AA=future_pmap(list(gsynth_ife_AA,boot_var="pct.effect", counterfac_eff="cf_pct.effect" ),compute_tot_se_jackknife)
#store gap plots
gsynth_ife_bootstrapped_pct_att_AA_plots=future_pmap(list(gsynth_ife_bootstrapped_pct_att_AA,
                                                          plot_title="Gsynth Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
gsynth_ife_pct_metric_byT_AA=future_pmap(list(gsynth_ife_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
gsynth_ife_pct_metric_AA=future_pmap(list(gsynth_ife_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
gsynth_ife_pct_coverage_AA=future_map(gsynth_ife_bootstrapped_pct_att_AA, compute_tot_coverage)


#SCDID PROCESS
scdid_results_AA=future_pmap(placebo_data_list, estimate_scdid_series,
                             id_var="entry",time_var="period", (treat_indicator="treatperiod_0"),  (outcome_var="target"),
                             (counterfac_var="counter_factual"),(pre_SDID = 0), (post_SDID = NULL),
                             (nn_SDID = NULL),  (scale_SDID = 100),  (period_SDID = 30)) 


#get bootstrapped CIs, counterfactual effects
scdid_bootstrapped_pct_att_AA=future_pmap(list(scdid_results_AA,boot_var="pct.effect", counterfac_eff="cf_pct.effect" ),compute_tot_se_jackknife)
#store gap plots
scdid_bootstrapped_pct_att_AA_plots=future_pmap(list(scdid_bootstrapped_pct_att_AA,
                                                     plot_title="SCDID Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
scdid_pct_metric_byT_AA=future_pmap(list(scdid_results_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
scdid_pct_metric_AA=future_pmap(list(scdid_results_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
scdid_pct_coverage_AA=future_map(scdid_bootstrapped_pct_att_AA, compute_tot_coverage)


#Gsynth IFE EM PROCESS
gsynth_ife_EM_AA=future_pmap(list(placebo_data_list, se_est=F, EM_flag=T), estimate_gsynth_series)

#get bootstrapped CIs, counterfactual effects
gsynth_ife_EM_bootstrapped_pct_att_AA=future_pmap(list(gsynth_ife_EM_AA,boot_var="pct.effect", counterfac_eff="cf_pct.effect" ),compute_tot_se_jackknife)
#store gap plots
gsynth_ife_EM_bootstrapped_pct_att_AA_plots=future_pmap(list(gsynth_ife_EM_bootstrapped_pct_att_AA,
                                                             plot_title="Gsynth EM Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
gsynth_ife_EM_pct_metric_byT_AA=future_pmap(list(gsynth_ife_EM_AA, metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
gsynth_ife_EM_pct_metric_AA=future_pmap(list(gsynth_ife_EM_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
gsynth_ife_EM_pct_coverage_AA=future_map(gsynth_ife_EM_bootstrapped_pct_att_AA, compute_tot_coverage)



#Gsynth MC PROCESS
gsynth_mc_AA=future_pmap(list(placebo_data_list, se_est=F, estimator_type="mc"), estimate_gsynth_series)

#get bootstrapped CIs, counterfactual effects
gsynth_mc_bootstrapped_pct_att_AA=future_pmap(list(gsynth_mc_AA,boot_var="pct.effect", counterfac_eff="cf_pct.effect" ),compute_tot_se_jackknife)
#store gap plots
gsynth_mc_bootstrapped_pct_att_AA_plots=future_pmap(list(gsynth_mc_bootstrapped_pct_att_AA,
                                                         plot_title="Gsynth MC Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
gsynth_mc_pct_metric_byT_AA=future_pmap(list(gsynth_mc_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
gsynth_mc_pct_metric_AA=future_pmap(list(gsynth_mc_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
gsynth_mc_pct_coverage_AA=future_map(gsynth_mc_bootstrapped_pct_att_AA, compute_tot_coverage)





#Causal Impact PRocess
#Gsynth MC PROCESS
causalimpact_results_AA=future_pmap(list(placebo_data_list), estimate_causalimpact_series)

#get bootstrapped CIs, counterfactual effects
causalimpact_bootstrapped_pct_att_AA=future_pmap(list(causalimpact_results_AA,boot_var="pct.effect", counterfac_eff="cf_pct.effect" ),compute_tot_se_jackknife)
#store gap plots
causalimpact_bootstrapped_pct_att_AA_plots=future_pmap(list(causalimpact_bootstrapped_pct_att_AA,
                                                            plot_title="Causal Impact Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
causalimpact_pct_metric_byT_AA=future_pmap(list(causalimpact_results_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
causalimpact_pct_metric_AA=future_pmap(list(causalimpact_results_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
causalimpact_pct_coverage_AA=future_map(causalimpact_bootstrapped_pct_att_AA, compute_tot_coverage)






#Placebo Ensemble Method

est_placebo_ensemble_weights=future_pmap(list(gsynth_ife_AA, scdid_results_AA, causalimpact_results_AA, pos_coef_constr=F), ensemble_placebo_weights)
placebo_ensemble_results=future_pmap(list(gsynth_ife_AB, scdid_results_AB, causalimpact_results_AB, est_placebo_ensemble_weights), ensembled_predictor)

#get bootstrapped CIs, counterfactual effects
placebo_ensemble_bootstrapped_pct_att_AB=future_pmap(list(placebo_ensemble_results,boot_var="pct.effect", counterfac_eff="cf_pct.effect" ),compute_tot_se_jackknife)
#store gap plots
placebo_ensemble_bootstrapped_pct_att_AB_plots=future_pmap(list(placebo_ensemble_bootstrapped_pct_att_AB,
                                                                plot_title="Ensemble Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
placebo_ensemble_pct_metric_byT_AB=future_pmap(list(placebo_ensemble_results,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
placebo_ensemble_pct_metric_AB=future_pmap(list(placebo_ensemble_results,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
placebo_ensemblee_pct_coverage_AB=future_map(placebo_ensemble_bootstrapped_pct_att_AB, compute_tot_coverage)



plot_tot_bias_per_t(placebo_ensemble_results[[1]], max_post_t = 1)

