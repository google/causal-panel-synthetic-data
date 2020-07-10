####################################################
# Start Example:

#First, need to load in the "causal_panel_benchmark_functions.R" file and functions
####################################################
library(here)
#still an issue about where you open from -- if the file itself in r_code, then don't need the r_code part.
source(here("/r_code/causal_panel_benchmark_functions_v2.R"))


set.seed(42)

library(tictoc)

#First step: ingest data into a list of DF.

# data_list <-list.files(here("Data/"), pattern = "^synthetic_data_v.*csv$") %>% paste(here(),"/Data/",., sep="") %>%
#   lapply(.,read_csv) %>% map(., ~ (.x %>% dplyr::select(-X1))) %>%
#   map(., ~ (.x %>% mutate(period=case_when(
#     min(period)==0~period+1,
#     TRUE~period
#   ))  )) %>% map(~clean_names(.))
# #Rename the list so that each element (df) is labelled properly
# names(data_list)=list.files(pattern = "^synthetic_data_v.*csv$")


data_list <-list.files(here("Data/"), pattern = "^synthetic_data_med.*csv$") %>% paste(here(),"/Data/",., sep="") %>%
  lapply(.,read_csv) %>% map(., ~ (.x %>% dplyr::select(-X1))) %>%
  map(., ~ (.x %>% mutate(period=case_when(
    min(period)==0~period+1,
    TRUE~period
  ))  )) %>% map(~clean_names(.))
#Rename the list so that each element (df) is labelled properly
names(data_list)=list.files(pattern = "^synthetic_data_med.*csv$")


# data_list <-list.files(here("Data/"), pattern = "^synthetic_data_small.*csv$") %>% paste(here(),"/Data/",., sep="") %>%
#   lapply(.,read_csv) %>% map(., ~ (.x %>% dplyr::select(-X1))) %>% 
#   map(., ~ (.x %>% mutate(period=case_when(
#     min(period)==0~period+1,
#     TRUE~period
#   ))  )) %>% map(~clean_names(.))
# #Rename the list so that each element (df) is labelled properly
# names(data_list)=list.files(pattern = "^synthetic_data_small.*csv$")


################################################
#Compute the AB results
################################################
plan(multisession, workers=availableCores()-1)

cluster_v1=makeClusterPSOCK(5)
plan(cluster, workers=cluster_v1)
#GSYNTH IFE PROCESS
#estimate the series of impacts
tic()
gsynth_ife_AB=future_pmap(list(data_list, se_est=F), estimate_gsynth_series)
toc()
#get bootstrapped CIs, counterfactual effects


gsynth_ife_bootstrapped_pct_att_AB=future_pmap(list(gsynth_ife_AB,  counterfac_var="counter_factual" ), compute_tot_se_jackknife)
#store gap plots
gsynth_ife_bootstrapped_pct_att_AB_plots=future_pmap(list(gsynth_ife_bootstrapped_pct_att_AB,
                                                          plot_title="Gsynth Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
gsynth_ife_pct_metric_byT_AB=future_pmap(list(gsynth_ife_AB, metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
gsynth_ife_pct_metric_AB=future_pmap(list(gsynth_ife_AB, metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
gsynth_ife_pct_coverage_AB=future_map(gsynth_ife_bootstrapped_pct_att_AB, compute_tot_coverage)


Sys.sleep(19)


#SCDID PROCESS
tic()
scdid_results_AB=future_map(data_list, estimate_scdid_series,
                            id_var="entry",time_var="period", (treat_indicator="treatperiod_0"),  (outcome_var="target"),
                            (counterfac_var="counter_factual"),(pre_SDID = 0), (post_SDID = NULL),
                            (nn_SDID = NULL),  (scale_SDID = 100),  (period_SDID = 30))

toc()




#get bootstrapped CIs, counterfactual effects
scdid_bootstrapped_pct_att_AB=future_pmap(list(scdid_results_AB, counterfac_var="counter_factual"),compute_tot_se_jackknife)
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
# gsynth_ife_EM_AB=future_map(data_list, estimate_gsynth_series, se_est=F, EM_flag=T)
# #get bootstrapped CIs, counterfactual effects
# gsynth_ife_EM_bootstrapped_pct_att_AB=future_pmap(list(gsynth_ife_EM_AB,counterfac_var="counter_factual" ),compute_tot_se_jackknife)
# #store gap plots
# gsynth_ife_EM_bootstrapped_pct_att_AB_plots=future_pmap(list(gsynth_ife_EM_bootstrapped_pct_att_AB),create_gap_ci_plot)
# #compute MSE and MAE bu post period
# gsynth_ife_EM_pct_metric_byT_AB=future_pmap(list(gsynth_ife_EM_AB,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
# #compute overall MSE, MAE
# gsynth_ife_EM_pct_metric_AB=future_pmap(list(gsynth_ife_EM_AB,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
# #compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
# gsynth_ife_EM_pct_coverage_AB=future_map(gsynth_ife_EM_bootstrapped_pct_att_AB, compute_tot_coverage)


Sys.sleep(19)
#Gsynth MC PROCESS
tic()
gsynth_mc_AB=future_map(data_list, estimate_gsynth_series, se_est=F, estimator_type="mc")
toc()

#get bootstrapped CIs, counterfactual effects
gsynth_mc_bootstrapped_pct_att_AB=future_pmap(list(gsynth_mc_AB,counterfac_var="counter_factual" ),compute_tot_se_jackknife)
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

Sys.sleep(19)
tic()
causalimpact_results_AB=future_map(data_list, estimate_causalimpact_series)
toc()


#get bootstrapped CIs, counterfactual effects
causalimpact_bootstrapped_pct_att_AB=future_pmap(list(causalimpact_results_AB, counterfac_var="counter_factual" ),compute_tot_se_jackknife)
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
Sys.sleep(19)
tic()
placebo_data_list=future_map(data_list, create_placebo_df)
toc()

#GSYNTH IFE PROCESS
#estimate the series of impacts
tic()
gsynth_ife_AA=future_map(placebo_data_list, estimate_gsynth_series, se_est=F)
toc()
#get bootstrapped CIs, counterfactual effects
gsynth_ife_bootstrapped_pct_att_AA=future_map(gsynth_ife_AA,compute_tot_se_jackknife,counterfac_var="counter_factual", stat_in="median" )
#store gap plots
gsynth_ife_bootstrapped_pct_att_AA_plots=future_map(gsynth_ife_bootstrapped_pct_att_AA,create_gap_ci_plot,
                                                          plot_title="Gsynth Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT",effect_var="jackknife_median_pct_tot", cf_var="median_pct_cf_tot")
#compute MSE and MAE bu post period
gsynth_ife_pct_metric_byT_AA=future_pmap(list(gsynth_ife_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
gsynth_ife_pct_metric_AA=future_pmap(list(gsynth_ife_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
gsynth_ife_pct_coverage_AA=future_map(gsynth_ife_bootstrapped_pct_att_AA, compute_tot_coverage)


#SCDID PROCESS
Sys.sleep(19)
tic()
scdid_results_AA=future_map(placebo_data_list, estimate_scdid_series,
                             id_var="entry",time_var="period", (treat_indicator="treatperiod_0"),  (outcome_var="target"),
                             (counterfac_var="counter_factual"),(pre_SDID = 0), (post_SDID = NULL),
                             (nn_SDID = NULL),  (scale_SDID = 100),  (period_SDID = 30)) 
toc()

#get bootstrapped CIs, counterfactual effects
scdid_bootstrapped_pct_att_AA=future_map(scdid_results_AA,compute_tot_se_jackknife, counterfac_var="counter_factual", stat_in="median")
#store gap plots
scdid_bootstrapped_pct_att_AA_plots=future_map(scdid_bootstrapped_pct_att_AA,create_gap_ci_plot,
                                                     plot_title="SCDID Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT", effect_var="jackknife_median_pct_tot")
#compute MSE and MAE bu post period
scdid_pct_metric_byT_AA=future_pmap(list(scdid_results_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
scdid_pct_metric_AA=future_pmap(list(scdid_results_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
scdid_pct_coverage_AA=future_map(scdid_bootstrapped_pct_att_AA, compute_tot_coverage)


#Gsynth IFE EM PROCESS
# gsynth_ife_EM_AA=future_map(placebo_data_list, estimate_gsynth_series, se_est=F, EM_flag=T)
# 
# #get bootstrapped CIs, counterfactual effects
# gsynth_ife_EM_bootstrapped_pct_att_AA=future_pmap(list(gsynth_ife_EM_AA,counterfac_var="counter_factual" ),compute_tot_se_jackknife)
# #store gap plots
# gsynth_ife_EM_bootstrapped_pct_att_AA_plots=future_pmap(list(gsynth_ife_EM_bootstrapped_pct_att_AA,
#                                                              plot_title="Gsynth EM Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
# #compute MSE and MAE bu post period
# gsynth_ife_EM_pct_metric_byT_AA=future_pmap(list(gsynth_ife_EM_AA, metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
# #compute overall MSE, MAE
# gsynth_ife_EM_pct_metric_AA=future_pmap(list(gsynth_ife_EM_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
# #compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
# gsynth_ife_EM_pct_coverage_AA=future_map(gsynth_ife_EM_bootstrapped_pct_att_AA, compute_tot_coverage)



#Gsynth MC PROCESS
Sys.sleep(19)
tic()
gsynth_mc_AA=future_map(placebo_data_list, estimate_gsynth_series, se_est=F, estimator_type="mc")
toc()

#get bootstrapped CIs, counterfactual effects
gsynth_mc_bootstrapped_pct_att_AA=future_pmap(list(gsynth_mc_AA,counterfac_var="counter_factual" ),compute_tot_se_jackknife)
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
Sys.sleep(19)
tic()
causalimpact_results_AA=future_map((placebo_data_list), estimate_causalimpact_series)
toc()

#get bootstrapped CIs, counterfactual effects
causalimpact_bootstrapped_pct_att_AA=future_pmap(list(causalimpact_results_AA,counterfac_var="counter_factual" ),compute_tot_se_jackknife)
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
Sys.sleep(19)
est_placebo_ensemble_weights=future_pmap(list(gsynth_ife_AA, scdid_results_AA, causalimpact_results_AA, gsynth_mc_AA), ensemble_placebo_weights, pos_coef_constr=F, intercept_allowed=T)
placebo_ensemble_results=future_pmap(list(gsynth_ife_AB, scdid_results_AB, causalimpact_results_AB, gsynth_mc_AB, est_placebo_ensemble_weights), ensembled_predictor)

#get bootstrapped CIs, counterfactual effects
placebo_ensemble_bootstrapped_pct_att_AB=future_map(placebo_ensemble_results,
                                                compute_tot_se_jackknife, counterfac_var="counter_factual")
#store gap plots
placebo_ensemble_bootstrapped_pct_att_AB_plots=future_pmap(list(placebo_ensemble_bootstrapped_pct_att_AB,
                                                                plot_title="Ensemble Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
placebo_ensemble_pct_metric_byT_AB=future_pmap(list(placebo_ensemble_results,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
placebo_ensemble_pct_metric_AB=future_pmap(list(placebo_ensemble_results,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
placebo_ensemblee_pct_coverage_AB=future_map(placebo_ensemble_bootstrapped_pct_att_AB, compute_tot_coverage)




plot_tot_bias_per_t(placebo_ensemble_results[[1]], max_post_t = 1, pct_eff_flag = T)
tsfeature_pc_by_treatment_plot(causalimpact_results_AB[[1]])
tsfeature_by_treat_df(data_list[[1]])








forc_weights=future_pmap(list(gsynth_ife_AA, scdid_results_AA, causalimpact_results_AA, gsynth_mc_AA),auto_forecomb_weights_static)
forcomb_pred=future_pmap(list(gsynth_ife_AB, scdid_results_AB, causalimpact_results_AB, gsynth_mc_AB, placebo_data_list, forc_weights),
                         forecomb_predictor_static)

forcomb_bootstrapped_pct_att_AB=future_map(forcomb_pred,
                                                    compute_tot_se_jackknife, counterfac_var="counter_factual")
#store gap plots
forcomb_bootstrapped_pct_att_AB_plots=future_pmap(list(forcomb_bootstrapped_pct_att_AB,
                                                                plot_title="Ensemble Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT"),create_gap_ci_plot)
#compute MSE and MAE bu post period
forcomb_pct_metric_byT_AB=future_pmap(list(forcomb_pred,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
forcomb_pct_metric_AB=future_pmap(list(forcomb_pred,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
forcomb_pct_coverage_AB=future_map(forcomb_bootstrapped_pct_att_AB, compute_tot_coverage)





