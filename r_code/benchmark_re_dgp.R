library(dplyr)
library(furrr)
library(here)
library(tictoc)

source(here("r_code/random_effects_dgp.R"))
source(here("r_code/panel_estimation.R"))
source(here("r_code/analysis_metrics.R"))
source(here("r_code/treatment_effect_bootstrap.R"))
cluster_v1=makeClusterPSOCK(5)
plan(cluster, workers=cluster_v1)
set.seed(1982)
n_seeds <- 10
seeds <- sample(1000:9999, size = n_seeds)


#change counter_factual to y0 and target to y1
#
tic("Starting DGP")
AA_data_no_sel=future_map(.x=seeds, .f=~gen_data_lubr(N=100, date_start="2018-01-01",
                         treat_start="2019-01-01", date_end="2020-07-04", 
                         tau_one_zero=0, T_freq = "monthly",
                         gamma=0, prop_treated=0.25, rho_y=0.3, type="random", seed=.x))
toc()

tic("Estimating Gsynth")
gsynth_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se_est=F)
gsynth_tot_AA=future_map(gsynth_AA, compute_tot_se_jackknife, counterfac_var="counter_factual")
toc()


#what would have happened in Q1 if everyone got treatment

