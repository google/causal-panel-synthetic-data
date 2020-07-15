library(CausalImpact)
library(Matrix)
library(tsfeatures)
library(tsibble)
library(ggfortify)
library(gsynth)
library(augsynth)
library(tidyverse)
library(panelView)
library(synthdid)
library(resample)
library(mvtnorm)
library(janitor)
library("qpcR")
library(dplyr)
library(furrr)
library(gridExtra)
library(quadprog)
library(rstatix)
library(ForecastComb)

# lots of edits:
# Fix: Include additional weight selection methods? How about adopting the SCDID code and doing double SCDID?
ensemble_placebo_weights <- function(method1_estimated_df, method2_estimated_df, method3_estimated_df, method4_estimated_df,
                                     time_var = "period", id_var = "entry", treat_time_var = "Treatment_Period",
                                     pred_var = "point.pred", counterfac_var = "counter_factual",
                                     constrained = T, intercept_allowed = T) {
  
  # Estimates the weights for an ensemble using linear regression on a placebo set
  # This method is, as of now, highly specific -- intended for a placebo data set
  # with fake treatment (but zero treatment effect). Given this structure, we have fit several SCM methods
  # on this placebo set to estimate the TE (hopefully close to zero)
  # This method ignores any pre-treat information and outputs weights on the post-treatment point predictions
  # from each method that yields the closest value to the truth. Thus, we must know the truth to compute these.
  
  
  # Args
  # method1_estimated_df: long-form dataframe of the estimated and counterfactual point predictions for a given method
  # method2_estimated_df, method3_estimated_df distinct methods of the same form.
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_time_var: string variable for the column indicating when that unit was treated
  # pred_var: string name of the point estimate column in the methodX dataframe
  # counterfac_var: name of var for the counterfactual value for the time series, if available (otherwise, NULL)
  # pos_coef_constr: Boolean flag for whether the weights should be non-negative
  # sum_to_one: Boolean flag for whether the coefficients (without intercept) sum to 1
  
  # Output
  # weights, as a 4x1 (num methods X 1) vector, from an unconstrained linear reg
  
  combined_methods_df <- method1_estimated_df %>%
    dplyr::select(
      !!as.name(time_var), !!as.name(id_var), !!as.name(treat_time_var),
      !!as.name(pred_var), !!as.name(counterfac_var)
    ) %>%
    rename(m1_pred = !!as.name(pred_var)) %>%
    left_join(
      method2_estimated_df %>%
        dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(pred_var)) %>%
        rename(m2_pred = !!as.name(pred_var)),
      by = c(id_var, time_var)
    ) %>%
    left_join(
      method3_estimated_df %>%
        dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(pred_var)) %>%
        rename(m3_pred = !!as.name(pred_var)),
      by = c(id_var, time_var)
    ) %>%
    left_join(
      method4_estimated_df %>%
        dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(pred_var)) %>%
        rename(m4_pred = !!as.name(pred_var)),
      by = c(id_var, time_var)
    )
  
  # Find the weights of the three methods that best fit the data in the post period
  # For now, we do this with simple linear regression and no constraints
  post_treat_combined_df <- combined_methods_df %>% filter(!!as.name(time_var) >= !!as.name(treat_time_var))
  
  x=as.matrix(cbind(
    combined_methods_df %>% pull(m1_pred),
    combined_methods_df %>% pull(m2_pred),
    combined_methods_df %>% pull(m3_pred),
    combined_methods_df %>% pull(m4_pred)
  ))
  
  if(intercept_allowed) x=cbind(1,x)
  
  
  r_inv <- solve(chol(t(x) %*% x))
  c=cbind(rep(1, 4), diag(4))
  
  if(intercept_allowed) c=t(cbind(0, rbind(1, diag(4))))
  
  b <- c(1, rep(0, 4))
  d <- t(combined_methods_df %>% pull(!!as.name(counterfac_var))) %*% x
  nn2 <- sqrt(norm(d, "2"))
  
  
  constr_weights_sol <- solve.QP(Dmat = r_inv * nn2, factorized = TRUE, dvec = d / (nn2^2), Amat = c, bvec = b, meq = 1)
  weight_vec <- constr_weights_sol$solution
  
  if(!constrained) weight_vec <- constr_weights_sol$unconstrained.solution
  
  
  return(weight_vec)
}



# Create a simple, alterative function for using the median as a naive ensemble... mean seems to do fairly well
ensembled_predictor <- function(method1_estimated_df, method2_estimated_df, method3_estimated_df, method4_estimated_df, est_weights,
                                time_var = "period", id_var = "entry", treat_time = "Treatment_Period",
                                pred_var = "point.pred", counterfac_var = "counter_factual",
                                outcome_var = "response") {
  # Estimates the ensemble using weights and predictions from 3 methods
  
  
  # Args
  # method1_estimated_df: long-form dataframe of the estimated and counterfactual point predictions for a given method
  # method2_estimated_df, method3_estimated_df distinct methods of the same form.
  # est_weights: 5 (num methods) by 1 vector of numeric weights to be placed on the predictions of each method
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_time_var: string variable for the column indicating when that unit was treated
  # pred_var: string name of the point estimate column in the methodX dataframe
  # counterfac_var: name of var for the counterfactual value for the time series, if available (otherwise, NULL)
  # outcome_var: the original y variable the methods were estimating
  
  # Output
  # Dataframe with the ensemble predictions by id, period for all time, as well as point effects and counterfactual effects
  # should work even if counter_fac is null
  if (length(est_weights) == 4) {
    ensemble_output <- method1_estimated_df %>%
      dplyr::select(
        !!as.name(time_var), !!as.name(id_var), !!as.name(treat_time),
        !!as.name(outcome_var)
      ) %>%
      mutate(
        point.pred = as.vector(as.matrix(
          cbind(
            method1_estimated_df %>% pull(!!as.name(pred_var)),
            method2_estimated_df %>% pull(!!as.name(pred_var)),
            method3_estimated_df %>% pull(!!as.name(pred_var)),
            method4_estimated_df %>% pull(!!as.name(pred_var))
          )
        ) %*% est_weights)
      ) %>%
      mutate(point.effect = !!as.name(outcome_var) - point.pred)
  }
  # if we have our 4 methods + intercept
  if (length(est_weights) == 5) {
    ensemble_output <- method1_estimated_df %>%
      dplyr::select(
        !!as.name(time_var), !!as.name(id_var), !!as.name(treat_time),
        !!as.name(outcome_var)
      ) %>%
      mutate(
        point.pred = as.vector(cbind(1, as.matrix(
          cbind(
            method1_estimated_df %>% pull(!!as.name(pred_var)),
            method2_estimated_df %>% pull(!!as.name(pred_var)),
            method3_estimated_df %>% pull(!!as.name(pred_var)),
            method4_estimated_df %>% pull(!!as.name(pred_var))
          )
        )) %*% est_weights)
      ) %>%
      mutate(point.effect = !!as.name(outcome_var) - point.pred)
  }
  
  
  if (!is.null(counterfac_var)) {
    ensemble_output <- ensemble_output %>%
      inner_join(
        method1_estimated_df %>%
          dplyr::select(id_var, time_var, counterfac_var),
        by = c(id_var, time_var)
      ) %>%
      mutate(
        cf_point.effect = (!!as.name(outcome_var) - !!as.name(counterfac_var)),
        cf_pct.effect = (!!as.name(outcome_var) / !!as.name(counterfac_var)) - 1
      )
  }
  
  return(ensemble_output)
}







forecomb_helper <- function(indiv_df, time_var = "period", id_var = "entry", treat_time_var = "Treatment_Period",
                            pred_var = "point.pred", counterfac_var = "counter_factual") {
  
  # helper function to estimate the weights for an ensemble using several methods on a placebo set
  # given a single "treated" (placebo) entry, find the weights on the post period estimates from each
  # method that best fit the counterfactual, and return these in a df
  
  
  # Args
  # indiv_df: dataframe formatted in the main function (auto_forecomb_weights_static)
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_time_var: string variable for the column indicating when that unit was treated
  # pred_var: string name of the point estimate column in the methodX dataframe
  # counterfac_var: name of var for the counterfactual value for the time series, if available (otherwise, NULL)
  
  # Output
  # weights, as a 5x1 (num methods+1 X 1) vector, with intercepts and optimal weights
  # if an error is thrown, a simple average with no intercept is used (can change this in the catch loop in forecomb_helper function)
  # Errors arise usually because some of the weights are "moving" (ex, the median of 4 methods will be different over time)
  
  # create the inputs to forecomb pacakge from a DF of individual data
  outcome_vec <- indiv_df %>% pull(!!as.name(counterfac_var))
  treatment_time <- indiv_df %>%
    distinct(!!as.name(treat_time_var)) %>%
    pull()
  
  # create a matrix of the model predictions
  model_pred_mat <- cbind(
    indiv_df %>% pull(m1_pred),
    indiv_df %>% pull(m2_pred),
    indiv_df %>% pull(m3_pred),
    indiv_df %>% pull(m4_pred)
  ) %>% as.matrix()
  
  # use the foreccomb package to structure the data for estimation
  
  # For placebo, where ground truth is known, just find the optimal combination for "post treat" -- where methods were not trained on
  forecomb_data <- foreccomb(observed_vector = outcome_vec[- (1:treatment_time)], prediction_matrix = model_pred_mat[- (1:treatment_time), ])
  
  # attempt to fit the model to find the best weights
  # potential errors here include too few post-treat time to fit, depending on methods.
  # in any case, the simple average is always computable, so if we have an error, we use that
  optim_model <- tryCatch({
    auto_combine(forecomb_data)
  },
  error = function(cond) {
    comb_SA(forecomb_data)
  }
  )
  
  # Several methods either return matrices of weights (for various metrics, weird format however)
  # or char vector informing us there are no static weights.
  # If that occurs, use simple average (no intercept)
  # if the weights exist and seem normal by dimension, use those
  if (is.vector(optim_model$Weights) & length(optim_model$Weights) == 4) {
    # TA, Med, WA, comb_CSR (this gives large numbers...) all do not return weight vectors
    optim_weights <- tibble(
      weight_1 = optim_model$Weights[1], weight_2 = optim_model$Weights[2],
      weight_3 = optim_model$Weights[3], weight_4 = optim_model$Weights[4],
      intercept = ifelse(is.null(optim_model$Intercept), 0, optim_model$Intercept)
    )
  }
  else { # doesn't seem to be a good way to remove specific methods from the grid search
    optim_model <- comb_SA(forecomb_data)
    optim_weights <- tibble(
      weight_1 = optim_model$Weights[1], weight_2 = optim_model$Weights[2],
      weight_3 = optim_model$Weights[3], weight_4 = optim_model$Weights[4],
      intercept = ifelse(is.null(optim_model$Intercept), 0, optim_model$Intercept)
    )
  }
  
  return(optim_weights)
}

# Fix allow for additional methods to foreccomb via ...
# Also, which metric should we be using?
auto_forecomb_weights_static <- function(method1_estimated_df, method2_estimated_df, method3_estimated_df, method4_estimated_df,
                                         time_var = "period", id_var = "entry", treat_time_var = "Treatment_Period",
                                         pred_var = "point.pred", counterfac_var = "counter_factual") {
  
  # Estimates the weights for an ensemble using several methods on a placebo set
  # This method is, as of now, highly specific -- intended for a placebo data set
  # with fake treatment (but zero treatment effect). Given this structure, we have fit several SCM methods
  # on this placebo set to estimate the TE (hopefully close to zero)
  # This method ignores any pre-treat information and outputs weights on the post-treatment point predictions
  # from each method that yields the closest value to the truth. Thus, we must know the truth to compute these.
  
  
  # Args
  # method1_estimated_df: long-form dataframe of the estimated and counterfactual point predictions for a given method
  # method2_estimated_df, method3_estimated_df, method4_estimated_df, distinct methods of the same form.
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_time_var: string variable for the column indicating when that unit was treated
  # pred_var: string name of the point estimate column in the methodX dataframe
  # counterfac_var: name of var for the counterfactual value for the time series, if available (otherwise, NULL)
  
  # Output
  # weights, as a 5x1 (num methods+1 X 1) vector, with intercepts and optimal weights
  # if an error is thrown, a simple average with no intercept is used (can change this in the catch loop in forecomb_helper function)
  
  
  # create a dataframe with the predictions from each of our 4 methods, along with
  # identifying information
  combined_methods_df <- method1_estimated_df %>%
    dplyr::select(
      !!as.name(time_var), !!as.name(id_var), !!as.name(treat_time_var),
      !!as.name(pred_var), !!as.name(counterfac_var)
    ) %>%
    rename(m1_pred = !!as.name(pred_var)) %>%
    left_join(
      method2_estimated_df %>%
        dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(pred_var)) %>%
        rename(m2_pred = !!as.name(pred_var)),
      by = c(id_var, time_var)
    ) %>%
    left_join(
      method3_estimated_df %>%
        dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(pred_var)) %>%
        rename(m3_pred = !!as.name(pred_var)),
      by = c(id_var, time_var)
    ) %>%
    left_join(
      method4_estimated_df %>%
        dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(pred_var)) %>%
        rename(m4_pred = !!as.name(pred_var)),
      by = c(id_var, time_var)
    )
  
  # for each ID (in parallel), send the dataframe over with the counterfactual (what we are trying to fit)
  # as well as the 4 method predictions, which we want to weight to get as close to the counterfactual as possible
  # store resulting weights and intercept in a df (num treated units x num_methods+1)
  weights_by_id_df <- combined_methods_df %>%
    split(.[[id_var]]) %>%
    furrr::future_map(~ forecomb_helper(.)) %>%
    do.call(bind_rows, .) %>%
    mutate(!!as.name(id_var) := method1_estimated_df %>% distinct(!!as.name(id_var)) %>% pull())
  
  return(weights_by_id_df)
}


# Fix allow for additional methods to foreccomb via ...
# Also, which metric should we be using?
forecomb_predictor_static <- function(method1_estimated_df, method2_estimated_df, method3_estimated_df, method4_estimated_df,
                                      placebo_data, est_indiv_weights,
                                      time_var = "period", id_var = "entry", treat_time_var = "Treatment_Period",
                                      pred_var = "point.pred", counterfac_var = "counter_factual",
                                      outcome_var = "response") {
  
  # Estimates the ensemble using weights and predictions from 4 methods
  
  
  # Args
  # method1_estimated_df: AB (true treated) long-form dataframe of the estimated and counterfactual point predictions for a given method
  # method2_estimated_df, method3_estimated_df, method4_estimated_df distinct methods of the same form, AB data!!
  # est_indiv_weights: 5 (num methods+intercept) by 1 vector of numeric weights to be placed on the predictions of each method
  # placebo_data: the full placebo dataset (output of create_placebo_df, or same format), to identify which True treated
  # observations the weights belong to.
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_time_var: string variable for the column indicating when that unit was treated
  # pred_var: string name of the point estimate column in the methodX dataframe
  # counterfac_var: name of var for the counterfactual value for the time series, if available (otherwise, NULL)
  # outcome_var: the original y variable the methods were estimating
  
  # Output
  # Dataframe with the ensemble predictions by id, period for all time, as well as point effects and counterfactual effects
  # should work even if counter_fac is null
  
  # Each "treated" entry in our placebo set now has an estimated weight.
  # We want to find the true treated unit that corresponds to this placebo treated unit, so
  # we can assign them the proper weight
  est_indiv_weights_match <- est_indiv_weights %>%
    left_join(
      placebo_data %>% dplyr::select(!!as.name(id_var), Treatment_Unit) %>% distinct(entry, .keep_all = T),
      by = id_var
    ) %>%
    select(-!!as.name(id_var)) %>%
    rename(!!as.name(id_var) := Treatment_Unit)
  
  # For each true treated prediction, combine according to the unit specific weights
  # create a dataframe that has relevant period/entry long form, append the weights by Id
  # then apply the weights over each method, and compute treatment effects
  ensemble_output <- method1_estimated_df %>%
    dplyr::select(
      !!as.name(time_var), !!as.name(id_var), !!as.name(treat_time_var),
      !!as.name(outcome_var)
    ) %>%
    left_join(
      est_indiv_weights_match,
      by = id_var
    ) %>%
    mutate(
      point.pred = method1_estimated_df %>% pull(!!as.name(pred_var)) * weight_1 +
        method2_estimated_df %>% pull(!!as.name(pred_var)) * weight_2 +
        method3_estimated_df %>% pull(!!as.name(pred_var)) * weight_3 +
        method4_estimated_df %>% pull(!!as.name(pred_var)) * weight_4 +
        intercept
    ) %>%
    mutate(
      point.effect = !!as.name(outcome_var) - point.pred,
      pct.effect = (!!as.name(outcome_var) / point.pred) - 1
    )
  
  # if coutnerfactual AB data exists, we can compute that here and add it to our output
  if (!is.null(counterfac_var)) {
    ensemble_output <- ensemble_output %>%
      inner_join(
        method1_estimated_df %>%
          dplyr::select(!!as.name(id_var), !!as.name(time_var), !!as.name(counterfac_var)),
        by = c(id_var, time_var)
      ) %>%
      mutate(
        cf_point.effect = (!!as.name(outcome_var) - !!as.name(counterfac_var)),
        cf_pct.effect = (!!as.name(outcome_var) / !!as.name(counterfac_var)) - 1
      )
  }
  
  return(ensemble_output)
}