pacman::p_load(dplyr, ggplot2, quadprog, purrr, furrr, tidyr, glue, tibble)

#' Reduce a list of predictions (different methods, same data) into a tibble.
#'
#' @param pred_list List of long form tibbles of estimated series, one per
#'    method to ensemble.
#' @inheritParams EstimateEnsemble
#' @param treat_time_var String name for the treated time from the method call.
#'    Fixed within the method call, so default option should not require change.
#' @param pred_var String name for the imputed outcome from the method call.
#'    Fixed within the method call, so default option should not require change.
#' @param response_var String name for the true outcome from the method call.
#'    Fixed within the method call, so default option should not require change.
#'
#' @return A single tibble of the prediction series for each method, long form.
.ReducePredList <- function(pred_list, time_var, id_var, counterfac_var,
                            treat_time_var = "Treatment_Period",
                            pred_var = "point.pred",
                            response_var = "response") {
  # Combine list of estimates into a single tibble by joining on true variables.
  # Rename duplicated point.pred columns to m_pred`#`, for each method number.
  tib_out <- pred_list %>%
    purrr::reduce(dplyr::left_join,
      by = c(id_var, time_var, treat_time_var, response_var)
    ) %>%
    dplyr::rename(
      !!purrr::set_names(
        tidyselect::contains(pred_var),
        c(glue::glue("m_pred{i}",
          i = seq_along(pred_list)
        ))
      )
    ) %>%
    dplyr::select(
      tidyselect::all_of(c(
        time_var, id_var, treat_time_var, response_var,
        c(glue::glue("m_pred{i}",
          i = seq_along(pred_list)
        ))
      ))
    )
  # If counterfactual is specified, add relevant columns to our combined tibble.
  if (!is.null(counterfac_var)) {
    tib_out <- tib_out %>%
      dplyr::inner_join(pred_list[[1]] %>%
        dplyr::select(tidyselect::all_of(c(
          counterfac_var,
          time_var, id_var,
          "cf_point.effect",
          "cf_pct.effect"
        ))),
      by = c(time_var, id_var)
      )
  }
  return(tib_out)
}

#' Estimates optimal weights of the given methods, constraints, and data.
#'
#' @param to_weight_tib Input data, either an individual series tibble or a full
#'    dataset with all units. Contains the various method predictions and true
#'    outcomes used to estimate.
#' @inheritParams .EnsemblePlaceboWeights
#'
#' @return Vector of the optimal combination of methods to predict the true
#'    outcome given the constraints.
.WeightSolver <- function(to_weight_tib, intercept_allowed, constrained,
                          method_names,
                          counterfac_var) {
  stopifnot(method_names %in%
              c("CausalImpact", "Gsynth", "SDID", "SDID_Uncon", "SCM"))
  stopifnot(length(method_names)>0)
  # Create a matrix of the prediction series by method.
  x <- to_weight_tib %>%
    dplyr::select(tidyselect::contains("m_pred")) %>%
    as.matrix()
  # Add intercept.
  if (intercept_allowed) x <- cbind(1, x)
  # Attempt to compute the Cholesky decomposition, throw "Failure" if error.
  r_inv <- tryCatch(
    {
      solve(chol(t(x) %*% x))
    },
    error = function(e) {
      "Failure"
    }
  )
  # Create the constraints for the optimization without an intercept.
  # Namely, positive weight on each method, and that they sum to 1.
  c <- cbind(rep(1, length(method_names)), diag(length(method_names)))
  # Amend the constraints if an intercept is allowed (can be negative).
  if (intercept_allowed) c <- t(cbind(0, rbind(1, diag(length(method_names)))))
  # Set up the remaining inputs for the quadratic programming solver.
  b <- c(1, rep(0, length(method_names)))
  d <- t(to_weight_tib %>% pull(!!as.name(counterfac_var))) %*% x
  nn2 <- sqrt(norm(d, "2"))
  # Attempt to solve for the optimal weights. If there is an error, return
  # the average of each of the methods and an intercept of 0.
  weight_vec <- tryCatch(
    {
      constr_weights_sol <- solve.QP(
        Dmat = r_inv * nn2, factorized = TRUE,
        dvec = d / (nn2^2), Amat = c, bvec = b,
        meq = 1
      )
      if (!constrained) {
        constr_weights_sol$unconstrained.solution
      } else {
        constr_weights_sol$solution
      }
    },
    error = function(e) {
      if (intercept_allowed) {
        error_out <- c(0, rep(1 / length(method_names), length(method_names)))
      } else {
        rep(1 / length(method_names), length(method_names))
      }
    }
  )
  return(weight_vec)
}

#' Computes the optimal weights for imputing placebo outcomes in post treatment.
#'
#' @inheritParams EstimateEnsemble
#' @param combined_methods_tib Long form tibble containing the predicted value
#'    and true value of the outcome of interest for each placebo treated unit,
#'    and each method input.
#' @param treat_time_var String name for the treated time from the method call.
#'    Fixed within the method call, so default option should not require change.
#'
#' @return Tibble with the weights associated with each method. Weights are
#'    determined as a function of the user inputted constraints, and estimated
#'    as the best combination of the methods for predicting all post treatment
#'    period outcomes in the placebo data.
.EnsemblePlaceboWeights <- function(method_names, 
                                    combined_methods_tib,
                                    indiv_weights, constrained,
                                    intercept_allowed, time_var,
                                    id_var, counterfac_var,
                                    treat_time_var = "Treatment_Period") {
  stopifnot(method_names %in%
              c("CausalImpact", "Gsynth", "SDID", "SDID_Uncon", "SCM"))
  # Store only the post treatment period of our predictions.
  post_treat_combined_tib <- combined_methods_tib %>%
    dplyr::filter(!!as.name(time_var) >= !!as.name(treat_time_var))

  if (indiv_weights) {
    # Estimate unique weights on each method per individual.
    # Split tibble into lists of series per ID, and estimate separate weights.
    list_to_ensemble <- post_treat_combined_tib %>% split(.[[id_var]])
    tib_result <- furrr::future_map(
      .x = list_to_ensemble,
      .f = ~ .WeightSolver(
        to_weight_tib = .x,
        intercept_allowed = intercept_allowed,
        constrained = constrained,
        method_names = method_names,
        counterfac_var = counterfac_var
      )
    ) %>%
      dplyr::bind_rows() %>%
      t() %>%
      tibble::as_tibble(rownames = id_var, .name_repair = "unique") %>%
      setNames(c("entry", paste0("V", seq_len(ncol(.)))))
  } else {
    # Estimate one set of weights common to all individuals.
    df_result <- .WeightSolver(
      to_weight_tib = post_treat_combined_tib,
      intercept_allowed = intercept_allowed,
      constrained = constrained,
      method_names = method_names,
      counterfac_var = counterfac_var
    ) %>%
      as.list() %>%
      data.frame()
    # Combine the data frame into a tibble and clean the names a bit.
    tib_result <- df_result %>%
      tibble::tibble(.name_repair = "universal") %>%
      setNames(c(paste0("V", seq_len(ncol(.)))))
  }
  # Tidy up the format of the weights tibble to return.
  if (intercept_allowed) {
    tib_result <- tib_result %>% dplyr::rename(Intercept = "V1")
  }
  tib_result <- tib_result %>%
    dplyr::rename(
      !!purrr::set_names(
        tidyselect::contains("V"),
        c(glue::glue("m_weight{i}",
          i = seq_along(method_names)
        ))
      )
    )
  return(tib_result)
}

#' Applies the estimator associated with function name to input data.
#'
#' @param method_name String name for the desired method to estimate.
#' @param data_name Dataset to apply the estimates to given the method.
#' @param ... Additional parameters to pass on to estimator, such as column
#'    column names for the given data set (id_var, time_var).
#'
#' @return A long form tibble of the estimated series with imputed potential
#'    untreated outcome for each period and treated unit.
.MethodEstimator <- function(method_name, 
                              data_name, ...) {
  stopifnot(method_name %in%
              c("CausalImpact", "Gsynth", "SDID", "SDID_Uncon", "SCM"))
  # Define the estimator call, and apply.
  # MC and SDID Unconstrained are exceptions to the rule (no unique function).
  if (method_name == "MC") {
    return(EstimateGsynthSeries(data_name, estimator = "mc", ...))
  }
  if (method_name == "SDID_Uncon") {
    return(EstimateSDIDSeries(data_name, constrained = F, ...))
  }
  # Define the call to the function based on method name, and cal.
  estimator_call <- paste0("Estimate", method_name, "Series")
  est_out <- do.call(estimator_call, list(data_name, ...))
  return(est_out)
}

#' Estimates ensemble outcomes combining methods with placebo optimal weights.
#'
#' @param method_names Vector of strings for the names of methods to ensemble.
#'    Names must be one of the potential methods, will be used to call
#'    Estimate"method_names"Series.
#' @param true_data Tibble of the formatted raw data to estimate ensemble for.
#'    This should be a true panel in that every unit has an observed outcome
#'    for every time period in the data set -- remove NAs beforehand.
#' @param pred_list List of tibbles containing the estimated series per method.
#' @param constrained Boolean for whether to constrain the ensemble to have
#'    non-negative weights that sum to one over the input methods.
#' @param intercept_allowed Boolean for whether to include an intercept
#'    adjustment in the ensemble.
#' @param indiv_weights Boolean for whether to estimate unit-specific weights.
#' @param time_var String name of the time variable in the raw data.
#' @param id_var String name of the unit variable in the raw data.
#' @param outcome_var String name of the outcome variable in the raw data.
#' @param treat_indicator String name of treatment indicator in the raw data.
#' @param counterfac_var String name of the counterfactual variable in the raw
#'    data. Can be NULL, as the placebo creates it's own counterfactual.
#'
#' @return A long form tibble with the predictions from each method, their
#'    associated weights, and the final ensembled point estimate as well as the
#'    effects.
# TODO(alexdkellogg): Potential improvement - focus ensemble weights on k
#    periods after treatment, rather than the full post-treat time period.
# TODO(alexdkellogg): Update the ensemble so that it averages over several
#    placebo datasets. Learn time vs accuracy trade-off.
EstimateEnsemble <- function(method_names, 
                             true_data, pred_list,
                             constrained = T, intercept_allowed = T,
                             indiv_weights = F, time_var = "period",
                             id_var = "entry", outcome_var = "target",
                             treat_indicator = "treatperiod_0",
                             counterfac_var = "counter_factual") {
  # Ensure the method names are appropriately specified.
  stopifnot(method_names %in%
              c("CausalImpact", "Gsynth", "SDID", "SDID_Uncon", "SCM"))
  # Check that dimension of the weight vector equals the number of methods.
  stopifnot(length(method_names) == length(pred_list))
  # Create a placebo dataset from true_data via matching.
  placebo_data <- CreatePlaceboData(true_data,
    id_var = id_var,
    time_var = time_var,
    outcome_var = outcome_var,
    counterfac_var = counterfac_var
  )
  # Store the mapping between true treated unit and its placebo treated unit.
  unit_mapping <- placebo_data %>%
    dplyr::filter(!is.na(Treatment_Period)) %>%
    dplyr::distinct(entry, Treatment_Unit)
  # Estimate the counterfactual series for each method.
  # Bug in future library requires calling function names to have them
  # recognized as functions in each parallel instance.
  EstimateBARTSeries
  EstimateGfooSeries
  EstimateGsynthSeries
  EstimateSDIDSeries
  EstimateCausalImpactSeries
  EstimateSCMSeries
  # Map the desired (default) estimators onto the placebo set.
  estimates_list <- furrr::future_map(
    .x = method_names, .f = ~.MethodEstimator( method_name = .x,
    data_name = placebo_data, id_var = id_var, time_var = time_var,
    outcome_var = outcome_var, treat_indicator = treat_indicator,
    counterfac_var = counterfac_var))
  names(estimates_list) <- method_names
  # Create the default counterfactual variable, which exists in placebo data.
  temp_cf <- ifelse(is.null(counterfac_var), "counter_factual", counterfac_var)
  # Combine the list of estimates (one per method) into a single tibble.
  combined_pred <- .ReducePredList(estimates_list,
    id_var = id_var,
    time_var = time_var,
    counterfac_var = temp_cf
  )
  # Estimate the weights on each of the methods
  method_weights <-
    .EnsemblePlaceboWeights(method_names,
      combined_methods_tib = combined_pred,
      indiv_weights = indiv_weights,
      constrained = constrained,
      intercept_allowed = intercept_allowed,
      time_var = time_var,
      id_var = id_var,
      counterfac_var = temp_cf
    )

  # Combine the existing estimates on true data into tibble to extract the
  # predictions from each method.
  combined_true_methods <- .ReducePredList(pred_list,
    id_var = id_var,
    time_var = time_var,
    counterfac_var = counterfac_var
  )
  if (indiv_weights) {
    # If unit specific weights were estimated, map the weights on placebo data
    # to the true data, and estimate the ensemble.
    matched_weights <- unit_mapping %>%
      dplyr::inner_join(method_weights %>% dplyr::mutate_all((as.numeric)),
        by = id_var
      ) %>%
      dplyr::select(-tidyselect::all_of(id_var)) %>%
      dplyr::rename(entry = "Treatment_Unit")

    combined_true_methods <- combined_true_methods %>%
      dplyr::left_join(matched_weights, by = id_var)
  } else {
    combined_true_methods <- combined_true_methods %>%
      dplyr::bind_cols(method_weights)
  }
  weighted_out <- .WeightPredictions(combined_true_methods, intercept_allowed)
  return(weighted_out)
}

#' Combines the estimated weights per method with method estimates.
#'
#' @param data Long form tibble with each row containing the estimates from
#'    each method and their optimal weight on the placebo unit.
#' @inheritParams EstimateEnsemble
#'
#' @return Tibble with the final, ensembled prediction per unit and time.
.WeightPredictions <- function(data, intercept_allowed) {
  # Extract the predictions from each method.
  method_preds <- data %>%
    dplyr::select(tidyselect::contains("m_pred")) %>%
    as.matrix()
  # Add intercept if needed.
  if (intercept_allowed) method_preds <- cbind(1, method_preds)
  # Extract the weights for each method.
  method_weights <- data %>%
    dplyr::select(tidyselect::matches("Intercept|m_weight")) %>%
    as.matrix()
  # Compute the linear combination of predictions and weights, compute effects.
  data <- data %>%
    dplyr::mutate(
      "point.pred" = c(rowSums(method_preds * method_weights)),
      "point.effect" = response - point.pred
    )
  return(data)
}
