
# TODO(alexdkellogg): Ensure compatibility with "MC", "SDID Unconstrained",
#    as well as any new methods (in the stopifnot and .EstimateMethod calls).
#    Still getting errors with "...", as well.

#' Compute the ATT of a given method after removing placebo estimated bias.
#' ATT = Average Treatment Effect on the Treated.
#'
#' @param raw_data Tibble of the formatted raw data to estimate ensemble for.
#'    This should be a true panel in that every unit has an observed outcome
#'    for every time period in the data set -- remove NAs beforehand.
#' @param method_name String name for the desired method to estimate.
#' @param num_placebos Number of placebo datasets used to estimate the placebo
#'    bias. The placebo is created once by matching, and the distinct datasets
#'    are drawn by adding random iid noise.
#' @param sd_placebo Standard deviation of the iid noise for the placebo data.
#' @param alpha_ci Confidence level for the ATT confidence bands.
#' @param post_treat_only Boolean variable for whether to only return post-
#'    treatment ATT estimates.
#' @param effect_type String name for the type of treatment effect to compute.
#'    Currently support one of "mean", "median", and "aggregate".
#' @param id_var String name for the unit ID variable.
#' @param time_var String name of the time column.
#' @param treat_indicator String name for the treatment time.
#' @param outcome_var String name of the true outcome variable.
#' @param counterfac_var String name for the true outcome.
#'
#' @return Tibble with the jackknife estimates of the desired debiased effect
#'    type by time to treatment, alongside the number of treated units in that
#'    time period and the confidence bounds (for both percent and absolute TE).
# TODO(alexdkellogg): Check whether I should rename id_var to be "entry" in here
#    so I can just call entry in the main, rather than as.name(id_var).
# TODO(alexdkellogg): May not need to call ComputeTE, don't need CI.
DebiasedATT <- function(raw_data, method_name, num_placebos = 50,
                        sd_placebo = 0.2, alpha_ci = 0.95, post_treat_only = F,
                        effect_type = "mean",
                        id_var = "entry", time_var = "period",
                        treat_indicator = "treatperiod_0",
                        outcome_var = "target",
                        counterfac_var = "counter_factual") {
  # Ensure the method names are appropriately specified.
  stopifnot(method_name %in%
    c("CausalImpact", "Gsynth", "SDID", "SDID_Uncon", "SCM"))
  # Use rules of thumb to try to create an independent placebo set.
  indep_data <- .CreateIndepPlacebo(raw_data,
    id_var = id_var,
    time_var = time_var,
    treat_indicator = treat_indicator,
    outcome_var = outcome_var,
    counterfac_var = counterfac_var
  )

  placebo_data <- indep_data[["placebo"]]
  true_subset <- indep_data[["true"]]

  # Bug in future library requires calling function names to have them
  # recognized as functions in each parallel instance.
  EstimateGsynthSeries
  EstimateSDIDSeries
  EstimateCausalImpactSeries
  EstimateSCMSeries

  # Estimate the imputed prediction series using helper function.
  pred_series_tib <- .MethodEstimator(
    method_name = method_name,
    data_name = true_subset,
    id_var = id_var, time_var = time_var,
    outcome_var = outcome_var,
    treat_indicator = treat_indicator,
    counterfac_var = counterfac_var
  )
  # Compute the Treatment Effect series given the estimated series.
  est_att_series <- ComputeTreatmentEffect(pred_series_tib,
    effect_type = effect_type,
    post_treat_only = post_treat_only,
    pred_var = "point.pred",
    alpha_ci = alpha_ci
  ) %>%
    dplyr::select(c(
      "post_period_t", "treated_n",
      tidyselect::matches("observed|cf")
    ))
  # Compute the ATT of the distinct placebo sets.
  placebo_att <-
    .EstimatePlaceboList(
      placebo_data = placebo_data, method_name = method_name,
      num_placebos = num_placebos, outcome_var = outcome_var,
      id_var = id_var, time_var = time_var,
      treat_indicator = treat_indicator,
      counterfac_var = counterfac_var,
      effect_type = effect_type, sd_placebo = sd_placebo,
      post_treat_only = post_treat_only
    )

  # Combine the placebo and true ATT and debias in this helper function.
  jk_comp_ci <- .JackknifeDebiasedATT(
    placebo_att = placebo_att,
    est_att_series = est_att_series,
    effect_type = effect_type,
    alpha_ci = alpha_ci
  )
  return(jk_comp_ci)
}

#' Heuristic creation of an independent placebo data set when data allows.
#'
#' @inheritParams DebiasedATT
#'
#' @return A list of two tibbles: the first (named placebo) stores the placebo
#'    data tibble in long form and the second (named true) stores the data to
#'    estimate. If there is sufficient data, these two tibbles will have no
#'    overlap.
.CreateIndepPlacebo <- function(raw_data,
                                id_var = "entry",
                                time_var = "period",
                                treat_indicator = "treatperiod_0",
                                outcome_var = "target",
                                counterfac_var = "counter_factual") {
  # Identify the treated entries.
  treated_ids <- .IdentifyTreated(raw_data,
    id_var = id_var,
    treat_indicator = treat_indicator
  )
  num_treated <- treated_ids %>%
    dplyr::distinct(!!as.name(id_var), .keep_all = T) %>%
    dplyr::filter(Treated == 1) %>%
    dplyr::count() %>%
    dplyr::pull()
  num_donors <- treated_ids %>%
    dplyr::distinct(!!as.name(id_var), .keep_all = T) %>%
    dplyr::filter(Treated == 0) %>%
    dplyr::count() %>%
    dplyr::pull()
  # Split the data cleanly if we have adequate data.
  true_full_ind <- NULL
  if (num_donors / num_treated >= 5 ) {
    num_slice <- 3 * num_treated
    frac_treated <- 1
  } else if (num_donors / num_treated >= 3) {
    num_slice <- num_treated
    frac_treated <- 0.5
  } else {
    # If we do not have adequate data, create a placebo with the full raw data.
    # Output warning to user as well.
    num_slice <- num_donors
    frac_treated <- 1
    # Set flag that true_subset should be the full data.
    true_full_ind <- 1
    warning("Too few donor units relative to treated units to create an
            independent placebo set. Debiasing may suffer as a result.")
  }
  for_placebo <- .SelectDonorIDs(raw_data = raw_data, 
                                 treated_ids = treated_ids,
                                 id_var = id_var, 
                                 num_slice = num_slice,
                                 frac_treated = frac_treated)
  placebo_temp <- CreatePlaceboData(for_placebo)
  if (!is.null(true_full_ind)) {
    true_subset <- raw_data
  } else{
    true_subset <- raw_data %>%
      dplyr::anti_join(placebo_temp, by=c(id_var, time_var))
  }
  list_out <- list(placebo_temp, true_subset)
  names(list_out) <- c("placebo", "true")
  return(list_out)
}

#' Selects a random set of IDs given the tibbles and sampling numbers.
#'
#' @inheritParams DebiasedATT
#' @param treated_ids A long form tibble with a column identifying treatment
#'    status, as well as the outcome series per unit and time combination.
#' @param num_slice The number of donor units to randomly sample for the placebo
#'    holdout set.
#' @param frac_treated The fraction of treated units to be matched on for the
#'    placebo treated set.
#'
#' @return A subset of the full data tibbles sampled randomly, with all relevant
#'    information for each selected observation over time.
.SelectDonorIDs <- function(raw_data, treated_ids, id_var, num_slice,
                            frac_treated) {
  # Slice sample of donor IDs, defaults to without replacement.
  set_aside_ids <- treated_ids %>%
    dplyr::filter(Treated == 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::slice_sample(n = num_slice) %>%
    dplyr::bind_rows(treated_ids %>%
    dplyr::filter(Treated == 1) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::slice_sample(prop = frac_treated))

  for_placebo <- set_aside_ids %>%
    dplyr::inner_join(raw_data)
  return(for_placebo)
}

#' Computes the debiased ATT given the list of placebo data, and formats CI.
#'
#' @param placebo_att List of the placebo estimated ATT (which is equal to the
#'    bias, since the true ATT should always be 0) for each draw.
#' @param est_att_series Tibble of the ATT for the true data.
#' @inheritParams DebiasedATT
#'
#' @return A formatted tibble containing the estimated (placebo debiased) ATT
#'   and the jackknifed CI bounds on the ATT of the given effect.
.JackknifeDebiasedATT <- function(placebo_att, est_att_series, effect_type,
                                  alpha_ci) {
  abs_predvar <- paste("observed", effect_type, "abs", "att", sep = "_")
  pct_predvar <- paste("observed", effect_type, "pct", "att", sep = "_")
  # Each ATT point estimate is a an estimate of the bias for the placebo.
  # Jackknife a debiased ATT using each of these estimates and the true data.
  combined_att_tib <- lapply(placebo_att, function(x) {
    x %>%
      dplyr::rename(
        placebo_abs_bias = abs_predvar,
        placebo_pct_bias = pct_predvar
      ) %>%
      dplyr::select(c(post_period_t, placebo_abs_bias, placebo_pct_bias)) %>%
      dplyr::inner_join(est_att_series, by = "post_period_t") %>%
      dplyr::mutate(
        new_abs_att = !!as.name(abs_predvar) - placebo_abs_bias,
        new_pct_att = !!as.name(pct_predvar) - placebo_pct_bias
      )
  }) %>% dplyr::bind_rows()
  # Split the tibble by post treatment period.
  effects_to_jk <- combined_att_tib %>%
    split(.[["post_period_t"]])
  # Compute jackknife estimates of the mean of the debiased ATT.
  jk_comp_abs <- effects_to_jk %>%
    furrr::future_map(~ .[["new_abs_att"]]) %>%
    furrr::future_map(~ resample::jackknife(., mean))
  jk_comp_pct <- effects_to_jk %>%
    furrr::future_map(~ .[["new_pct_att"]]) %>%
    furrr::future_map(~ resample::jackknife(., mean))

  # Compute the relevant upper and lower CI bounds for the debiased ATT.
  # Gather the bounds on the confidence intervals and format as tibble.
  # Join in counterfactual information from the original series, if it exits.
  post_range <- c(
    min(combined_att_tib$post_period_t),
    max(combined_att_tib$post_period_t)
  )
  jk_comp_ci <- .FormatJackknife(
    jk_est = jk_comp_abs, ci_in = alpha_ci, col_name =
      paste0(effect_type, "_abs_att"), range = post_range
  ) %>%
    dplyr::inner_join(
      .FormatJackknife(
        jk_est = jk_comp_pct, ci_in = alpha_ci, col_name =
          paste0(effect_type, "_pct_att"), range = post_range
      ),
      by = "post_period_t"
    ) %>%
    dplyr::inner_join(
      est_att_series %>%
        dplyr::select(c(
          "post_period_t", "treated_n",
          tidyselect::matches("cf")
        )),
      by = "post_period_t"
    )

  # Reorder the columns, and return.
  jk_comp_ci <- jk_comp_ci %>%
    dplyr::select(c(post_period_t, treated_n), dplyr::everything())
  return(jk_comp_ci)
}

#' Creates and estimates the ATT series over a list of placebo draws.
#'
#' @param placebo_data Tibble (long form) of the placebo treated and donor units
#'    with target and counterfactual outcome (equal) for each unit and time.
#' @inheritParams DebiasedATT
#'
#' @return A list of tibbles, each containing the estimated ATT (here, also
#'    equal to bias) for each placebo draw.
.EstimatePlaceboList <- function(placebo_data, method_name, num_placebos = 50,
                                 sd_placebo = 0.2, outcome_var = "target",
                                 id_var = "entry", time_var = "period",
                                 treat_indicator = "treatperiod_0",
                                 counterfac_var = "counter_factual",
                                 effect_type = "mean",
                                 post_treat_only = F) {
  # Create a list of placebos and add random noise.
  # Note that matching from treated unit to placebo treated unit is done just
  # once, so the 50 draws are for the same set of placebo treated units.
  list_placebos <- lapply(seq_len(num_placebos), function(x) placebo_data)
  noised_placebos <- furrr::future_map(
    .x = list_placebos,
    .f = function(z) {
      z %>%
        dplyr::mutate(
          !!as.name(outcome_var) :=
            !!as.name(outcome_var) +
              stats::rnorm(
                dplyr::n(), 0,
                sd_placebo
              ),
          "counter_factual" =
            !!as.name(outcome_var)
        )
    }
  )
  # Estimate the outcome series for each placebo in the list.
  # Map the desired (default) estimators onto the placebo set.
  placebo_estimates <- furrr::future_map(
    .x = noised_placebos, .f = ~ .MethodEstimator(
      data_name = .x,
      method_name = method_name, id_var = id_var, time_var = time_var,
      outcome_var = outcome_var, treat_indicator = treat_indicator,
      counterfac_var = counterfac_var
    )
  )
  # Estimate the ATT over the placebo sets.
  placebo_att <- furrr::future_map(
    .x = placebo_estimates,
    .f = ~ ComputeTreatmentEffect(
      pred_series_tib = .x,
      effect_type = effect_type,
      post_treat_only = post_treat_only,
      pred_var = "point.pred"
    )
  )
  return(placebo_att)
}

#' Compares desire method performance (RMSE) on placebo data and computes
#'    estimates for the winner.
#'
#' @param method_names Vector of string names for desired methods to estimate.
#' @param horizon Range of periods over which to display the output. Should be
#'    non-negative here because several methods perfectly fit pre-treatment
#'    periods, which would result in more desirable RMSE properties.
#' @inheritParams DebiasedATT
#'
#' @return A long form tibble with the winning method (given input) estimates of
#'    the outcome variable (missing outcome imputation) as well as the point
#'    and percent effects, and a message informing the user of the winner.
# TODO(alexdkellogg): Could create a similar function that instead returns
# a table with Bias, RMSE, coverage of each method on placebo draws?
# TODO(alexdkellogg): Append with "break" functionality if the data fails 
# diagnostic tests (t-test table, or KL metric on feature density).
AutoEstimator <- function(raw_data, method_names,
                          num_placebos = 50, sd_placebo = 0.2,
                          horizon = c(0, 4), effect_type = "mean",
                          id_var = "entry",
                          time_var = "period",
                          treat_indicator = "treatperiod_0",
                          outcome_var = "target",
                          counterfac_var = "counter_factual") {
  # Ensure the method names are entered appropriately.
  stopifnot(method_names %in%
    c("CausalImpact", "Gsynth", "SDID", "SDID_Uncon", "SCM"))
  # Make sure horizon is positive. Some methods have 0 training error, so 
  # including negative horizon is neither useful nor comparable across methods.
  stopifnot(0 <= horizon[1] & horizon[1] <= horizon[2])
  # Bug in future library requires calling function names to have them
  # recognized as functions in each parallel instance.
  EstimateGsynthSeries
  EstimateSDIDSeries
  EstimateCausalImpactSeries
  EstimateSCMSeries

  # Use rules of thumb to try to create an independent placebo set.
  indep_data <- .CreateIndepPlacebo(raw_data,
    id_var = id_var,
    time_var = time_var,
    treat_indicator = treat_indicator,
    outcome_var = outcome_var,
    counterfac_var = counterfac_var
  )
  placebo_data <- indep_data[["placebo"]]
  true_subset <- indep_data[["true"]]

  # Get the ATT (or bias, since this is placebo) for each method over
  # several seed (num_placebos).
  placebo_att <-
    purrr::map(.x = c(method_names), .f = ~ .EstimatePlaceboList(
      method_name = .x,
      placebo_data = placebo_data, num_placebos = num_placebos,
      sd_placebo = sd_placebo,
      outcome_var = outcome_var, id_var = id_var, time_var = time_var,
      treat_indicator = treat_indicator,
      counterfac_var = counterfac_var, effect_type = effect_type,
      post_treat_only = F
    ))
  # Combine the average bias and metrics over the seeds into a tibble.
  bias_tib <- lapply(placebo_att[seq_along(method_names)], ComputeBias,
    horizon = horizon
  ) %>%
    purrr::reduce(dplyr::left_join, by = "post_period_t")
  metrics_tib <- lapply(placebo_att[seq_along(method_names)], ComputeMetrics,
    horizon = horizon
  ) %>%
    purrr::reduce(dplyr::left_join, by = "post_period_t") %>%
    dplyr::select(c("post_period_t", tidyselect::contains("rmse")))
  # TODO(alexdkellogg): Is there a better composite score to use? Issue with
  # bias is it can be negative and positive (cancels), so take abs or square.
  best_method <-
    colMeans(metrics_tib %>% dplyr::select(-post_period_t)) %>% which.min()
  # Define the winning method and call it.
  auto_method <- method_names[best_method]

  # Map the desired estimator onto the placebo set.
  pred_series_auto <- .MethodEstimator(
    method_name = auto_method,
    data_name = raw_data,
    outcome_var = outcome_var,
    id_var = id_var, time_var = time_var,
    treat_indicator = treat_indicator,
    counterfac_var = counterfac_var
  )
  # TODO(alexdkellogg): Is there a better way to pass this information to user?
  message(paste(
    auto_method, "was RMSE minimizing estimator on", num_placebos,
    "draws of the placebo given a horizon from", horizon[1],
    "to", horizon[2], "relative to treatment."
  ))
  return(pred_series_auto)
}
