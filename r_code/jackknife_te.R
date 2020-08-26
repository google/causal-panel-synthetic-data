pacman::p_load(dplyr, furrr, stringr, tidyr, tibble, resample)

#' Compute the true treatment effects for the data by post treatment time.
#'
#' @inheritParams ComputeTreatmentEffect
#'
#' @return Tibble with the post treatment time period and the true treatment
#'    effect (percent and absolute) of the type given by effect_type.
.ComputeCounterfactualTE <- function(pred_series_tib,
                                     effect_type = c(
                                       "mean", "median",
                                       "aggregate"
                                     ), time_var,
                                     treat_period_var,
                                     outcome_var,
                                     counterfac_var) {
  effect_type <- match.arg(effect_type)
  # Column names for output.
  cf_abs_col <- paste0(effect_type, "_abs_cf_att")
  cf_pct_col <- paste0(effect_type, "_pct_cf_att")

  # Computing the relevant effect and storing as tibble.
  cf_att_df <- pred_series_tib %>%
    dplyr::mutate(post_period_t = !!as.name(time_var) -
      !!as.name(treat_period_var)) %>%
    dplyr::group_by(post_period_t) %>%
    dplyr::summarise(
      !!as.name(cf_abs_col) := dplyr::case_when(
        effect_type == "aggregate" ~ sum(!!as.name(outcome_var) -
          !!as.name(counterfac_var)),
        effect_type == "mean" ~ mean(!!as.name(outcome_var) -
          !!as.name(counterfac_var)),
        effect_type == "median" ~ median(!!as.name(outcome_var) -
          !!as.name(counterfac_var))
      ),
      !!as.name(cf_pct_col) := dplyr::case_when(
        effect_type == "aggregate" ~ 
          ifelse(sum(!!as.name(counterfac_var))==0, NA_real_,
                 (sum(!!as.name(outcome_var)) /
                    sum(!!as.name(counterfac_var))) - 1  
                 ),
        effect_type == "mean" ~
          ifelse(!!as.name(counterfac_var)==0, NA_real_,
                 mean((!!as.name(outcome_var) /
                         !!as.name(counterfac_var)) - 1)
                 ),
        effect_type == "median" ~ 
          ifelse(!!as.name(counterfac_var)==0, NA_real_,
                 median((!!as.name(outcome_var) /
                           !!as.name(counterfac_var)) - 1)               
                 )
      )
    ) %>%
    dplyr::ungroup()

  return(cf_att_df)
}

#' Compute jackknife estimates of treatment effects and confidence intervals.
#'
#' @param pred_series_tib Tibble of predicted and true series, long form. That
#'    is, each row should represent a unit in a particular time, and have the
#'    data on their true outcome as well as a prediction from some method.
#' @param time_var String name of the time column.
#' @param treat_period_var String name for the treatment time.
#' @param pred_var String name of the predicted outcome from some method.
#' @param outcome_var String name of the true outcome variable.
#' @param effect_type String name for the type of treatment effect to compute.
#'    Currently support one of "mean", "median", and "aggregate".
#' @param alpha_ci Confidence level desired for the jackknife bounds.
#' @param counterfac_var String name for the true outcome.
#' @param post_treat_only Boolean indicator for whether to ignore pre-treatment.
#'
#' @return Tibble with the jackknife estimates of the desired effect type by
#'    time to treatment, alongside the number of treated units in that time
#'    period and the confidence bounds (for both percent and absolute TE).
# TODO(alexdkellogg): I think the jackknifing of the "aggregate" TE is not doing
#    what we want. Fix this -- for now, restricting functionality.
ComputeTreatmentEffect <- function(pred_series_tib,
                                   time_var = "period",
                                   treat_period_var = "Treatment_Period",
                                   pred_var = "point.pred",
                                   outcome_var = "response",
                                   effect_type = c("mean", "median"),
                                   alpha_ci = 0.95,
                                   counterfac_var = "counter_factual",
                                   post_treat_only = F) {
  effect_type <- match.arg(effect_type)
  # Generate a post treatment time indicator -- time since treatment.
  pred_series_tib <- pred_series_tib %>%
    mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var))
  # Potentially drop the pretreatment periods.
  if (post_treat_only) {
    pred_series_tib <- pred_series_tib %>%
      dplyr::filter(post_period_t >= 0)
  }
  # Compute the relevant effect using the counterfactual data.
  if (!is.null(counterfac_var)) {
    cf_att_tib <- .ComputeCounterfactualTE(
      pred_series_tib = pred_series_tib,
      effect_type = effect_type, time_var = time_var,
      treat_period_var = treat_period_var,
      outcome_var = outcome_var,
      counterfac_var = counterfac_var
    )
  }
  # Generate a list of observed effects by time-to-treatment.
  effects_to_jk <- pred_series_tib %>%
    dplyr::mutate(
      abs_eff = !!as.name(outcome_var) - !!as.name(pred_var),
      pct_eff = ifelse(!!as.name(pred_var)==0, NA_real_,
                       (!!as.name(outcome_var) / !!as.name(pred_var)) - 1)
    ) %>%
    split(.[["post_period_t"]])

  # Jackknife the relevant effect (e.g. mean TE) for absolute terms.
  jk_comp_abs <- switch(effect_type,
    "mean" = effects_to_jk %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ resample::jackknife(., mean)),
    "median" = effects_to_jk %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ resample::jackknife(., median)),
    "aggregate" = effects_to_jk %>%
      furrr::future_map(~ .[c(outcome_var,pred_var)]) %>%
      furrr::future_map(~ resample::jackknife(., function(x) {
       sum(x[[outcome_var]]) - sum(x[[pred_var]]) }))
    )
  # Jackknife the relevant effect (e.g. mean TE) for percentage terms.
  jk_comp_pct <- switch(effect_type,
    "mean" = effects_to_jk %>%
      furrr::future_map(~ .[["pct_eff"]]) %>%
      furrr::future_map(~ resample::jackknife(., mean)),
    "median" = effects_to_jk %>%
      furrr::future_map(~ .[["pct_eff"]]) %>%
      furrr::future_map(~ resample::jackknife(., median)),
    "aggregate" = effects_to_jk %>%
      furrr::future_map(~ .[c(outcome_var,pred_var)]) %>%
      furrr::future_map(~ resample::jackknife(., function(x) {
        out <- ifelse(sum(x[[pred_var]])==0, NA_real_,
               (sum(x[[outcome_var]]) / sum(x[[pred_var]])) - 1)
        return(out)
      }))
  )
  # Gather the bounds on the confidence intervals and format as tibble.
  post_range <- c(
    min(pred_series_tib$post_period_t),
    max(pred_series_tib$post_period_t)
  )
  jk_comp_ci <- .FormatJackknife(
    jk_est = jk_comp_abs, ci_in = alpha_ci, col_name =
      paste(effect_type, "_abs_att", sep = ""), range = post_range
  ) %>%
    dplyr::inner_join(
      .FormatJackknife(
        jk_est = jk_comp_pct, ci_in = alpha_ci, col_name =
          paste(effect_type, "_pct_att", sep = ""), range = post_range
      ),
      by = "post_period_t"
    )
  # Add the number of treated units per period, to give perspective on bounds.
  jk_ci_by_t <- pred_series_tib %>%
    dplyr::group_by(post_period_t) %>%
    dplyr::summarise(treated_n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(post_period_t, treated_n) %>%
    dplyr::inner_join(jk_comp_ci, by = "post_period_t")
  # Add corresponding counterfactual effects if relevant.
  if (!is.null(counterfac_var)) {
    jk_ci_by_t <- jk_ci_by_t %>% 
      dplyr::left_join(cf_att_tib, by = "post_period_t")
  }
  return(jk_ci_by_t)
}

#' Format jackknifed estimates into a tibble with relevant effect and bounds.
#'
#' @param jk_est Jackknife estimate output in list form per time period.
#' @param ci_in Numeric level specified for the confidence band.
#' @param col_name Name for the type of the desired output effect column.
#' @param range Range of time periods to compute and format effects for.
#'
#' @return Tibble with the jackknifed estimates and confidence intervals
#'    by post treatment time period.
.FormatJackknife <- function(jk_est, ci_in, col_name, range) {
  # Name of the output columns this function is formatting.
  statname <- paste("jackknife_", col_name, sep = "")
  obs_statname <- paste("observed_", col_name, sep = "")
  # Map the jackknife estimates to a confidence interval, then extract the
  # components related to the treatment effect and the confidence intervals.
  # Store these estimates in a tibble, and create columns from the input list.
  jk_est_ci <-
    furrr::future_map(jk_est, CI.t,
      probs = c(0.5 - ci_in / 2, 0.5 + ci_in / 2)
    ) %>%
    do.call(rbind, .) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      "post_period_t" = seq(range[1], range[2], 1),
      !!as.name(statname) := jk_est %>%
        furrr::future_map(~ .[["stats"]]) %>%
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist(),
      !!as.name(obs_statname) := jk_est %>%
        furrr::future_map(~ .[["stats"]]) %>%
        furrr::future_map(~ .[["Observed"]]) %>%
        unlist()
    )
  # Rename the tibble for the output.
  names(jk_est_ci)[1:2] <- c(
    paste("jackknife_lb_", col_name, sep = ""),
    paste("jackknife_ub_", col_name, sep = "")
  )

  return(jk_est_ci)
}
