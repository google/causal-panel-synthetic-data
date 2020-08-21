pacman::p_load(dplyr, furrr, tidyr, resample, tibble)

#' Computes the bias of the ATT.
#'
#' @param tot_list List of tibbles, each tibble provides the estimated ATT from
#'    a given simulated sample from the same DGP. Will work if the input is just
#'    a single tibble, but this represents just one noisy estimate.
#' @param horizon Range of periods over which to display the output.
#' @param time_var String name for the post period column in the input tibble.
#' @param pred_var String name for the outcome to compute Bias over (e.g. ATT).
#' @param counterfac_var String name of the true outcome.
#'
#' @return A tibble with the period, the absolute bias and the percent bias
#'   estimated for that period over each DGP draw.
ComputeBias <- function(tot_list, horizon = c(0,4),
                        time_var="post_period_t",
                        pred_var="observed_mean_abs_tot",
                        counterfac_var="mean_abs_cf_tot") {
  stopifnot(length(horizon)==2 & horizon[1]<= horizon[2] )
  # Keep only the periods within the horizon.
  tot_tib <- tot_list %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!!as.name(time_var) <= horizon[2],
                  !!as.name(time_var) >= horizon[1])
  # For each period compute the average observed ATT less the average true ATT.
  bias_tib <- tot_tib %>%
    dplyr::group_by(!!as.name(time_var)) %>%
    dplyr::summarise(
      "bias" = mean(!!as.name(pred_var), na.rm = T) -
        mean(!!as.name(counterfac_var), na.rm = T)
    )
  
  return(bias_tib)
}

#' Computes the variance of the ATT.
#'
#' @inheritParams ComputeBias
#'
#' @return A tibble with the period, the variance of the absolute ATT and the
#'    variance of the percent ATT in that period.
ComputeVariance <- function(tot_list, horizon = c(0,4),
                            time_var="post_period_t",
                            pred_var="observed_mean_abs_tot") {
  stopifnot(length(horizon)==2 & horizon[1]<= horizon[2] )
  # Keep only the periods within the horizon.
  tot_tib <- tot_list %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!!as.name(time_var) <= horizon[2],
                  !!as.name(time_var) >= horizon[1])
  # For each period compute the variance of the att estimate.
  var_tib <- tot_tib %>%
    dplyr::group_by(!!as.name(time_var)) %>%
    dplyr::summarise("var_att" = stats::var(!!as.name(pred_var), na.rm = T))
  
  return(var_tib)
}

#' Computes the coverage of the ATT.
#' @param tot_list List of tibbles, each tibble provides the estimated ATT from
#'    a given simulated sample from the same DGP. Will run if the input is a
#'    single tibble, but coverage will vary wildly (either at 0% or 100%).
#' @inheritParams ComputeBias
#'
#' @return A tibble with the period and the estimate coverage of the ATT CI for
ComputeCoverage <- function(tot_list, horizon = c(0,4),
                            time_var="post_period_t",
                            counterfac_var="mean_abs_cf_tot",
                            ub_var="jackknife_ub_mean_abs_tot",
                            lb_var= "jackknife_lb_mean_abs_tot") {
  stopifnot(length(horizon)==2 & horizon[1]<= horizon[2] )
  # Keep only the periods within the horizon.
  tot_tib <- tot_list %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!!as.name(time_var) <= horizon[2],
                  !!as.name(time_var) >= horizon[1])
  # For each period, see how many of the sample CIs contain the truth.
  coverage_out <- tot_tib %>%
    dplyr::group_by(!!as.name(time_var)) %>%
    dplyr::mutate(
      contains_truth = as.numeric(!!as.name(counterfac_var) <=
                                    !!as.name(ub_var) &
                                    !!as.name(counterfac_var) >=
                                    !!as.name(lb_var))) %>%
    dplyr::summarise(
      coverage = mean(contains_truth),
    )
  
  return(coverage_out)
}
#' Compute the metrics associated with the ATT.
#'
#' @inheritParams ComputeBias
#'
#' @return A tibble with the RMSE and MAE of the ATT over DGP seeds,
#'    computed per period, for both absolute and percent ATT.
ComputeMetrics <- function(tot_list, horizon = c(0,4),
                           time_var="post_period_t",
                           pred_var="observed_mean_abs_tot",
                           counterfac_var="mean_abs_cf_tot") {
  stopifnot(length(horizon)==2 & horizon[1]<= horizon[2] )
  # Keep only the periods within the horizon.
  tot_tib <- tot_list %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!!as.name(time_var) <= horizon[2],
                  !!as.name(time_var) >= horizon[1])
  # For each period compute the mse and mae of the att estimate.
  metrics_tib <- tot_tib %>%
    dplyr::group_by(!!as.name(time_var)) %>%
    dplyr::summarise(
      "rmse_att" = sqrt(
        mean((!!as.name(pred_var) - !!as.name(counterfac_var))^2,
             na.rm = T
        )
      ),
      "mae_att" = mean(
        abs(!!as.name(pred_var) - !!as.name(counterfac_var)),
        na.rm = T
      )
    )
  
  return(metrics_tib)
}
