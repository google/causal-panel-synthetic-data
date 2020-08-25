pacman::p_load(dplyr, ggplot2, ggpubr)
#################################################################
# Plots
#################################################################

#' Create a Gap Plot of the ATT over time, along with CI bounds.
#'
#' ATT = Average Treatment Effect on the Treated. Stand-in for other aggregation
#'     methods, such as median or aggregate treatment effect on the treated.
#' CI = Confidence interval, defaults to 95%. 
#' @param att_tib Tibble of ATTs for a single dataset, with CI bounds.
#' @param time_var String name for the time relative to treatment column.
#' @param horizon Number range for the x axis of the plot.
#' @param pred_var String name for the main outcome to plot (e.g. ATT estimate).
#' @param counterfac_var String name for the true outcome (e.g. true ATT).
#' @param ub_var String name for the upper CI bound of the main outcome.
#' @param lb_var String name for the lower CI bound of the main outcome.
#' @param plot_title Optional string for the title of the plot.
#' @param plot_y_lab Optional string for the name of "pred_var" as y-axis label.
#'
#' @return A gap plot, with time to treatment on the x-axis and the outcome of
#'    interest on the y-axis, typically the ATT.
GapPlot <- function(att_tib, time_var = "post_period_t", horizon = c(-10, 30),
                    pred_var = "observed_mean_abs_att",
                    counterfac_var = "mean_abs_cf_att",
                    ub_var = "jackknife_ub_mean_abs_att",
                    lb_var = "jackknife_lb_mean_abs_att",
                    plot_title = NULL,
                    plot_y_lab = NULL) {
  stopifnot(length(horizon) == 2 & horizon[1] < horizon[2])
  # Keep the relevant columns and time periods from the outcome tibble.
  att_tib <- att_tib %>%
    dplyr::select(tidyselect::all_of(
      c(time_var, pred_var, ub_var, lb_var, counterfac_var)
    )) %>%
    dplyr::filter(post_period_t <= horizon[2], post_period_t >= horizon[1])
  # Create a base plot with the outcome of interest and the CI bands.
  plot_out <- att_tib %>%
    ggplot2::ggplot(aes(
      x = !!as.name(time_var), y = !!as.name(pred_var),
      color = "Estimate"
    )) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(aes(
      ymin = !!as.name(lb_var),
      ymax = !!as.name(ub_var)
    ),
    alpha = 0.3, color = NA
    ) +
    ggplot2::ggtitle(ifelse(is.null(plot_title), paste("Gap Plot of", pred_var),
      plot_title
    )) +
    ggplot2::labs(
      x = "Time Relative to Treatment",
      y = ifelse(is.null(plot_y_lab), pred_var, plot_y_lab)
    )
  # Add the true outcome of interest if it's available.
  if (!is.null(counterfac_var)) {
    plot_out <- plot_out +
      ggplot2::geom_line(aes(y = !!as.name(counterfac_var), color = "True"))
  }
  # Finalize some formatting from output.
  plot_out <- plot_out +
    ggplot2::scale_colour_manual(
      name = "",
      values = c("Estimate" = "black", "True" = "red"),
      labels = c("Estimate", "True")
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "blue")
  ggplot2::theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

  return(plot_out)
}

#' Plot the prediction and truth for a grid of 6 individual time series.
#'
#' @param pred_series_tib Tibble of predicted and true series, long form. That
#'    is, each row should represent a unit in a particular time, and have the
#'    data on their true outcome as well as a prediction from some method. 
#'    Since estimators will typically fail with missing data, there should be 
#'    no missing data in this post-estimation input.
#' @param time_var String name of the time column, from 1 to max period.
#' @param id_var String name of the id column.
#' @param outcome_var String name of the true outcome variable.
#' @param pred_var String name of the predicted outcome from some method.
#' @param treat_period_var String name for the treatment time.
#' @param method_name String name for the method used (e.g. "gsynth", "sdid").
#'
#' @return A grid of plots for the 3 biggest and 3 additional random time series
#'    in the data, displaying the outcome process and its prediction.
AssortedSeriesPlot <- function(pred_series_tib, time_var = "period", 
                               id_var = "entry",
                               outcome_var = "response",
                               pred_var = "point.pred",
                               treat_period_var = "Treatment_Period",
                               method_name = NULL) {
  # Get the ids for the top 3 units by their first observed period size (by id). 
  big_ids <- pred_series_tib %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var)) %>%
    dplyr::distinct(!!as.name(id_var), .keep_all=T) %>%
    dplyr::slice_max(n = 3, order_by = !!as.name(outcome_var)) %>%
    dplyr::select(tidyselect::all_of(c(id_var, treat_period_var)))
  # Get the ids for the 3 random units.
  random_ids <- pred_series_tib %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var)) %>%
    dplyr::distinct(!!as.name(id_var), .keep_all=T) %>%
    dplyr::filter(!!as.name(id_var) %in%
      setdiff(!!as.name(id_var), big_ids %>%
        dplyr::pull(!!as.name(id_var)))) %>%
    dplyr::slice_sample(n = 3) %>%
    dplyr::select(tidyselect::all_of(c(id_var, treat_period_var)))
  # Store these ids as vectors for passing into a helper function for plotting.
  indiv_ids <- big_ids %>%
    dplyr::bind_rows(random_ids) %>%
    dplyr::pull(!!as.name(id_var))
  # Call helper to plot the 6 ids.
  plot_out <- furrr::future_map(
    .x = indiv_ids,
    .f = ~ IndividualPlotter(
      pred_series_tib = pred_series_tib, id_inp = .x,
      time_var = time_var,
      id_var = id_var, outcome_var = outcome_var,
      pred_var = pred_var, base_title = F
    )
  )
  # Arrange into a grid for output.
  grid_out <- ggpubr::ggarrange(
    plotlist = plot_out, ncol = 2, nrow = 3,
    common.legend = TRUE, legend = "bottom"
  )
  grid_out <-
    ggpubr::annotate_figure(grid_out, top = text_grob(
      paste(method_name, "Individual Plots")
    ))
  return(grid_out)
}

#' Creates prediction and outcome plot for a given unit.
#'
#' @inheritParams IndividualSeriesPlot
#' @param id_inp Numeric ID for the particular series to plot.
#' @param base_title Boolean indicating whether the output should contain a
#'    default base title describing the plot. Set False when using with
#'    AssortedSeriesPlot.
#'
#' @return A ggplot for the particular time series, formatted.
IndividualPlotter <- function(pred_series_tib, id_inp,
                              time_var = "period", id_var = "entry",
                              outcome_var = "response",
                              pred_var = "point.pred",
                              treat_period_var = "Treatment_Period",
                              base_title = T) {
  # Identify the treatment time for the given unit.
  treat_time <- pred_series_tib %>%
    dplyr::filter(!!as.name(id_var) == id_inp) %>%
    dplyr::pull(!!as.name(treat_period_var))
  # Stop if there is no treat time: covers ID not in the data.
  if (identical(treat_time, numeric(0))) {
    stop("No such ID with a treatment period.")
  }
  # Identify the true and predicted time series of the ID and plot.
  plot_out <- pred_series_tib %>%
    dplyr::filter(!!as.name(id_var) == id_inp) %>%
    ggplot2::ggplot(aes(x = !!as.name(time_var), y = !!as.name(outcome_var))) +
    ggplot2::geom_line(aes(color = "Outcome")) +
    ggplot2::geom_line(aes(y = !!as.name(pred_var), color = "Predicted")) +
    ggplot2::geom_vline(xintercept = max(treat_time), color = "blue") +
    ggplot2::scale_color_manual(name = "Type", 
                                values = c("black", "darkgray")) +
    ggplot2::theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = "Time", y = "Outcome",
      title = ifelse(base_title,
        paste(
          "Counterfactual vs Outcome Series, ID=",
          id_inp
        ),
        paste("ID=", id_inp)
      )
    )
  return(plot_out)
}