pacman::p_load(
  dplyr, tibble, CausalImpact, gsynth, tidyr, Matrix, quadprog,
  janitor, stats, readr, bigstatsr, purrr, furrr
)

# TODO(alexdkellogg): How to handle negative point predictions (which mess up
#    computation of percent effects). For synthetic data, just set to 0 as all
#    data is positive. If data can be negative, perhaps just shift up data.
# TODO(alexdkellogg): Rename the output of the estimates. Keep target instead of
#    switching to "response". Get rid of "point.pred" notation, switch to "_"?

#' Estimate the causal impact for a single unit.
#'
#' @param treat_data Individual treated unit time series.
#' @param control_data Matrix of control time series.
#' @param pre_range Vector input of the pretreatment period start and end date
#'    for a given unit.
#' @param post_range Vector input of the posttreatment period start and end date
#'    for a given unit.
#'
#' @return Output from the CausalImpact function for the effect of treatment on
#'    a single treated unit.
# TODO(alexdkellogg): Figure out how to introduce flexibility with ... using
#    arguments that causalimpact can handle.
.CausalImpactHelper <- function(treat_data,
                                control_data,
                                pre_range,
                                post_range) {
  causal_est <- CausalImpact::CausalImpact(
    cbind(treat_data,
          control_data,
          deparse.level = 0
          ),
    pre_range,
    post_range
  )
  return(causal_est)
}

#' Estimate the causalimpact point estimates of the outcome for a full dataset.
#'
#' @param data_full Input data, long form, with each row containing columns for
#'    the treatment status, true outcome, and potentially counterfactual outcome
#'    of a particular unit-time combination. There shouldn't be any missing data
#'    and the format should be a true panel, with each unit having data for each
#'    time period in the data.
#' @param id_var String name for the unit identifier column.
#' @param time_var String name for the time identifier column.
#' @param treat_indicator String name for the treatment indicator column.
#' @param outcome_var String name for the outcome variable.
#' @param counterfac_var String name for the counterfactual outcome
#'    (e.g. potential untreated), can be NULL.
#'
#' @return A long form tibble with the causal impact method estimates of the
#'    outcome variable (missing outcome imputation), the point and percent
#'    effects, and a lower and upper bound from causal impact.
EstimateCausalImpactSeries <- function(data_full,
                                       id_var = "entry",
                                       time_var = "period",
                                       treat_indicator = "treatperiod_0",
                                       outcome_var = "target",
                                       counterfac_var = "counter_factual") {
  # Identify the ids for treated and control entries.
  tr_entries <- data_full %>%
    dplyr::filter(!!as.name(treat_indicator) > 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull() %>%
    sort()
  ct_entries <- setdiff(data_full %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull(), tr_entries)
  # Create a dataframe of the subset of control units.
  cd <- data_full %>% dplyr::filter(!!as.name(id_var) %in% ct_entries)
  cd <- merge(cd, data.frame(entry = ct_entries, rank = seq_along(ct_entries)))
  # Construct the control data matrix, i=time, j=unit, x=outcome.
  donor_outcome_matrix <- as.matrix(Matrix::sparseMatrix(
    x = cd[[outcome_var]],
    i = cd[[time_var]],
    j = cd$rank
  ))
  # Create a tibble of only treated entries.
  td_all <- data_full %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries) %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var))

  # Split into a list of tibbles, one element is one treated unit.
  split_treat_data <- td_all %>% split(.[[id_var]])
  # Grab the outcome variable from each list, store in it's own list.
  list_of_input_data <- split_treat_data %>%
    furrr::future_map(~ .[[outcome_var]])

  # Store the pre treatment range for each treated list into its own list.
  list_of_pretreat_ranges <- split_treat_data %>%
    furrr::future_map(~ .[[treat_indicator]]) %>%
    furrr::future_map(~ range(1, sum(. == 0)))
  # Store the post treatment range for each treated list into it's own list.
  list_of_posttreat_ranges <- split_treat_data %>%
    furrr::future_map(~ .[[treat_indicator]]) %>%
    furrr::future_map(~ range(sum(. == 0) + 1, sum(. == 0) + sum(. == 1)))
  # Map these input lists into parallelized computation of the causal impact
  # estimator using the helper function (wrapper).
  list_helper_out <- furrr::future_pmap(list(
    list_of_input_data,
    list_of_pretreat_ranges,
    list_of_posttreat_ranges
    ),
    .CausalImpactHelper,
    control_data = donor_outcome_matrix
  )
  # Grab the relevant output from causal impact -- the series of predictions.
  list_of_causalimpact_series <- list_helper_out %>%
    furrr::future_map(~ .$series) %>%
    furrr::future_map(~ as.data.frame(.))
  # Create a mapping for each treated id and time period, used as row names.
  period_entry_rowlabs_df <- data.frame(
    tempid = rep(sort(tr_entries),
      each = max(data_full[[time_var]])
    ),
    temp_t = rep(seq(max(data_full[[time_var]])),
      times = length(tr_entries)
    )
  )
  names(period_entry_rowlabs_df) <- c(id_var, time_var)
  # Arrange the causal impact series output list into a data frame.
  # Join with the row labels, and use these to merge into the input data.
  causalimpact_series_output <- list_of_causalimpact_series %>%
    do.call(bind_rows, .) %>%
    dplyr::select(c(
      response, point.pred, point.effect, point.effect.lower,
      point.effect.upper
    )) %>%
    cbind(., period_entry_rowlabs_df, deparse.level = 0) %>%
    dplyr::inner_join(td_all %>%
      dplyr::group_by(!!as.name(id_var)) %>%
      dplyr::mutate(
        Treatment_Period =
          length(!!as.name(treat_indicator)) -
            sum(!!as.name(treat_indicator)) + 1
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(tidyselect::all_of(c(
        id_var, time_var,
        "Treatment_Period"
      ))),
    by = c(id_var, time_var)
    ) %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var)) %>%
    dplyr::mutate(pct.effect = ifelse(point.pred==0, NA, 
                                      (response / point.pred) - 1)) %>%
    tibble::as_tibble()
  # If the input data has counterfactuals, compute counterfactual effects.
  if (!is.null(counterfac_var)) {
    causalimpact_series_output <- causalimpact_series_output %>%
      dplyr::inner_join(
        td_all %>%
          dplyr::select(tidyselect::all_of(c(
            id_var, time_var,
            counterfac_var
          ))),
        by = c(id_var, time_var)
      ) %>%
      dplyr::mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = ifelse(!!as.name(counterfac_var)==0, NA, 
                               (response / !!as.name(counterfac_var)) - 1)
      )
  }
  rownames(causalimpact_series_output) <- NULL
  # Reorder columns for consistency across methods.
  causalimpact_series_output <- causalimpact_series_output %>%
    dplyr::select(
      tidyselect::all_of(c(
        time_var, id_var, "response",
        "point.pred", "point.effect",
        "pct.effect", "Treatment_Period"
      )),
      dplyr::everything()
    )
  return(causalimpact_series_output)
}

#' Estimate the gsynth point estimates of the outcome.
#'
#' @inheritParams EstimateCausalImpact
#' @param ... Additional parameters to pass to gsynth::gsynth function. These
#'    include the estimator (="mc", "ife"), X variables, and se (among others).
#'
#' @return A long form tibble with the gsynth method estimates of the
#'    outcome variable (missing outcome imputation), as well as point and
#'    percent effects.
EstimateGsynthSeries <- function(data_full,
                                 id_var = "entry",
                                 time_var = "period",
                                 treat_indicator = "treatperiod_0",
                                 outcome_var = "target",
                                 counterfac_var = "counter_factual", ...) {
  unnamed_args <- list(...)
  se <- ifelse(is.null(unnamed_args$se), F, unnamed_args$se)
  inference <- ifelse(is.null(unnamed_args$inference),
    "nonparametric", unnamed_args$inference
  )
  # Need to turn parallel off as it's on by default.
  parallel <- ifelse(is.null(unnamed_args$parallel), F, unnamed_args$parallel)
  # Call the gsynth function on the given input.
  gsynth_agg_te_all_t <- gsynth::gsynth(
    Y = outcome_var, D = treat_indicator, data = data_full,
    index = c(id_var, time_var), parallel = parallel, ...
  )
  # Determine whether any uncertainty parameters were called.
  # If user specified se and inference output, provide it.
  if (se & inference == "parametric") {
    # Get all individual effects for all time periods.
    # Renames the multi-dimensional array output for easy pivotting.
    gsynth_indiv_te_series <-
      gsynth_agg_te_all_t$est.ind[, , seq((gsynth_agg_te_all_t$Ntr))] %>%
      tibble::as_tibble() %>%
      dplyr::mutate(!!as.name(time_var) := seq_len(nrow(.))) %>%
      janitor::clean_names() %>%
      stats::setNames(ifelse(stringr::str_count(names(.), "_") > 1,
        stringr::str_replace(names(.), "[_]", ""),
        names(.)
      ))
    # Pivot the wide (one row per time) tibble of effects (and CI) per entry
    # into long data, where each period-entry combination has a column for
    # effect, CI.
    gsynth_series_output <- gsynth_indiv_te_series %>%
      tidyr::pivot_longer(
        -!!as.name(time_var),
        names_to = c(".value", id_var),
        names_sep = "_"
      ) %>%
      dplyr::mutate(!!as.name(id_var) := as.numeric(!!as.name(id_var))) %>%
      dplyr::rename(point.effect = eff)
  }
  else {
    # If no uncertainty estimates were desired, pivot the output without CI.
    gsynth_series_output <-
      (gsynth_agg_te_all_t$Y.tr - gsynth_agg_te_all_t$Y.ct) %>%
      tibble::as_tibble() %>%
      tibble::rownames_to_column(var = time_var) %>%
      tidyr::pivot_longer(
        -!!as.name(time_var),
        names_to = id_var,
        values_to = "point.effect"
      ) %>%
      dplyr::mutate(
        !!as.name(id_var) := as.numeric(!!as.name(id_var)),
        !!as.name(time_var) := as.numeric(!!as.name(time_var))
      )
  }
  # Compute the effects and join the output from gsynth with the original data.
  gsynth_series_output <- gsynth_series_output %>%
    dplyr::inner_join(data_full %>%
      dplyr::group_by(!!as.name(id_var)) %>%
      dplyr::mutate(
        Treatment_Period =
          length(!!as.name(treat_indicator)) -
            sum(!!as.name(treat_indicator)) + 1
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(Treatment_Period)) %>%
      dplyr::select(tidyselect::all_of(c(
        id_var, time_var,
        outcome_var,
        "Treatment_Period"
      ))),
    by = c(id_var, time_var)
    ) %>%
    dplyr::rename(response = outcome_var) %>%
    dplyr::mutate(
      point.pred = response - point.effect,
      pct.effect = ifelse(point.pred==0, NA, 
                          (response / point.pred) - 1)
    )
  # Add counterfactual data if specified.
  if (!is.null(counterfac_var)) {
    gsynth_series_output <- gsynth_series_output %>%
      dplyr::inner_join(
        data_full %>%
          dplyr::select(tidyselect::all_of(c(
            id_var, time_var,
            counterfac_var
          ))),
        by = c(id_var, time_var)
      ) %>%
      dplyr::mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = ifelse(!!as.name(counterfac_var)==0, NA, 
                               (response / !!as.name(counterfac_var)) - 1)
        )
  }
  # Reorder columns for consistency across methods.
  gsynth_series_output <- gsynth_series_output %>%
    dplyr::select(
      tidyselect::all_of(c(
        time_var, id_var, "response",
        "point.pred", "point.effect",
        "pct.effect", "Treatment_Period"
      )),
      dplyr::everything()
    )
  return(gsynth_series_output)
}

#' Estimate the weights for the synthetic control method.
#'
#' @param m_mat Matrix of the relevant control outcomes per period.
#' @param target True outcome of the treated unit per period.
#' @param zeta Penalty term in optimization.
#' @param constrained Boolean indicator for whether the weights should both sum
#'    to one and be non-negative.
#'
#' @return Optimal weights given the control data in m_mat and the data we are
#'    trying to have it equal.
.SyntheticControlWeight <- function(m_mat, target, zeta = 1.0, 
                                    constrained = T) {
  if (nrow(m_mat) != length(target) || length(target) == 0) {
    stop("invalid dimensions")
  }
  
  # Adjust penalty because solve.QP cannot have 0 penalty for quadratic term.
  if (zeta == 0) {
    zeta <- 1e-06
  }
  # Solve a QP with parameters [weights, imbalance]
  # where we use an equality constraint to impose that
  # imbalance = M * weights - target. Objective is encoded as
  # zeta*||weights||^2 + || imbalance ||^2 / length(target)
  # = [ weights, imbalance]' * [zeta*I, 0; 0, (1/length(target)) I]
  # * [weights, imbalance]. In our call to solve.QP, the parameter
  # Dmat is this block-diagonal matrix, and we pass dvec=0 because we
  # have no linear term.
  d_mat <- diag(c(rep(zeta, ncol(m_mat)), rep(1 / length(target), nrow(m_mat))))
  dvec <- rep(0, ncol(m_mat) + nrow(m_mat))
  # Our first nrow(M)+1 constraints are equality constraints.
  # The first nrow(M) impose that M*weights - imbalance = target.
  # The next imposes that sum(weights)=1,
  # and the remaining constraints impose the positivity of our weights.
  
  # Number of equality constrains.
  meq <- nrow(m_mat) + as.numeric(constrained)
  # Each row of a_t is a donor unit, with the constraint that the combination
  # of donors and weights equals the target added column-wise.
  a_t <- cbind(m_mat, diag(1, nrow(m_mat)))
  # Set up the constraints if desired.
  if (constrained) {
    # Add on extra constraints, that weights sum to 1 and are non-negative.
    a_t <- rbind(
      a_t,
      c(rep(1, ncol(m_mat)), rep(0, nrow(m_mat))),
      cbind(diag(1, ncol(m_mat)), matrix(0, ncol(m_mat), nrow(m_mat)))
    )
  }
  # Define the target outcomes (treated unit).
  bvec <- target
  if (constrained) {
    # Add constraints on weights summing to 1 and donor weights (of which there
    # are ncol(m_mat)) non-negative.
    bvec <- c(bvec, 1, rep(0, ncol(m_mat)))
  }
  # Solve the quadratic programming problem for optimal weights.
  gamma <- tryCatch(
    {
      soln <- quadprog::solve.QP(d_mat, dvec, t(a_t), bvec, meq = meq)
      soln$solution[seq_len(ncol(m_mat))]
    },
    error = function(e) {
      # If there is an error, weight each donor equally.
      rep(1 / ncol(m_mat), seq_len(ncol(m_mat)))
    }
  )
  # Return the relevant portion of the solution.
  return(gamma)
}

#' Gets SDID prediction for one treated unit in one period.
#'
#' @param y_mat Matrix of control entries and one treated entry over all
#'  time period prior to the treatment period and one post treatment period
#'  for which the counterfactual prediction is made.
#' @param t_0 Time of treatment period.
#' @param pre_periods Number of periods right before treatment period excluded
#'  from synthetic control.
#' @param post_periods Number of periods past treatment used in counterfactual
#'  prediction. Default = NULL means all periods after are predicted.
#' @inheritParams .SyntheticControlWeight
#' @return a scalar estimate of the counterfactual prediction for one treated
#'  entry and the specified post treatment periods.
.SingleSDIDPredict <- function(y_mat, t_0, pre_periods,
                               post_periods, zeta = var(as.numeric(y_mat)),
                               constrained = T) {
  # Extract the number of entries, and store the last column (treated unit).
  # y_mat is a number of neighbors + treated unit X period matrix.
  nn <- nrow(y_mat)
  pre <- y_mat[nn, ]
  # Use all post periods if not specified.
  if (is.null(post_periods)) {
    post_periods <- ncol(y_mat) - t_0
  }
  # Relevant start and end time periods for prediction.
  end_t <- min(ncol(y_mat), t_0 + post_periods)
  start_t <- max(t_0 - pre_periods, 1)
  # Unit weight per neighbor -- vector of length equal to number of neighbors.
  omega_weight <- .SyntheticControlWeight(t(y_mat[-nn, seq_len(start_t - 1)]),
                                          y_mat[nn, seq_len(start_t - 1)],
                                          zeta = zeta,
                                          constrained = constrained
  )
  # Computation of the time weights, for each post treatment period.
  # Keep time weight constrained for interpretability and avoiding overfitting.
  for (t in seq(start_t, end_t)) {
    # Extract the columns for all our pre-treatment periods and the target 
    # time period, which ranges from the end of the training period to the end
    # of the post-treatment period.
    y_t <- y_mat[, c(seq_len(start_t - 1), t)]
    tt <- ncol(y_t)
    # Find the combination of time weights (from the training period) that best
    # fits the target period, t.
    lambda_weight <- .SyntheticControlWeight(y_t[-nn, -tt], y_t[-nn, tt],
                                             zeta = zeta,
                                             constrained = T
    )
    # Compute the weights * training periods to compute the time-weighted point
    # prediction for the treated unit and target period, t.
    sc_transpose_est <- sum(lambda_weight * y_t[nn, -tt])
    # Compute the weighted combination of donors over all time to form the point
    # prediction using just unit weights.
    sc_est <- sum(omega_weight * y_t[-nn, tt])
    # Find the interaction of the weights.
    interact_est <- omega_weight %*% y_t[-nn, -tt] %*% lambda_weight
    # Compute the final point prediction for treatment in target period t. 
    # Append to vector of treatment unit with training data.
    pre[t] <- sc_est + sc_transpose_est - interact_est
  }
  
  return(pre)
}


#' Main function for computing the SDID predictions of missing outcomes.
#'
#' @param y_df A dataframe of values before treatment, and treatment time.
#' @param ct_mat Matrix with time series entries for control units.
#' @param pre_periods Number of periods excluded from training before the
#'  treatment period.
#' @param post_periods Number of periods for which prediction is done after the
#'  treatment period. When it is NULL, prediction is done for all periods.
#' @param nnsize Nearest neighbor size to be selected for each treated entry.
#'  default value is NULL.
#' @param scale Scaling the entries of the matrix by a constant value to help
#'  the optimization problem as it often fails to encompass large values (so
#'  first scale down to smaller values and then after computation scale
#'  up to large values).
#' @param period Training period, ignores earliest periods in finding nearest 
#'    neighbors based on euclidean distance, to serve as the donors.
#'  @inheritParams .SyntheticControlWeight
#' @return Imputed outcome prediction for the given entry.
.SDIDPredictionHelper <- function(y_df,
                                  ct_mat,
                                  pre_periods = 0,
                                  post_periods = NULL,
                                  nnsize = NULL,
                                  scale = 100,
                                  period = 52,
                                  constrained = T) {
  # TODO(alexkellogg): Consider reducing the duplicated code also in 
  # EstimateSCMSeries.
  # Separate the treatment time and the outcome process.
  treatperiod <- y_df[1, 2]
  stopifnot(treatperiod > 2)
  y <- y_df[seq(1, (treatperiod - 1)), 1]
  # Rescale the outcomes for control and treated.
  stopifnot(scale > 0)
  y_con <- ct_mat / scale
  nc <- ncol(y_con)
  np <- nrow(y_con)
  y_pre <- c(y / scale, rep(0, np - treatperiod + 1))
  
  # tv is a linear weight vector for time. Periods before treatperiod - period
  # use zero weights.
  # nnsize neighbours are identified using weighted distances before the
  # treatment period.
  start_period <- max(1, (treatperiod - 1) - period - 1)
  tv <- c(rep(0, start_period), seq_len(np - start_period))
  # Matrix of the linear period weights, T X N. First start_period rows have
  # 0 weight, remaining rows have weight 1:N.
  tmc <- t(matrix(tv, nrow = nc, ncol = length(tv), byrow = TRUE))
  # Treatment periods have non-zero weight for control, but 0 for treatment.
  wt_y_con <- sqrt(tmc) * y_con
  wt_y_pre <- sqrt(tv) * y_pre
  
  # If nnsize is NULL, the number of the neighbors are chosen to be close
  # to the number of periods.
  if (is.null(nnsize)) {
    nnsize <- treatperiod
  }
  # y_r is a T X N matrix of the time weighted treatment values (0 everywhere 
  # except relevant pre-periods).
  y_r <- t(matrix(wt_y_pre, nrow = nc, ncol = np, byrow = TRUE))
  # Compute the distance between the weighted control and treatment time series 
  # in the relevant pre-periods, and order them by distance.
  e <-
    (y_r[seq_len(treatperiod - 1), ] - wt_y_con[seq_len(treatperiod - 1), ])**2
  
  es <- colSums(e)
  es_order <- order(es)
  # Combine the nearest neighbors and the treated unit; compute the predictions.
  # y_input is a Period X (Num Neighbors + treated unit) matrix.
  y_input <- cbind(y_con[, es_order[seq_len(min(nnsize, nc))]], y_pre)
  pred <- .SingleSDIDPredict(t(y_input), treatperiod,
                             pre_periods = pre_periods,
                             post_periods = post_periods,
                             constrained = constrained
  )
  out <- pred * scale
  return(out)
}

#' Compute the SDID estimates of the outcome series.
#'
#' @inheritParams .SDIDPredictionHelper
#' @inheritParams EstimateCausalImpactSeries
#'
#' @return A long form tibble with the SDID method estimates of the
#'    outcome variable (missing outcome imputation), as well as point and
#'    percent effects.
# TODO(alexdkellogg): Replace/Compare with the 2 way fixed effects, full panel
#    version of SDID as in EQ 2.3 in Athey et al 2019.
EstimateSDIDSeries <- function(data_full,
                               id_var = "entry",
                               time_var = "period",
                               treat_indicator = "treatperiod_0",
                               outcome_var = "target",
                               counterfac_var = "counter_factual",
                               pre_periods = 0,
                               post_periods = NULL,
                               nnsize = NULL,
                               scale = 100,
                               period = 52,
                               constrained = T) {
  # Identify the treated and control ids.
  tr_entries <- data_full %>%
    dplyr::filter(!!as.name(treat_indicator) > 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull() %>%
    sort()
  ct_entries <- setdiff(data_full %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull(), tr_entries)

  # Create control data separate from the treatment.
  control_data <- data_full %>% dplyr::filter(!!as.name(id_var) %in% ct_entries)
  # Number of control entries, just the number of unique entries in cd.
  n0 <- control_data %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    nrow()

  # Create a treated data set.
  treat_data <- data_full %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries) %>%
    dplyr::mutate(new_id = as.integer(as.factor(!!as.name(id_var))) + n0) %>%
    dplyr::arrange(!!as.name(time_var), new_id) %>%
    dplyr::group_by(new_id) %>%
    dplyr::mutate(
      Treatment_Period =
        length(!!as.name(treat_indicator)) -
          sum(!!as.name(treat_indicator)) + 1
    ) %>%
    dplyr::ungroup()

  # Transform the control data into matrix form for passing to helper functions.
  control_matrix <- tidyr::spread(
    control_data %>%
      dplyr::select(tidyselect::all_of(c(time_var, id_var, outcome_var))),
    !!as.name(time_var), !!as.name(outcome_var)
  ) %>%
    dplyr::select(-tidyselect::all_of(id_var)) %>%
    as.matrix() %>%
    t()
  # Split the treatment data into a list of individual treated series.
  split_treat_data <- treat_data %>% split(.[[id_var]])
  # For each treated unit, find and store when it was treated.
  list_of_treat_times <- split_treat_data %>%
    lapply(., function(x) {
      x %>%
        dplyr::select(Treatment_Period) %>%
        dplyr::first()
    })
  # Grab just the outcome process for each treated list.
  list_of_treat_data <- split_treat_data %>% lapply(., function(x) {
    x %>%
      dplyr::select(tidyselect::all_of(outcome_var)) %>%
      as.matrix()
  })
  # Bind the elements of the list together for passing to the helper.
  list_inputs <- Map(cbind, list_of_treat_data, list_of_treat_times)
  list_of_scdid_series <- furrr::future_map(list_inputs, .SDIDPredictionHelper,
    control_matrix,
    pre_periods = pre_periods,
    post_periods = post_periods,
    nnsize = nnsize, scale = scale,
    period = period,
    constrained = constrained
  )

  # Transform the output of the SDID individual estimators (list of estimates)
  # into a single tibble, pivoting back to long form and reformatting.
  scdid_series_tib <- list_of_scdid_series %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = time_var) %>%
    dplyr::mutate(!!as.name(time_var) := as.numeric(!!as.name(time_var))) %>%
    tidyr::pivot_longer(
      cols = starts_with("X"),
      names_to = "temp_id",
      names_transform = list(temp_id = readr::parse_number),
      values_to = "point.pred"
    ) %>%
    dplyr::rename(!!as.name(id_var) := temp_id) %>%
    dplyr::inner_join(
      treat_data %>%
        dplyr::select(tidyselect::all_of(c(
          id_var, time_var, outcome_var,
          "Treatment_Period"
        ))) %>%
        dplyr::rename(response = outcome_var),
      by = c(id_var, time_var)
    ) %>%
    dplyr::mutate(
      point.effect = response - point.pred,
      pct.effect = ifelse(point.pred==0, NA, 
                          (response / point.pred) - 1)
    )
  # Join in the counterfactual data and effects if desired.
  if (!is.null(counterfac_var)) {
    scdid_series_tib <- scdid_series_tib %>%
      dplyr::left_join(
        treat_data %>%
          dplyr::select(tidyselect::all_of(c(
            id_var, time_var,
            counterfac_var
          ))),
        by = c(id_var, time_var)
      ) %>%
      dplyr::mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = ifelse(!!as.name(counterfac_var)==0, NA, 
                               (response / !!as.name(counterfac_var)) - 1)
        )
  }
  # Reorder columns for consistency across methods.
  scdid_series_tib <- scdid_series_tib %>%
    dplyr::select(
      tidyselect::all_of(c(
        time_var, id_var, "response",
        "point.pred", "point.effect",
        "pct.effect", "Treatment_Period"
      )),
      dplyr::everything()
    )
  return(scdid_series_tib)
}

#' Compute the SCM estimates of the outcome series.
#'
#' @inheritParams EstimateCausalImpactSeries
#'
#' @return A long form tibble with the SCM method estimates of the
#'    outcome variable (missing outcome imputation), as well as point and
#'    percent effects.
# TODO(alexdkellogg): Lots of redundant set up code (e.g. first 30 lines here
#    and in SDID). Figure out a way to tie it together, despite later use.
EstimateSCMSeries <- function(data_full, id_var = "entry",
                              time_var = "period",
                              treat_indicator = "treatperiod_0",
                              outcome_var = "target",
                              counterfac_var = "counter_factual") {
  # Identify the ids of the treated and control units.
  tr_entries <- data_full %>%
    dplyr::filter(!!as.name(treat_indicator) > 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull() %>%
    sort()
  ct_entries <- setdiff(data_full %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull(), tr_entries)
  # Create a control dataset.
  control_data <- data_full %>% dplyr::filter(!!as.name(id_var) %in% ct_entries)
  # Create a treated data set.
  treat_data <- data_full %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries) %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var)) %>%
    dplyr::group_by(!!as.name(id_var)) %>%
    dplyr::mutate(
      Treatment_Period =
        length(!!as.name(treat_indicator)) -
          sum(!!as.name(treat_indicator)) + 1
    ) %>%
    dplyr::ungroup()
  # Transform the control data into matrix form for passing to helper functions.
  control_matrix <- tidyr::spread(
    control_data %>%
      dplyr::select(tidyselect::all_of(c(time_var, id_var, outcome_var))),
    !!as.name(time_var), !!as.name(outcome_var)
  ) %>%
    dplyr::select(-tidyselect::all_of(id_var)) %>%
    as.matrix() %>%
    t()
  # Split the treated data, taking only the pretreatment periods for training.
  split_treat_data_pre <- treat_data %>%
    dplyr::filter(!!as.name(time_var) < Treatment_Period) %>%
    dplyr::select(tidyselect::all_of(c(id_var, outcome_var))) %>%
    dplyr::group_by(!!as.name(id_var)) %>%
    dplyr::group_split(.keep = F) %>%
    lapply(., function(x) x %>% dplyr::pull())

  # For each treated unit, find when it was treated.
  list_of_treat_times <- treat_data %>%
    split(.[[id_var]]) %>%
    lapply(., function(x) {
      x %>%
        dplyr::select(Treatment_Period) %>%
        dplyr::slice(1) %>%
        dplyr::pull()
    })
  # Parallel computation of the SCM estimates, unit by unit.
  list_imputed_y <- purrr::map2(
    .x = split_treat_data_pre,
    .y = list_of_treat_times,
    .f = ~ .SCMImpute(
      treat_data = .x,
      control_mat = control_matrix,
      treat_time = .y
    )
  )
  # Combine the list of estimated series into a dataframe, compute effects.
  flex_scm_series <- Map(
    dplyr::bind_cols, treat_data %>% split(.[[id_var]]),
    lapply(list_imputed_y, function(x) {
      return(tibble::tibble("point.pred" = x))
    })
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::rename(response = !!as.name(outcome_var)) %>%
    dplyr::mutate(
      point.effect = response - point.pred,
      pct.effect = ifelse(point.pred==0, NA, 
                          (response / point.pred) - 1)
    )
  # Compute the counterfactual effects if relevant.
  if (!is.null(counterfac_var)) {
    flex_scm_series <- flex_scm_series %>%
      dplyr::mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = ifelse(!!as.name(counterfac_var)==0, NA, 
                               (response / !!as.name(counterfac_var)) - 1)
      )
  }

  # Reorder the columns for consistency across methods.
  flex_scm_series <- flex_scm_series %>%
    dplyr::select(
      tidyselect::all_of(c(
        time_var, id_var, "response",
        "point.pred", "point.effect",
        "pct.effect", "Treatment_Period"
      )),
      dplyr::everything()
    )

  return(flex_scm_series)
}

#' Compute the unconstrained SCM weights for a given treated units.
#'
#' @param treat_data Treated series for a unit, over all time.
#' @param control_mat Matrix of control series, over all time.
#' @param treat_time Time when the treatment is assigned, splits pre and post.
#'
#' @return Predicted series of the outcome for the treated unit (over all time)
#'    based on Elastic Net Penalty trained on the pre period data.
# TODO(alexdkellogg): Change up bigSpLinReg to glmentUtils::cva function
# due to as_FBM() not being compatible with looping in the benchmarking.
.SCMImpute <- function(treat_data, control_mat, treat_time) {
  # Define the control data matrix for the training period.
  control_pre <- control_mat[seq_len(treat_time - 1), ]
  # Fit a speedy Cross-Model Selection and Averaging Grid for ENP.
  fit_enp <- bigstatsr::big_spLinReg(bigstatsr::as_FBM(control_pre),
    treat_data,
    alphas = c(1e-4, 0.2, 0.5, 0.8, 1),
    warn = F
  )
  # summary(fit_enp, best.only=T)
  # Compute the predictions using the best cross-validated model.
  imputed_y <- predict(fit_enp, bigstatsr::as_FBM(control_mat))
  return(imputed_y)
}