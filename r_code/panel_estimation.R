pacman::p_load(dplyr, tibble, CausalImpact,gsynth,tidyr,Matrix, quadprog,
               janitor, stats,readr)


############################################
# Functions for estimating Treatment Effects

# Potential Improvements/Options:
# Could add a units_to_estimate parameter to control directly how mnay units to run through
# Fix: How should we handle negative point predictions when they don't make sense? Set to 0?
# this issue percolates through since, if response>0 all t, % impacts don't mean anything for these sign changers
###########################################


causalimpact_helper <- function(treat_data,
                                control_data,
                                pre_range,
                                post_range) {
  return(CausalImpact::CausalImpact(
    cbind(treat_data,
          control_data,
          deparse.level = 0
    ),
    pre_range,
    post_range
  ))
}

estimate_causalimpact_series <- function(data_full,
                                         id_var = "entry",
                                         time_var = "period",
                                         treat_indicator = "treatperiod_0",
                                         outcome_var = "target",
                                         counterfac_var = "counter_factual") {
  # Estimates CausalImpact treatment effects given a long form data set, outputting a dataframe
  # consisting of a series of treatments effects for each id_var by time_var in all  periods
  
  # Args
  # data_full: long-form dataframe with both treated and control entries, sorted by period and entry (ascending).
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_indicator:column name of binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)
  
  # Higher level description of data_full:
  # rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  # each row should have a treatment indicator (treat_indicator), a period number (time_var),
  # an individual ID (id_var), and an outcome (outcome_var)
  # for associated with that period-ID combination
  
  
  # Output
  # Dataframe containing the full series of the outcome (all T), as well as predicted (counterfactual) outcome
  # and the associated effects by id_Var, time_var
  
  
  # Split the dataset based on whether they are ever treated
  tr_entries <- data_full %>%
    dplyr::filter(!!as.name(treat_indicator) > 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull() %>%
    sort()
  ct_entries <- setdiff(data_full %>%
                          dplyr::distinct(!!as.name(id_var)) %>%
                          dplyr::pull(), tr_entries)
  # Create a dataframe of the subset of control units
  cd <- data_full %>% dplyr::filter(!!as.name(id_var) %in% ct_entries)
  cd <- merge(cd, data.frame(entry = ct_entries, rank = seq_along(ct_entries)))
  
  
  
  # construct control data matrix
  # For all the control data: Row of matrix (index i) indicates the time, Column (j) indicates observation
  # matrix entry represents the Target
  donor_outcome_matrix <- as.matrix(Matrix::sparseMatrix(x = cd[[outcome_var]],
                                                 i = cd[[time_var]],
                                                 j = cd$rank))
  
  
  
  # To output Parallelized process, we gather the inputs as a series of lists
  # first, we want the full treated dataframe to manipulate
  td_all <- data_full %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries) %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var)) # NEW, MIGHT CAUSE ERROR
  
  # Parallelized computation of the data matrix, by id_var, required for causal impact
  # specifically,  append treated data as the first column in the
  # donor_outcome_matrix matrix (described above, col=obs, row=time, entry=outcome)
  # matrix is now TREATED_DATA in col 1 for all of time, and then all control data
  split_treat_data <- td_all %>% split(.[[id_var]])
  
  
  list_of_input_data <- split_treat_data %>%
    furrr::future_map(~ .[[outcome_var]]) # %>% furrr::future_map(~cbind(.,donor_outcome_matrix, deparse.level = 0))
  
  # Parallelized computation of the pre-treatment range, by id_var
  list_of_pretreat_ranges <- split_treat_data %>%
    furrr::future_map(~ .[[treat_indicator]]) %>%
    furrr::future_map(~ range(1, sum(. == 0)))
  # Parallelized computation of the post-treatment range, by id_var
  list_of_posttreat_ranges <- split_treat_data %>%
    furrr::future_map(~ .[[treat_indicator]]) %>%
    furrr::future_map(~ range(sum(. == 0) + 1, sum(. == 0) + sum(. == 1)))
  
  # We want to create a "causal impact helper function" that only takes the control matrix 1, and then binds the particular treated unit to the control matrix
  # and then itself calls causal impact. We will call this function in parallel, and thus only need one 1 repition of the large control matrix per call
  # Parallelized computation of the causal impact
  list_of_causalimpact_series <- furrr::future_pmap(list(
    list_of_input_data,
    list_of_pretreat_ranges,
    list_of_posttreat_ranges
  ),
  causalimpact_helper,
  control_data = donor_outcome_matrix
  ) %>%
    furrr::future_map(~ .$series) %>%
    furrr::future_map(~ as.data.frame(.))

  period_entry_rowlabs_df <- data.frame(
    tempid = rep(sort(tr_entries),
                 each = max(data_full[[time_var]])
    ),
    temp_t = rep(1:max(data_full[[time_var]]),
                 times = length(tr_entries)
    )
  )
  names(period_entry_rowlabs_df) <- c(id_var, time_var)
  
  causalimpact_series_output <- list_of_causalimpact_series %>%
    do.call(bind_rows, .) %>%
    dplyr::select(c(response, point.pred, point.effect, point.effect.lower, point.effect.upper)) %>%
    cbind(., period_entry_rowlabs_df, deparse.level = 0) %>%
    dplyr::inner_join(td_all %>%
                 dplyr::group_by(!!as.name(id_var)) %>%
                 dplyr::mutate(Treatment_Period = length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1) %>%
                 dplyr::ungroup() %>%
                 dplyr::select(id_var, time_var, Treatment_Period), by = c(id_var, time_var)) %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var))
  
  
  if (!is.null(counterfac_var)) {
    causalimpact_series_output <- causalimpact_series_output %>%
      dplyr::inner_join(
        td_all %>%
          dplyr::select(id_var, time_var, counterfac_var),
        by = c(id_var, time_var)
      ) %>%
      dplyr::mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = (response / !!as.name(counterfac_var)) - 1
      )
  }
  
  # add a column with relative (pct) effect
  causalimpact_series_output <- causalimpact_series_output %>% dplyr::mutate(
    pct.effect = (response / point.pred) - 1
  )
  
  rownames(causalimpact_series_output) <- NULL
  
  
  
  
  return(causalimpact_series_output)
}





# Potential Improvement: can I use ... as an argument so that any other gsynth arguments (ex, k=5) can be passed? right now, no way to do that for user
estimate_gsynth_series <- function(data_full, id_var = "entry", time_var = "period", treat_indicator = "treatperiod_0", outcome_var = "target", x_in = NULL,
                                   counterfac_var = "counter_factual", ...) {
  # Estimates Gsynth treatment effects given a long form data set, outputting a dataframe
  # consisting of a series of treatments effects for each id_var by time_var in all periods
  
  # Args
  # data_full: long-form dataframe with both treated and control entries, sorted by period and entry (ascending).
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # outcome_var: the y var for the time series
  # treat_indicator:column name of binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)
  # x_in: time-varying covariates
  # note: a panel reg formula is implicitly defined as outcome_var~treat_indicator+x_in
  # for the remaining parameters, these are options into Gsynth (use ?gsynth for info)
  
  
  # Higher level description of data_full:
  # rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  # each row should have a treatment indicator (treat_indicator), a period number (time_var),
  # an individual ID (id_var), and an outcome (outcome_var)
  # for associated with that period-ID combination
  
  
  # Output
  
  # estimate the panel SCM
  gsynth_agg_te_all_t <- gsynth::gsynth(
    Y = outcome_var, D = treat_indicator, data = data_full, 
    index = c(id_var, time_var), X = x_in, ...
  )
  
  unnamed_args=list(...)
  #default in gsynth is F
  se=ifelse(is.null(unnamed_args$se), F,unnamed_args$se )
  #Default in GSynth is nonparametric
  inference=ifelse(is.null(unnamed_args$inference),
                   "nonparametric",unnamed_args$inference )
  
  if (se & inference == "parametric") {
    # gets all indiv effects for all T
    # renames the multi-d array output so that we can more easily pivot it
    gsynth_indiv_te_series <- gsynth_agg_te_all_t$est.ind[, , 1:(gsynth_agg_te_all_t$Ntr)] %>%
      tibble::as_tibble() %>%
      dplyr::mutate(!!as.name(time_var) := seq_len(nrow(.))) %>%
      janitor::clean_names() %>%
      stats::setNames(ifelse(stringr::str_count(names(.), "_") > 1,
                      stringr::str_replace(names(.), "[_]", ""),
                      names(.)))
    
    
    # Pivot the wide dataframe of effects (and CI) per entry with one row per period into a long data set
    # where each period-entry combination has a column for effect, CI.
    # Then, for each period, summarize over all the entries
    gsynth_series_output <- gsynth_indiv_te_series %>%
      tidyr::pivot_longer(
        -!!as.name(time_var),
        names_to = c(".value", id_var),
        names_sep = "_"
      ) %>%
      dplyr::mutate(!!as.name(id_var) := as.numeric(!!as.name(id_var))) %>% 
      dplyr::inner_join(data_full %>%
                   dplyr::group_by(!!as.name(id_var)) %>% 
                   dplyr::mutate(Treatment_Period = length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1) %>%
                   dplyr::ungroup() %>%
                   dplyr::filter(!is.na(Treatment_Period)) %>%
                   dplyr::select(id_var, time_var, outcome_var, Treatment_Period), by = c(id_var, time_var)) %>%
      dplyr::rename(response = outcome_var, point.effect = eff) %>%
      dplyr::mutate(point.pred = response - point.effect)
  }
  else {
    gsynth_series_output <- (gsynth_agg_te_all_t$Y.tr - gsynth_agg_te_all_t$Y.ct) %>%
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
      ) %>%
      dplyr::inner_join(data_full %>%
                   dplyr::group_by(!!as.name(id_var)) %>% 
                   dplyr::mutate(Treatment_Period = length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1) %>%
                   dplyr::ungroup() %>%
                   dplyr::filter(!is.na(Treatment_Period)) %>%
                   dplyr::select(id_var, time_var, outcome_var, Treatment_Period), by = c(id_var, time_var)) %>%
      dplyr::rename(response = outcome_var) %>%
      dplyr::mutate(point.pred = response - point.effect)
  }
  
  if (!is.null(counterfac_var)) {
    gsynth_series_output <- gsynth_series_output %>%
      dplyr::inner_join(
        data_full %>%
          dplyr::select(id_var, time_var, counterfac_var),
        by = c(id_var, time_var)
      ) %>%
      dplyr::mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = (response / !!as.name(counterfac_var)) - 1
      )
  }
  
  # add a column with relative (pct) effect
  gsynth_series_output <- gsynth_series_output %>% dplyr::mutate(
    pct.effect = (response / point.pred) - 1
  )
  
  return(gsynth_series_output)
}






synthetic_control_weight <- function(m_mat, target, zeta = 1) {
  if (nrow(m_mat) != length(target)) {
    stop("invalid dimensions")
  }
  
  # solve.QP cannot have 0 penalty for quadratic term
  if (zeta == 0) {
    zeta <- 1e-06
  }
  # we solve a QP with parameters [weights, imbalance]
  # where we use an equality constraint to impose that
  # imbalance = M * weights - target. Our objective is encoded as
  # zeta*||weights||^2 + || imbalance ||^2 / length(target)
  # = [ weights, imbalance]' * [zeta*I, 0; 0, (1/length(target)) I]
  # * [weights, imbalance] in our call to solve.QP, the parameter
  # Dmat is this block-diagonal matrix, and we pass dvec=0 because we
  # have no linear term
  
  d_mat <- diag(c(rep(zeta, ncol(m_mat)), rep(1 / length(target), nrow(m_mat))))
  dvec <- rep(0, ncol(m_mat) + nrow(m_mat))
  
  # our first nrow(M)+1 constraints are equality constraints
  # the first nrow(M) impose that M*weights - imbalance = target
  # the next imposes that sum(weights)=1
  # and the remaining constraints impose the positivity of our weights
  
  meq <- nrow(m_mat) + 1
  a_t <- rbind(
    cbind(m_mat, diag(1, nrow(m_mat))),
    c(rep(1, ncol(m_mat)), rep(0, nrow(m_mat))),
    cbind(diag(1, ncol(m_mat)), matrix(0, ncol(m_mat), nrow(m_mat)))
  )
  bvec <- c(target, 1, rep(0, ncol(m_mat)))
  soln <- quadprog::solve.QP(d_mat, dvec, t(a_t), bvec, meq = meq)
  gamma <- soln$solution[seq_len(ncol(m_mat))]
  
  return(gamma)
}



#' Gets SDID prediction for one treated advertiser and one post treatment
#' period.
#'
#' @param Y the matrix of control entries and one treated entry over all
#'  time period prior to the treatment period and one post treatment period
#'  for which the counterfactual prediction is made.
#' @param T_0 the treatment period.
#' @param pre_periods Number of periods right before treatment period exlcuded
#'  from sunthetic control.
#' @param post_periods Number of periods past adoption used in counterfactual
#'  prediction. Default = NULL. means all periods after are predicted.
#' @return a scalar estimate of the counterfactual prediction for one treated
#'  entry and the specified post treatment periods.
#' @export

scdid_predict <- function(y_mat, t_0, pre_periods, post_periods, zeta = var(as.numeric(y_mat))) {
  
  # The unit weights are only estimated once, but time weights are estimated
  # for each period.
  nn <- nrow(y_mat)
  pre <- y_mat[nn, ]
  
  if (is.null(post_periods)) {
    post_periods <- ncol(y_mat) - t_0
  }
  
  end_t <- min(ncol(y_mat), t_0 + post_periods)
  start_t <- max(t_0 - pre_periods, 1)
  
  
  omega_weight <- synthetic_control_weight(t(y_mat[-nn, seq_len(start_t - 1)]), y_mat[nn, seq_len(start_t - 1)], zeta = zeta)
  
  for (t in start_t:end_t) {
    y_t <- y_mat[, c(seq_len(start_t - 1), t)]
    tt <- ncol(y_t)
    lambda_weight <- synthetic_control_weight(y_t[-nn, -tt], y_t[-nn, tt], zeta = zeta)
    sc_transpose_est <- sum(lambda_weight * y_t[nn, -tt])
    sc_est <- sum(omega_weight * y_t[-nn, tt])
    interact_est <- omega_weight %*% y_t[-nn, -tt] %*% lambda_weight
    pre[t] <- sc_est + sc_transpose_est - interact_est
  }
  
  return(pre)
}



#' Main function for predictions by NSDID.
#'
#' @param y a vector of history values before treatment.
#' @param ct_mat matrix with time series entries
#' @param treatperiod treatment period.
#' @param pre_periods number of periods excluded from training before the
#'  treatment period.
#' @param post_periods number of periods for which prediction is done after the
#'  treatment period. When it is NULL, prediction is done for all periods.
#' @param nnsize Nearest neighbour size to be selected for each treated entry.
#'  default value is NULL.
#' @param scale scaling the entries of the matrix by a constant value to help
#'  the optimization problem as it often fails to encompass large values (so
#'  first scale down to smaller values and then after computation scale
#'  up to large values).
#' @param period training period as history before treatperiod - period is
#'  ignored.
#' @return prediction for the given entry.
#' @export

nsdid_prediction <- function(y_df, 
                             ct_mat,
                             pre_periods = 0,
                             post_periods = 20,
                             nnsize = NULL, 
                             scale = 100, 
                             period = 52) {
  
  treatperiod <- y_df[1, 2]
  y <- y_df[1:(treatperiod-1), 1]
  
  y_con <- ct_mat / scale
  nc <- ncol(y_con)
  np <- nrow(y_con)
  y_pre <- c(y / scale, rep(0, np - treatperiod + 1))
  
  # tv is a linear weight vector. Periods before treatperiod - period use
  # zero weights.
  # nnsize neighbours are identified using weighted distances before the
  # treatment period.
  
  start_period <- max(1, (treatperiod - 1) - period - 1)
  tv <- c(rep(0, start_period), seq_len(np - start_period))
  tmc <- t(matrix(tv, nrow = nc, ncol = length(tv), byrow = TRUE))
  wt_y_con <- sqrt(tmc) * y_con
  wt_y_pre <- sqrt(tv) * y_pre
  
  
  
  # If nnsize is NULL, the number of the neighbours are chosen to be close
  # to the number of periods.
  
  if (is.null(nnsize)) {
    nnsize <- treatperiod
  }
  
  y_r <- t(matrix(wt_y_pre, nrow = nc, ncol = np, byrow = TRUE))
  e <- (y_r[seq_len(treatperiod - 1), ] - wt_y_con[seq_len(treatperiod - 1), ])**2
  
  es <- colSums(e)
  es_order <- order(es)
  
  y_input <- cbind(y_con[, es_order[seq_len(min(nnsize, nc))]], y_pre)
  pred <- scdid_predict(t(y_input), treatperiod, pre_periods = pre_periods, post_periods = post_periods)
  return(pred * scale)
}




estimate_scdid_series <- function(data_full,
                                  id_var = "entry",
                                  time_var = "period",
                                  treat_indicator = "treatperiod_0",
                                  outcome_var = "target",
                                  counterfac_var = "counter_factual",
                                  pre_sdid = 0,
                                  post_sdid = NULL,
                                  nn_sdid = NULL,
                                  scale_sdid = 100,
                                  period_sdid = 30) {
  # Estimates SCDID treatment effects given a long form data set, outputting a dataframe
  # consisting of a series of treatments effects for each id_var by time_var in all  periods
  
  # Args
  # data_full: long-form dataframe with both treated and control entries, sorted by period and entry (ascending).
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_indicator:column name of binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)
  # outcome_var: the y var for the time series
  # counterfac_var: the counterfactual value for the time series, if available (otherwise, NULL)
  # remaining values are from nsdid_prediction function
  # pre_sdid: number of periods excluded from training before the treatment period.
  # post_sdid: number of periods for which prediction is done after the treatment period. When it is NULL, prediction is done for all periods.
  # nn_sdid: Nearest neighbour size to be selected for each treated entry. default value is NULL.
  # scalesdid: scaling the entries of the matrix by a constant value to help the optimization problem as it often fails to
  # encompass large values (so first scale down to smaller values and then after computation scale up to large values).
  # period_sdid: number determining how many pre-treat periods to weight in the estimation (for synthetic control)
  
  # Higher level description of data_full:
  # rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  # each row should have a treatment indicator (treat_indicator), a period number (time_var),
  # an individual ID (id_var), and an outcome (outcome_var)
  # for associated with that period-ID combination
  
  
  # Output
  # Dataframe containing the full series of the outcome (all T), as well as predicted (counterfactual) outcome
  # and the associated effects by id_Var, time_var
  
  # Split the dataset based on whether they are ever treated
  tr_entries <- data_full %>%
    dplyr::filter(!!as.name(treat_indicator) > 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull() %>%
    sort()
  ct_entries <- setdiff(data_full %>% dplyr::distinct(!!as.name(id_var)) %>% dplyr::pull(), tr_entries)
  
  
  # create control data frame, with a new id for the sake of ordering observations later
  control_data <- data_full %>% dplyr::filter(!!as.name(id_var) %in% ct_entries)
  # n0, number of control entries, is just the number of unique entries in cd
  n0 <- control_data %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    nrow()
  
  # In the loop, we also want to compute the SCDID Estimates as they require a single observation at a time
  # Because SCDID cannot handle staggered adoption, introduce one treated unit at a time
  # THIS GOES WITHIN Causal Impact For Loop
  treat_data <- data_full %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries) %>%
    dplyr::mutate(new_id = as.integer(as.factor(!!as.name(id_var))) + n0) %>%
    dplyr::arrange(!!as.name(time_var), new_id) %>%
    dplyr::group_by(new_id) %>%
    dplyr::mutate(Treatment_Period = length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1) %>%
    dplyr::ungroup()
  

  # create the control matrix once, which is an input to sdid estimator
  control_matrix <- tidyr::spread(
    control_data %>% dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(outcome_var)),
    !!as.name(time_var), !!as.name(outcome_var)
  ) %>%
    dplyr::select(-!!as.name(id_var)) %>%
    as.matrix() %>%
    t()
  
  split_treat_data <- treat_data %>% split(.[[id_var]])

  # for each treated unit, find when it was treated
  list_of_treat_times <- split_treat_data %>%
    lapply(., function(x) {
      x %>%
        dplyr::select(Treatment_Period) %>%
        dplyr::first()
    })

  # split the treated units into individual vectors
  list_of_treat_data <- split_treat_data %>% lapply(., function(x) {
    x %>%
      dplyr::select(!!as.name(outcome_var)) %>%
      as.matrix()
  })

  list_inputs <- Map(cbind, list_of_treat_data, list_of_treat_times)


  list_of_scdid_series <- furrr::future_map(list_inputs, nsdid_prediction, control_matrix, pre_sdid, post_sdid, nn_sdid, scale_sdid, period_sdid)

  
  # compute the TE by subtracting the matrix of predictions (list_of_scdid_series) from the outcome_var
  # Reformat so that the output is the same as the other functions, namely,
  # each row represents a period-unit combination, with the outcome_var, the prediction, the effect (and the treatment time/indicator)
  df_scdid_series <- list_of_scdid_series %>%
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
        dplyr::select(!!as.name(id_var), !!as.name(time_var), !!as.name(outcome_var), Treatment_Period) %>%
        dplyr::rename(response = outcome_var),
      by = c(id_var, time_var)
    ) %>%
    dplyr::mutate(point.effect = response - point.pred)
  
  
  if (!is.null(counterfac_var)) {
    df_scdid_series <- df_scdid_series %>%
      dplyr::left_join(
        treat_data %>%
          dplyr::select(id_var, time_var, counterfac_var),
        by = c(id_var, time_var)
      ) %>%
      dplyr::mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = (response / !!as.name(counterfac_var)) - 1
      )
  }
  
  # add a column with relative (pct) effect
  df_scdid_series <- df_scdid_series %>% dplyr::mutate(
    pct.effect = (response / point.pred) - 1
  )
  
  return(df_scdid_series)
}