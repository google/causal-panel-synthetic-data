pacman::p_load(dplyr, furrr, tidyr, stats,
               caret, glmnet, bigstatsr)







#' Using nearest neighbour method to do the prediction.
#'
#' @param y a vector of history values.
#' @param ct.m sparse matrix of all the histories.
#' @param treatperiod treat period of the given entry.
#' @param neighbour.method method for prediction once the list of neighbours
#'   are found. "median" skips zero entries when computing the median;
#'   "mean" used mean of the neighbours, "none" disables scaling.
#' @param neighbour.weighting boolean value of weighting option. If "TRUE",
#'   more recent history is weighted higher using linear weight. If
#'   "FALSE", no weight is used.
#' @param scaling.option if "centered", we use the standard scaling to treat
#'   histories: (x-mean(x))/std(x); if "mean", we use x/mean(x). If
#'   "none", no scaling is used.
#' @param neighbour.periods maximum number of periods used for prediction.
#' @return A vector with prediction values that also include fit of
#'   pre-treatment periods.
#' @export
NearestNeighbourPrediction <- function(y,
                                       ct.m,
                                       treatperiod,
                                       num.neighbours = 8,
                                       neighbour.method =
                                         c("median", "mean"),
                                       neighbour.weighting = TRUE,
                                       scaling.option =
                                         c("centered", "mean", "none"),
                                       neighbour.periods = 20,
                                       allow.zero.start = FALSE,
                                       trend.upper.bound = 20) {
  neighbour.method <- match.arg(neighbour.method)
  scaling.option <- match.arg(scaling.option)
  
  # Start period starts from non-zero period.
  if (!allow.zero.start) {
    start.period <- min(which(y > 0))
  } else {
    start.period <- 1
  }
  if (length(y) - start.period > neighbour.periods) {
    start.period <- length(y) - neighbour.periods
  }
  y <- y[seq(from = start.period, to = length(y), by = 1)]
  if (sum(y) == 0) {
    pred <- rep(0, nrow(ct.m))
    return(pred)
  }
  x <- ct.m[seq(from = start.period, to = (treatperiod - 1), by = 1), ]
  
  # Make sure there is no zero-sum rows
  col.sumx <- colSums(x)
  x <- x[, which(col.sumx > 0)]
  ct.m <- ct.m[, which(col.sumx > 0)]
  
  # Scale the histories if needed.
  yx <- cbind(y, x)
  if (scaling.option == "centered") {
    yx <- scale(cbind(y, x))
  } else if (scaling.option == "mean") {
    yx <- t(t(yx) / colSums(yx)) * nrow(yx)
  }
  yx[, which(is.na(yx[1, ]))] <- 0
  
  # Set the weights for distance calculation.
  if (neighbour.weighting) {
    weights <- seq_len(nrow(yx))
  } else {
    weights <- rep(1, nrow(yx))
  }
  distance <- colSums(((yx[, seq(from = 2, to = ncol(yx), by = 1)]
                        - yx[, 1]) * weights)^2)
  
  # id <- order(distance)[seq_len(num.neighbours)]
  # Modification: We choose 1.5 times as many neighbours but only use those
  # with closest volume
  id.init <- order(distance)[seq_len(floor(num.neighbours * 1.5))]
  xx <- ct.m[seq(from = start.period, to = (treatperiod - 1), by = 1), id.init]
  yxx <- cbind(y, xx)
  distance <- colSums((yxx[, seq(from = 2, to = ncol(yxx), by = 1)]
                       - yxx[, 1])^2)
  id <- id.init[order(distance)[seq_len(num.neighbours)]]
  
  # All the secret sauces
  if (scaling.option == "centered") {
    z.var <- apply(x[, id], 2, var)
    z.id <- which(z.var == 0)
    if (length(z.id) > 0) z.var[which(z.var == 0)] <- 1
    z <- t(t(ct.m[, id] - matrix(rep(colMeans(x[, id]),
                                     each = nrow(ct.m)
    ), nrow = nrow(ct.m))) /
      sqrt(z.var))
  } else if (scaling.option == "mean") {
    deno <- colSums(ct.m[seq(
      from = start.period, to = (treatperiod - 1),
      by = 1
    ), id])
    # Apply an upper bound on the trends to prevent blowup
    z <- mean(y) * nrow(yx) * apply(
      t(t(ct.m[, id]) / deno), c(1, 2),
      function(x) {
        return(min(x, trend.upper.bound))
      }
    )
  } else {
    z <- t(t(ct.m[, id]))
  }
  
  # Estimate the predicted values using different method.
  # Skip zero period when computing median.
  if (neighbour.method == "mean") {
    pred <- apply(z, 1, mean)
  } else if (neighbour.method == "median") {
    tmedian <- function(x) {
      if (sum(x) == 0) {
        return(0)
      }
      return(median(x[x != 0]))
    }
    pred <- apply(z, 1, tmedian)
  }
  
  # Scale the prediction back.
  if (scaling.option == "centered") {
    pred <- pred * sqrt(var(y)) + mean(y)
  }
  
  if (start.period > 1) {
    pred[seq_len(start.period - 1)] <- 0
  }
  return(pred)
}

#' Generate prediction for seasonal entry.
#'
#' @param y a vector of history values before treatment.
#' @param ct.m matrix with time series entries
#' @param treatperiod treatment period.
#' @param num.periods total period length that includes both history and
#'   prediction.
#' @param frequency frequency of the seasonality.
#' @param neighbour.method method for prediction once the list of neighbours
#'   are found. "median" skips zero entries when computing the median;
#'   "mean" used mean of the neighbours, "none" desables scaling.
#' @param neighbour.weighting boolean value of weighting option. If "TRUE",
#'   more recent history is weighted higher using linear weight. If
#'   "FALSE", no weight is used.
#' @param scaling.option if "centered", we use the standard scaling to treat
#'   histories: (x-mean(x))/std(x); if "mean", we use x/mean(x). If
#'   "none", no scaling is used.
#' @param neighbour.periods maximum number of periods used for prediction.
#' @param max.pretreatment.period maximum number of periods used for seasonal
#'   prediction.
#' @return A vector of prediction values.
#' @export
SeasonalPrediction <- function(y,
                               ct.m,
                               treatperiod,
                               num.periods,
                               frequency,
                               num.neighbours,
                               neighbour.method,
                               neighbour.weighting = TRUE,
                               scaling.option,
                               neighbour.periods,
                               max.pretreatment.period) {
  original.y <- y
  start.period <- min(which(y > 0))
  if (length(y) - start.period > max.pretreatment.period) {
    start.period <- length(y) - max.pretreatment.period
  }
  y <- y[seq(from = start.period, to = length(y), by = 1)]
  if (sum(y) == 0) {
    return(rep(0, num.periods))
  }
  
  # If the history is shorter than 2 frequency cycles or there are zero history
  # periods, use airline model. Otherwise we apply stl to de-seasonalize the
  # history into trend and seasonal components. We use nearest neighbour method
  # to get the trend prediction and then multiply the seasonality back.
  if (length(y) < 2 * frequency || min(y) <= 0) {
    fit <- try(forecast::Arima(ts(y),
                     order = c(0, 1, 1),
                     seasonal = list(
                       order = c(0, 1, 1),
                       period = frequency,
                       lambda = 0
                     )
    ), silent = TRUE)
    
    # The fitted values or ARIMA model is the one-step forecast.
    if (class(fit) != "try-error") {
      pred <- c(fitted(fit), forecast::forecast(fit,
                                      h = num.periods - treatperiod + 1
      )$mean)
      if (start.period > 1) {
        pred <- c(original.y[seq_len(start.period - 1)], pred)
      }
    } else {
      pred <- NearestNeighbourPrediction(y,
                                         ct.m,
                                         treatperiod,
                                         num.neighbours = num.neighbours,
                                         neighbour.method = neighbour.method,
                                         neighbour.weighting =
                                           neighbour.weighting,
                                         scaling.option = scaling.option,
                                         neighbour.periods =
                                           neighbour.periods
      )
    }
  } else {
    fit.stl <- try(stats::stl(ts(log(y), frequency = frequency), "per"), silent = TRUE)
    if (class(fit.stl) == "try-error") {
      pred <- NearestNeighbourPrediction(y,
                                         ct.m,
                                         treatperiod,
                                         num.neighbours = num.neighbours,
                                         neighbour.method = neighbour.method,
                                         neighbour.weighting =
                                           neighbour.weighting,
                                         scaling.option = scaling.option,
                                         neighbour.periods =
                                           neighbour.periods
      )
    } else {
      detrend.y <- exp(fit.stl$time.series[, 2] + fit.stl$time.series[, 3])
      season <- exp(fit.stl$time.series[, 1])
      pred <- NearestNeighbourPrediction(detrend.y,
                                         ct.m,
                                         treatperiod,
                                         num.neighbours = num.neighbours,
                                         neighbour.method = neighbour.method,
                                         neighbour.weighting =
                                           neighbour.weighting,
                                         scaling.option = scaling.option,
                                         neighbour.periods =
                                           neighbour.periods
      )
      pred <- pred * season[mod(seq_len(length(pred)) - 1, frequency) + 1]
    }
  }
  return(pred)
}





















#AK implementation of doudchenko Imbens style estimators

flexible_scm <- function(data_full, id_var = "entry", time_var = "period", treat_indicator = "treatperiod_0", outcome_var = "target",
                         counterfac_var = "counter_factual"){
  tr_entries <- data_full %>%
    dplyr::filter(!!as.name(treat_indicator) > 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull() %>%
    sort()
  ct_entries <- setdiff(data_full %>% dplyr::distinct(!!as.name(id_var)) %>% dplyr::pull(), tr_entries)
  
  
  # create control data frame, with a new id for the sake of ordering observations later
  control_data <- data_full %>% dplyr::filter(!!as.name(id_var) %in% ct_entries)


  treat_data <- data_full %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries) %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var)) %>%
    dplyr::group_by(!!as.name(id_var)) %>%
    dplyr::mutate(Treatment_Period = length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1) %>%
    dplyr::ungroup()
  
  
  # create the control matrix once, which is an input to sdid estimator
  #NxT
  control_matrix <- tidyr::spread(
    control_data %>% dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(outcome_var)),
    !!as.name(time_var), !!as.name(outcome_var)
  ) %>%
    dplyr::select(-!!as.name(id_var)) %>%
    as.matrix() %>%
    t()
  
  split_treat_data_pre <- treat_data %>% 
    dplyr::filter(!!as.name(time_var)<Treatment_Period) %>%
    dplyr::select(tidyselect::all_of(c(id_var, outcome_var))) %>%
    dplyr::group_by(!!as.name(id_var)) %>%
    dplyr::group_split( .keep = F) %>%
    lapply(., function(x) x %>% dplyr::pull())

    # for each treated unit, find when it was treated
  list_of_treat_times <- treat_data %>%
    split(.[[id_var]]) %>%
    lapply(., function(x) {
      x %>% dplyr::select(Treatment_Period) %>% dplyr::slice(1) %>% dplyr::pull()
    })
  
  list_inputed_y=furrr::future_map2(.x=split_treat_data_pre,
                                    .y=list_of_treat_times,
                                    .f=~scm_imputation(treat_data = .x,
                                                   control_mat = control_matrix,
                                                   treat_time = .y) )
  
  flex_scm_series=Map(dplyr::bind_cols, treat_data %>% split(.[[id_var]]),
                      lapply(list_inputed_y, function(x){
                        return(tibble::tibble("point.pred"=x))
                      } )) %>% 
    dplyr::bind_rows() %>% 
    dplyr::rename(response=!!as.name(outcome_var)) %>%
    dplyr::mutate(point.effect = response - point.pred)

    if (!is.null(counterfac_var)) {
    flex_scm_series <- flex_scm_series %>%
      dplyr::mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = (response / !!as.name(counterfac_var)) - 1
      )
  }
  
  # add a column with relative (pct) effect
  flex_scm_series <- flex_scm_series %>% 
    dplyr::mutate(pct.effect = (response / point.pred) - 1 ) %>%
    dplyr::select(tidyselect::all_of(c(time_var, id_var,"point.pred" ,"response",
                                       "Treatment_Period", "point.effect", 
                                       counterfac_var,"cf_point.effect",
                                       "cf_pct.effect", "pct.effect"))) %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var))
  
}

scm_imputation <- function(treat_data, control_mat, treat_time){
  #define the control data matrix
  control_pre=control_mat[seq_len(treat_time-1),]
  #fit a speedy Cross-Model Selection and Averaging Grid for ENP
  fit_enp = bigstatsr::big_spLinReg(bigstatsr::as_FBM(control_pre), 
                                    treat_data, alphas = c(1e-4,0.2,0.5,0.8,1), 
                                    warn = F)
  
  #summary(fit_enp, best.only=T)
  imputed_y=predict(fit_enp, bigstatsr::as_FBM(control_mat))
  return(imputed_y)
}



flexible_sdid <- function(data_full, id_var = "entry", time_var = "period", treat_indicator = "treatperiod_0", outcome_var = "target",
                         counterfac_var = "counter_factual"){
  tr_entries <- data_full %>%
    dplyr::filter(!!as.name(treat_indicator) > 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull() %>%
    sort()
  ct_entries <- setdiff(data_full %>% dplyr::distinct(!!as.name(id_var)) %>% dplyr::pull(), tr_entries)
  
  
  # create control data frame, with a new id for the sake of ordering observations later
  control_data <- data_full %>% dplyr::filter(!!as.name(id_var) %in% ct_entries)
  
  
  treat_data <- data_full %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries) %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var)) %>%
    dplyr::group_by(!!as.name(id_var)) %>%
    dplyr::mutate(Treatment_Period = length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1) %>%
    dplyr::ungroup()
  
  
  # create the control matrix once, which is an input to sdid estimator
  #NxT
  control_matrix <- tidyr::spread(
    control_data %>% dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(outcome_var)),
    !!as.name(time_var), !!as.name(outcome_var)
  ) %>%
    dplyr::select(-!!as.name(id_var)) %>%
    as.matrix() %>%
    t()
  
  split_treat_data_pre <- treat_data %>% 
    dplyr::filter(!!as.name(time_var)<Treatment_Period) %>%
    dplyr::select(tidyselect::all_of(c(id_var, outcome_var))) %>%
    dplyr::group_by(!!as.name(id_var)) %>%
    dplyr::group_split( .keep = F) %>%
    lapply(., function(x) x %>% dplyr::pull())
  
  # for each treated unit, find when it was treated
  list_of_treat_times <- treat_data %>%
    split(.[[id_var]]) %>%
    lapply(., function(x) {
      x %>% dplyr::select(Treatment_Period) %>% dplyr::slice(1) %>% dplyr::pull()
    })
  
  browser()
  #gonna need weights from this method
  #Then, need time weights (similar regression except now the outcome var
  # is the first Treatment Period for each control unit, and the we predict using
  #all the past values for that particular unit?)
  list_inputed_y=furrr::future_map2(.x=split_treat_data_pre,
                                    .y=list_of_treat_times,
                                    .f=~scm_weights(treat_data = .x,
                                                    control_mat = control_matrix,
                                                    treat_time = .y) )
  
  
  flex_scm_series=Map(dplyr::bind_cols, treat_data %>% split(.[[id_var]]),
                      lapply(list_inputed_y, function(x){
                        return(tibble::tibble("point.pred"=x))
                      } )) %>% 
    dplyr::bind_rows() %>% 
    dplyr::rename(response=!!as.name(outcome_var)) %>%
    dplyr::mutate(point.effect = response - point.pred)
  
  if (!is.null(counterfac_var)) {
    flex_scm_series <- flex_scm_series %>%
      dplyr::mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = (response / !!as.name(counterfac_var)) - 1
      )
  }
  
  # add a column with relative (pct) effect
  flex_scm_series <- flex_scm_series %>% 
    dplyr::mutate(pct.effect = (response / point.pred) - 1 ) %>%
    dplyr::select(tidyselect::all_of(c(time_var, id_var,"point.pred" ,"response",
                                       "Treatment_Period", "point.effect", 
                                       counterfac_var,"cf_point.effect",
                                       "cf_pct.effect", "pct.effect"))) %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var))
  
}




scm_weights <- function(treat_data, control_mat, treat_time){
  #define the control data matrix
  control_pre=control_mat[seq_len(treat_time-1),]
  #fit a speedy Cross-Model Selection and Averaging Grid for ENP
  fit_enp = bigstatsr::big_spLinReg(bigstatsr::as_FBM(control_pre), 
                                    treat_data, alphas = c(1e-4,0.2,0.5,0.8,1), 
                                    warn = F)
  
  #summary(fit_enp, best.only=T)
  imputed_y=predict(fit_enp, bigstatsr::as_FBM(control_mat))
  return(imputed_y)
}








############################################
# Functions to create a placebo dataset,
# assigning treatment based on matching.
###########################################
nearest_ts_euclidean <- function(ts_tomatch, ts_rest) {
  # Helper function for matching_without_replacement
  # Measures the euclidean distance between ts_tomatch
  # and each of the rows in the tibble ts_rest.
  
  # Args
  # ts_tomatch: treated time series tibble, rows identify ID and columns Time
  # ts_rest: donor time series tibble, rows identify ID and columns Time.
  
  # Output
  # Vector of indices indicating the row in ts_rest that are best match
  # (min L2 norm distance) to the row in ts_tomatch.
  # (Vector length equal to number of rows in ts_tomatch, ie treated entries)
  # If ts_tomatch has more than 1 entry, the matching is done with replacement
  # (multiple treated can have the same match).
  
  return(apply(ts_tomatch, 1, function(ts_tomatch) {
    which.min(
      apply(
        ts_rest, 1, function(ts_rest, ts_tomatch) {
          stats::dist(rbind(ts_rest, ts_tomatch))
        },
        ts_tomatch
      )
    )
  }))
}


matching_without_replacement <- function(treated_block,
                                         control_block,
                                         id_var,
                                         treat_period) {
  # finds the nearest match for each treated unit (treated_block)
  # among the donor pool (control_block) without replacement by calling helper
  # to nearest_ts_euclidean.
  
  # Args
  # treated_block: treated time series tibble, rows are ID and columns Time
  # control_block: donor time series tibble, rows are ID and columns Time.
  
  # Output
  # df_toreturn, a tibble containing a column for the placebo-treated unit 
  # ID numbers, the treated unit it was the nearest match to, 
  # and the time that treated unit was actually
  # treated (num_rows of the dataframe equal to num_rows of treated_block).
  
  # initialize an empty vector for the donor IDs that match
  already_matched <- c()
  # Store the time of treatment and treated ID for the true treated units
  placebo_treat_period <- treated_block %>% dplyr::pull(tidyselect::all_of(treat_period))
  treatment_unit <- treated_block %>% dplyr::pull(!!as.name(id_var))
  for (i in seq_len(nrow(treated_block))) {
    #Find nearest euclidean match among unmatched controls for 
    #each treated observation
    temp_match <- nearest_ts_euclidean(
      treated_block %>%
        dplyr::slice(i) %>%
        dplyr::select(-tidyselect::all_of(c(id_var,treat_period))),
      control_block %>%
        dplyr::filter(!!as.name(id_var) %in% setdiff(!!as.name(id_var), 
                                                     already_matched)) %>% 
        dplyr::select(-tidyselect::all_of(id_var)))
    #Update the vector of already matched donors
    already_matched[i] <- control_block %>%
      dplyr::filter(!!as.name(id_var) %in% 
                      setdiff(!!as.name(id_var), already_matched)) %>%
      dplyr::slice(temp_match) %>%
      dplyr::pull(!!as.name(id_var))
  }
  # Store the resulting vectors in a tibble for output
  df_toreturn <- tibble::tibble(temp_id = already_matched, 
                                temp_treattime = placebo_treat_period, 
                                Treatment_Unit = treatment_unit)
  
  return(df_toreturn %>% 
           dplyr::rename(!!as.name(id_var) := temp_id,
                         !!as.name(treat_period) := temp_treattime))
}



#attempting to match based on TS features and compare the performance (seems worse)
#TODO(alexdkellogg): account for multiple treatments, treatment end dates
create_placebo_df <- function(data_full, id_var = "entry",
                              time_var = "period", 
                              treat_indicator = "treatperiod_0",
                              outcome_var = "target", 
                              counterfac_var = "counter_factual",
                              match_type="feature") {
  # Generates a placebo tibble, using matching methods to select
  # placebo-treated entries as those most similar to truly-treated.
  
  # Args
  # data_full: long-form dataframe with both treated and control entries.
  #            rows of the df represent period-entry combinations 
  #            (eg N (total num of entry) rows for period t).
  #            each row should have a treatment indicator (treat_indicator), 
  #            a period number (time_var),
  #            an individual ID (id_var), and an outcome (outcome_var)
  #            for associated with that period-ID combination
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, 
  #          in increasing order (eg 0 is the first time, 120 is the last)
  # treat_indicator:column name of binary (0, 1) indicator for whether id_var i 
  #                 in time_var t was treated. 
  #                 Once treated, must always be treated.
  
  
  # Output
  # placebo_df_long, a tibble of the same format as data_full, entirely
  # consisting of donor units, some of which are placebo-treated 
  # (based on matching).
  
  # Split the dataset based on whether the unit is ever treated
  tr_entries <- data_full %>%
    dplyr::filter(!!as.name(treat_indicator) > 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull()
  
  ct_entries <- setdiff(data_full %>%
                          dplyr::distinct(!!as.name(id_var)) %>%
                          dplyr::pull(), tr_entries)
  # Create a dataframe pf the subset of control units
  cd <- data_full %>% dplyr::filter(!!as.name(id_var) %in% ct_entries)
  cd <- merge(
    cd,
    data.frame(
      entry = ct_entries,
      rank = seq_along(ct_entries)
    )
  )
  
  # Pivot the long cd dataframe to wide. Each row will represent an id_var, 
  #with columns for the outcome_var at each time period
  cd_for_match <- tidyr::pivot_wider(
    data = cd %>%
      dplyr::arrange(!!as.name(time_var), !!as.name(id_var)),
    names_from = !!as.name(time_var),
    id_cols = c(!!as.name(id_var)),
    values_from = c(!!as.name(outcome_var))
  )
  
  # Store treated data in a tibble
  treated_to_match <- data_full %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries)
  
  #Create a variable indicating the time treatment is assigned by unit
  treated_to_match <- treated_to_match %>%
    dplyr::group_by(!!as.name(id_var)) %>%
    dplyr::mutate(
      Treatment_Period =
        length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1
    ) %>%
    dplyr::ungroup()
  
  # for pivoting, potential issues arise if we have several time varying 
  # covariates -- added as TODO.
  # we'd have to take the values_from each of them, 
  #and for any constant args we'd presumably have to add them to id_cols
  #Pivot the data wide, so each row has all the time series data for 
  #a given unit.
  data_wide_m <- tidyr::pivot_wider(
    data = treated_to_match,
    names_from = !!as.name(time_var),
    id_cols = c(!!as.name(id_var), Treatment_Period),
    values_from = c(!!as.name(outcome_var))
  ) # %>% as.matrix()
  
  
  list_of_ts <- data_full %>%
    dplyr::select(!!as.name(id_var), !!as.name(outcome_var)) %>%
    split(.[[id_var]]) %>%
    furrr::future_map(~ .[[outcome_var]]) %>%
    furrr::future_map(~ stats::ts(.))
  df_feat <- tsfeatures::tsfeatures(list_of_ts) %>% 
    dplyr::mutate(!!as.name(id_var):=seq_len(dplyr::n())) 
  
  treated_ts_features= df_feat %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries) %>%
    dplyr::left_join(
      treated_to_match %>%
        dplyr::select(tidyselect::all_of(c(id_var, "Treatment_Period"))) %>%
        dplyr::distinct(), by=c(id_var) )
  
  control_ts_features= df_feat %>%
    dplyr::filter(!!as.name(id_var) %in% ct_entries)
  
  if(match_type=="ts"){
    # Determine control units to assign to placebo treatment group via matching.
    matched_placebo_df_temp <- matching_without_replacement(data_wide_m, 
                                                            cd_for_match, 
                                                            id_var,
                                                            "Treatment_Period")
  }
  
  if(match_type=="feature"){
    matched_placebo_df_temp<- matching_without_replacement(treated_ts_features, 
                                                           control_ts_features, 
                                                           id_var,
                                                           "Treatment_Period")
  }

  return(matched_placebo_df_temp)
}


pair_distance_helper=function(pair_map, data_full){
  control_series=pair_map %>% 
    dplyr::select(entry) %>%
    dplyr::inner_join(data_full, by="entry") %>%
    dplyr::select(entry, period, target)
  
  treatment_series=pair_map %>% 
    dplyr::select(Treatment_Unit) %>%
    dplyr::rename(entry="Treatment_Unit") %>%
    dplyr::inner_join(data_full, by="entry") %>%
    dplyr::select(entry, period, target) %>%
    dplyr::rename(treat_unit="entry",
                  target_treat="target")
  
  matched_series=control_series %>%
    dplyr::left_join(treatment_series, by="period") %>%
    dplyr::mutate(dif_sq=(target-target_treat)^2)
  
  distance_result=matched_series %>%
    dplyr::summarise(sum_sq_dif=sum(dif_sq)) %>%
    dplyr::pull()
  
  return(distance_result)
  
}
evaluate_matches <- function(data_full, match_map){
  #PLAN: grab each match pair, turn it into a single dataframe
  #of series 1 and series2
  #compute the L2 distance
  #apply over each row in the mapping
  
  split_pairs=match_map %>% split(.[["entry"]])
  
  dist_list=furrr::future_map(.x=split_pairs,
                             .f=~pair_distance_helper(data_full=data_full,
                                                  pair_map=.x))
  
  return(dist_list %>% unlist() %>% sum() )
}





#Messing around with IV's and TSCS data
set.seed(1982)
N <- 1000 # Number of observations
T=20
Treat_Time=T/2
PX <- 5 # Number of exogenous variables
errors <- MASS::mvrnorm(N*T, rep(0, 2), rockchalk::lazyCor(X = 0.5, d = 2))
X_exog <- MASS::mvrnorm(N, rep(0, PX), cor(matrix(rnorm(PX*(PX+5)),PX+5,PX)))

# Instrument
Z <- sample(0:1, N, T)

time_grid <-tidyr::expand_grid(time = 1:T,
                               unit = 1:N) %>%
  dplyr::bind_cols(errors1=errors[,1],errors2=errors[,2]) %>% 
  dplyr::inner_join(
    tibble::tibble(time=1:T,time_re=rnorm(T)), by="time") 
    
    

# First stage
unobserved_adoption_driver=rnorm(N, 0,2)
true_gammas <- rnorm(5, 0, 0.1)
u <- -0.5 + X_exog %*% true_gammas + 1.5* Z - unobserved_adoption_driver + rnorm(N,0,0.5)
p <- 1/(1+exp(-u))
u <- runif(N)
X_endog <- dplyr::case_when(p > u ~ 1,TRUE ~ 0)

X_tib=tibble(unit=1:N,
             X1 = X_exog[,1],
             X2 = X_exog[,2],
             X3 = X_exog[,3],
             X4 = X_exog[,4],
             X5 = X_exog[,5],
             Z=Z,
             X_endog=X_endog,
             Treat_Time=Treat_Time,
             Unobs=unobserved_adoption_driver) 

fake_data= X_tib %>% 
  dplyr::inner_join(time_grid, by="unit") %>%
  dplyr::mutate(X_endog=(time>=Treat_Time)*X_endog,
                Post=(time>=Treat_Time)) %>%
  dplyr::group_by(unit) %>%
  dplyr::mutate(Treated=max(X_endog)) %>%
  dplyr::ungroup()

X_tib %>% group_by(Z) %>% summarise(Treat_Prob_z0=mean(X_endog))
X_tib %>% group_by(X_endog) %>% summarise(unobs_mean=mean(Unobs))


# Second Stage
true_deltas <- rnorm(5, 0, 1)

fake_data=fake_data %>% 
  dplyr::mutate(
  int=0.2+ (as.matrix(dplyr::select(fake_data, X1:X5)) %*% true_deltas)[,1],
  Y_outcome=int+0.89*X_endog+1.5*Unobs+rnorm(dplyr::n(),0,0.5)
)

  
#saveRDS(synthetic_data, 'synthetic_iv_truth_089.RDS')
sd_y <- sd(fake_data$Y_outcome)


# naive_reg=lm(Y_outcome~X_endog+X1+X2+X3+X4+X5, data=synthetic_data)
# iv_reg= AER::ivreg(Y_outcome~X_endog+X1+X2+X3+X4+X5|Z+X1+X2+X3+X4+X5, data=synthetic_data)

did_reg=lm(Y_outcome~X_endog+time+(X_endog*time), data=fake_data)

itt <- lm(X_endog ~ Z + X1+X2+X3+X4+X5, data=fake_data %>% 
            dplyr::filter(time==Treat_Time))

num=lm(Y_outcome ~ Z + X1+X2+X3+X4+X5, data=fake_data %>% 
               dplyr::filter(time==Treat_Time))


#Causal Impact gets it right b/c selection is not time-varying
aa=estimate_causalimpact_series(fake_data,id_var = "unit", time_var = "time",
                                         treat_indicator = "X_endog",
                                         outcome_var = "Y_outcome",
                                         counterfac_var = NULL)

aa_tot=compute_tot_se_jackknife(aa, time_var = "time", treat_period_var = "Treatment_Period",
                                     pred_var = "point.pred", outcome_var = "response",
                                     stat_in = "mean", alpha_ci = 0.95,
                                     compute_cf_eff = F, counterfac_var = NULL,
                                     post_treat_only=T)


set.seed(1982)
N <- 10000 # Number of observations
PX <- 5 # Number of exogenous variables
errors <- MASS::mvrnorm(N, rep(0, 2), rockchalk::lazyCor(X = 0.5, d = 2))
X_exog <- MASS::mvrnorm(N, rep(0, PX), cor(matrix(rnorm(PX*(PX+5)), PX+5, PX)))

Z <- sample(0:1, N, T) # Instrument
# First stage
true_gammas <- rnorm(5, 0, 0.1)
u <- -0.5 + X_exog %*% true_gammas + 0.9 * Z + errors[,1]
p <- 1/(1+exp(-u))
u <- runif(N)
X_endog <- dplyr::case_when(p > u ~ 1,
                            
                            TRUE ~ 0)
# Second Stage
true_deltas <- rnorm(5, 0, 1)
Y_outcome <- 0.2 + X_exog %*% true_deltas + 0.89 * X_endog + errors[,2]
  # tibble
synthetic_data <- tibble::tibble(Y_outcome = Y_outcome,
                                 
                                 X_endog = X_endog,
                                 Z = Z,
                                 X1 = X_exog[,1],
                                 X2 = X_exog[,2],
                                 X3 = X_exog[,3],
                                 X4 = X_exog[,4],
                                 X5 = X_exog[,5])

#saveRDS(synthetic_data, 'synthetic_iv_truth_089.RDS')
sd_y <- sd(synthetic_data$Y_outcome)

naive_reg=lm(Y_outcome~X_endog+X1+X2+X3+X4+X5, data=synthetic_data)
iv_reg= AER::ivreg(Y_outcome~X_endog+X1+X2+X3+X4+X5|Z+X1+X2+X3+X4+X5, data=synthetic_data)

itt <- lm(X_endog ~ Z + X1+X2+X3+X4+X5, data=synthetic_data)

num=lm(Y_outcome ~ Z + X1+X2+X3+X4+X5, data=synthetic_data)



#explore why data is lacking overlap where it shouldn't
load(here::here("Data", "Variations","aa_noisy_factors.RData" ))
random_data=sample(1:n_seeds,1)
tsfeature_by_treat_df(formatted_data[[random_data]])
tsfeature_pc_by_treatment_plot(formatted_data[[random_data]])

data_requested %>%
  group_by(treated) %>%
  summarise(mean_rho=mean(autocorr),
            mean_int=mean(intercept),
            mean_load1=mean(loading1),
            mean_load2=mean(loading2),
            mean_load3=mean(loading3),
            mean_load4=mean(loading4),
            mean_load5=mean(loading5),
            mean_load6=mean(loading6))
  



AA_data_re=gen_data(N = 800,  N_time = 60, treat_start=15,
                    tau_one_zero = 0, tau_decay=0.9, gamma = 0.0009,
                    seed = 1982, prop_treated = 0.5,  rho_y = 0.4, 
                    type = c("random"))



n_seeds <- 5
seeds <- sample(1000:9999, size = n_seeds)

AA_data_list <- future_map(.x = seeds, .f = ~gen_data(N = 800, prop_treated = 0.2, tau_one_zero = 0.0, gamma = 0.0, rho_y = 0.7, N_time = 100, type = "random", seed = .x))


library(here)
#still an issue about where you open from -- if the file itself in r_code, then don't need the r_code part.
source(here("/r_code/panel_estimation.R"))
source(here("r_code/treatment_effect_bootstrap.R"))
source(here("r_code/analysis_metrics.R"))
source(here("r_code/analysis_visualizations.R"))

set.seed(42)

library(tictoc) 
tic()
gsynth_ife_AA=future_pmap(list(AA_data_list, se_est=F), estimate_gsynth_series)
toc()
#get bootstrapped CIs, counterfactual effects
gsynth_ife_bootstrapped_pct_att_AA=future_map(gsynth_ife_AA,compute_tot_se_jackknife,counterfac_var="counter_factual", stat_in="mean" )
#store gap plots
gsynth_ife_bootstrapped_pct_att_AA_plots=future_map(gsynth_ife_bootstrapped_pct_att_AA,create_gap_ci_plot,
                                                    plot_title="Gsynth Mean ToT by Post Treat Period, Jackknife CI", plot_x_lab="Post Treat Period", plot_y_lab="Mean ToT",effect_var="jackknife_median_pct_tot", cf_var="median_pct_cf_tot")
#compute MSE and MAE bu post period
gsynth_ife_pct_metric_byT_AA=future_pmap(list(gsynth_ife_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric_per_t)
#compute overall MSE, MAE
gsynth_ife_pct_metric_AA=future_pmap(list(gsynth_ife_AA,  metric_str="both", pct_eff_flag=T),compute_avg_metric)
#compute "coverage" -- not sure if that's exactly it, but the frac of time the true TE is in the bounds
gsynth_ife_pct_coverage_AA=future_map(gsynth_ife_bootstrapped_pct_att_AA, compute_tot_coverage)



# #change counter_factual to y0 and target to y1
# #
# tic("Starting DGP")
# AA_data_no_sel=future_map(.x=seeds, .f=~gen_data_lubr(N=100, date_start="2018-01-01",
#                          first_treat="2019-01-01", date_end="2020-07-04", 
#                          tau_one_zero=0, T_freq = "weekly",
#                          gamma=0, prop_treated=0.25, rho_y=0.3, type="random", seed=.x))
# toc()

tic("Estimating Gsynth")
gsynth_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se_est=F, outcome_var="y", counterfac_var="y0")
gsynth_tot_AA=future_map(gsynth_AA, compute_tot_se_jackknife, counterfac_var="counter_factual")
toc()








#propensity score

aa_prop_score=gen_data_lubr(N=100, date_start="2018-01-01",
                  first_treat ="2019-01-01", date_end="2020-07-04", 
                  tau_one_zero=0, T_freq = "monthly",
                  gamma=0, prop_treated=0.25, rho_y=0.3, type="observables", seed=19)

treat_reg_data=aa_prop_score %>% select(contains("v"), treated, entry) %>% distinct(entry, .keep_all = T) %>% select(-entry)
temp_prop_reg=lm(treated~ ., data=treat_reg_data)
summary(temp_prop_reg)

aa_prop_score_unobs=gen_data_lubr(N=100, date_start="2018-01-01",
                            first_treat ="2019-01-01", date_end="2020-07-04", 
                            tau_one_zero=0, T_freq = "monthly",
                            gamma=0, prop_treated=0.25, rho_y=0.3, type="unobservables", seed=19)

treat_reg_data=aa_prop_score_unobs %>% select(contains("v"), treated, entry) %>% distinct(entry, .keep_all = T) %>% select(-entry)
temp_prop_reg=lm(treated~ ., data=treat_reg_data)
summary(temp_prop_reg)






#Attempt at Daily data
gen_dataaktest <- function(N = 800, N_months = 60, T_freq = c("monthly", "weekly", "daily"), treat_start = 15,
                           tau_one_zero = 0.02, tau_decay = 0.9, gamma = 0.0009,
                           seed = 1982, prop_treated = 0.5, rho_y = 0.4,
                           type = c("random", "observables")) {
  
  # Generates tibble of long form panel data along with a number of x variables, the counterfactual outcome, and the target.
  # Each row is time period x unit unique combnation.
  
  # Args
  # N: number of units to be generated
  # N_time: number of total time periods
  # treat_start: first treatment period. For now, N_time-treat_start>5
  # tau_one_zero: initial period treatment impact
  # tau_decay: fraction multiplied by tau_one_zero for each post treatment period after the impact (for all time)
  # thus, if set to 0.9, the impact is originally 1, then 0.9, then 0.81 in the subsequent periods
  # gamma: multiplier on the impact based on the total number of pitches (treatment intensity) -- impact+gamma*log(total_pitches)
  # seed: random number seed
  # prop_treated: proportion of the N entries that should receive treatment, roughly
  
  # rho_y: scales up the sd of the first period TE (normal sd_0=sd_e/(sqrt(1 - rho_y ^ 2)))
  # also determines how much the center of the normal dist for y_it moves
  # given prior errors --  mu_it <- rho_y * (y_t0[i, (t - 1)] - u[i, (t - 1)]) +u[i, t]
  # type: one of "random" or "observables" (to fix -- unobservables?) determining assignment to treatment
  
  
  # Higher level description of data_full:
  # rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  # each row should have a treatment indicator (treat_indicator), a period number (time_var),
  # an individual ID (id_var), and an outcome (outcome_var)
  # for associated with that period-ID combination
  
  
  # Output
  # placebo_df_long, a dataframe of the same format as data_full, but now entirely
  # consists of donor units, some of which are placebo-treated (based on matching).
  set.seed(seed)
  if (missing(type)) {
    type <- "random"
  } else {
    type <- match.arg(type)
  }
  
  if (missing(T_freq)) {
    T_freq <- "monthly"
  } else {
    T_freq <- match.arg(T_freq)
  }
  
  stopifnot(prop_treated > 0 & prop_treated <= 1)
  
  # heuristics of 4 weeks per month, and 30 days per month
  T_freq_mapping <- tibble(freqs = c("monthly", "weekly", "daily"), multiplier = c(1, 4, 30))
  T_adjustment <- T_freq_mapping$multiplier[which(T_freq_mapping$freqs == T_freq)]
  N_time <- N_months * T_adjustment
  
  
  
  stopifnot(N_months - treat_start > 5)
  # Generate time-invariant covariates
  
  # these are x variables drawn from a number of different distributions
  # there are a total of 24 variables, so that the Xs are Nx24
  Sigma <- matrix(c(10, 3, 3, 2), 2, 2)
  vs <- MASS::mvrnorm(n = N, rep(0, 2), Sigma)
  v3 <- rbeta(
    n = N,
    shape1 = 2,
    shape2 = 2
  )
  v4 <- rbinom(n = N, size = 10, prob = 0.2)
  v5 <- matrix(rnorm(N * 20), N, 20)
  # colnames(v5) <- glue::glue("v{5:24}") #redundant
  xs <- cbind(vs, v3, v4, v5)
  colnames(xs) <- glue::glue("v{1:24}")
  xs <- as_tibble(xs) %>%
    mutate(division_country = 1:N)
  
  
  # Vertical
  # creating a tibble with 12 rows (one for each vertical),
  # with columns for the name of the vertical,
  # it's probability (ie prob that an observations falls into that vertical)
  # and the random effect shock associated with that variable
  
  # Q for Ignacio: why this rather than creating a tibble outright?
  vertical_re <-
    structure(
      list(
        # name of verticals
        vertical = c(
          "A", "B", "C", "D", "E", "F", "G",
          "H", "I", "J", "K", "L"
        ),
        # probability of belonging in each vertical
        p = c(
          0.0629,
          0.039,
          0.0343,
          0.1943,
          0.0848,
          0.0895,
          0.12,
          0.0343,
          0.2352,
          0.001,
          0.0619,
          0.0428
        ),
        # random effect associated with each verical
        vertical_re = c(
          `vertical[1]` = 0.2212,
          `vertical[2]` = -0.3044,
          `vertical[3]` = 0.3183,
          `vertical[4]` = -0.6042,
          `vertical[5]` = -0.0863,
          `vertical[6]` = 0.5209,
          `vertical[7]` = -0.3082,
          `vertical[8]` = -0.0151,
          `vertical[9]` = 0.0912,
          `vertical[10]` = 0.0544,
          `vertical[11]` = 0.2699,
          `vertical[12]` = -0.1517
        )
      ),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -12L)
    )
  
  # crosswalk that creates a tibble with a row for each unit (N)
  # and the vertical they have been randomly assigned to
  # assignment is done by sampling within creation, with replacement of course.
  # and according to the pmf defined in vertical_re$p
  verticalxwalk <- tibble(
    division_country = 1:N,
    vertical = sample(
      x = LETTERS[1:12],
      size = N,
      replace = TRUE,
      prob = vertical_re$p
    )
  )
  
  # each unit assigned a vertical, we map that back to their row in our data (xs) using join
  # Then, given the vertical they are assigned, we map on their vertical random effect
  # (ASK Ignacio -- dropped immediately, not used for anything? -- sol: merged later)
  # then, we create a new variable, v25, that is 0 unless they are in vertical E or K
  xs <- inner_join(xs, verticalxwalk, by = "division_country") %>%
    inner_join(vertical_re, by = "vertical") %>%
    mutate(v25 = case_when(
      vertical == "E" ~ -3.1,
      vertical == "K" ~ -2.02,
      TRUE ~ 0
    )) %>%
    select(-vertical, -p, -vertical_re)
  
  
  
  # Treatment
  # Type determines whether the selection into treatement is random or not
  if (type == "random") {
    # if randome assignment of treatment
    # create a tibble with N rows (one per unit)
    # and sample treatment (0 if never, 1 if ever) using argument prop_treated
    treated <- tibble(
      division_country = 1:N,
      treated = sample(
        x = c(0, 1),
        size = N,
        prob = c(1 - prop_treated, prop_treated),
        replace = TRUE
      )
    )
  } else { # if treatment is based on observables,
    # select a subset of the variables (v1,3,5,6,8,9,11,12,13,14,16,18,22,23) -- why these (ask ignacio)?
    # and store them in a matrix
    z <- xs %>%
      select(
        -division_country,
        -v2,
        -v4, -v7,
        -v10,
        -v15,
        -v17, -v20,
        -v21,
        -v24
      ) %>%
      as.matrix()
    # draw a vec of standard normal rvs, with length of the num of variables in z
    beta <- rnorm(
      n = ncol(z),
      mean = 0,
      sd = 1
    )
    # compute the "z_hat" by multiplying the variables for each invidual by their beta
    # this will be N (num units) by 1 dimensional
    zbeta <- z %*% beta
    
    colnames(zbeta) <- "zbeta"
    
    # zbeta can range from negative to positive, so we rescale it using the
    # logistic function to map it into the probability space (0,1)
    p <- zbeta %>%
      as_tibble() %>%
      mutate(p = 1 / (1 + exp(zbeta)))
    
    # Ask ignacio: could we just take the top "prop_treated" from the logistic transformation above?
    # or, to have more randomness, subtract a random uniform from each and take the top frac?
    
    # while the average probability is different than our ultimate goal from prop_treated
    # adjust each individual p by the difference.
    # if the average p is below our goal, we will wind up adding to p. if above, we will subtract
    while (abs(mean(p$p) - prop_treated) > 0.0001) {
      add_p <- prop_treated - mean(p$p)
      p <- p %>%
        rowwise() %>%
        mutate(p = min(1, p + add_p))
    }
    # to actually assign treatment, we take our adjusted p from above
    # we compare it to a random uniform (for each indiv)
    # if it's larger, we say they are treated, else not treated
    # store the division (just the row, or unit) and the treated status
    treated <- p %>%
      ungroup() %>%
      mutate(
        u = runif(n = N, min = 0, max = 1),
        treated = case_when(
          p > u ~ 1,
          p <= u ~ 0
        ),
        division_country = 1:N
      ) %>%
      select(division_country, treated)
  }
  
  # count up the num control and treated
  N_treated <- sum(treated$treated)
  N_control <- N - N_treated
  
  # First period is a tibble that assigns, for each unit
  # the first period of treatment. output is an Nx2 tibble.
  # Samples periods between 5 and Total_T-5 for each unit, and assigned
  # 9999L if they are in fact never treated
  # ask ignacio: can we parameterize the first treat (5)?
  # Also thinking of using N-1, but the frac of treated for only 1 period would be too high
  # thoughts on sampling from a dist like Alex P?
  first_period <-
    tibble(
      division_country = 1:N,
      first = sample(
        x = (treat_start * T_adjustment):(N_months * T_adjustment - 5 * T_adjustment),
        size = N,
        replace = TRUE
      )
    ) %>%
    inner_join(treated, by = "division_country") %>%
    mutate(first = case_when(
      treated == 1 ~ first,
      TRUE ~ 9999L
    )) %>%
    select(-treated)
  
  
  
  
  # a_i
  # individual random effect for each unit, drawn from normal
  # output is tibble with ID and the indiv specific effect
  a <- tibble(
    division_country = 1:N,
    a_i = rnorm(
      n = N, mean = 0.1748,
      sd = 1.4678
    )
  )
  
  # create a tibble of unit, unit random effect,
  # unit vertical (with prob of assignment and the vertical Random Effect)
  # and unit ever treated status
  fake_data <-inner_join(a, verticalxwalk, by = "division_country") %>%
    inner_join(vertical_re, by = "vertical") %>%
    inner_join(treated, by = "division_country")
  
  # d_it
  # creates a NxT by 2 tibble; each row is a unique entry (unit) x time combo
  time_grid <-
    tidyr::expand_grid(
      time = 1:N_time,
      division_country = fake_data$division_country
    )
  
  mu_d <- 0.01
  
  ## countries
  # like the verticals, create a tibble of countries (9 total),
  # with their associated probability of assignment for each unit,
  # and the random effect associated with that country
  countries <- structure(
    list(
      country = 1:9,
      country_re = c(
        -0.9231,
        0.3533,
        0.6499,
        -0.0572,
        -1.1556,
        -0.8078,
        0.4037,
        0.0955,
        1.434
      ),
      p = c(
        0.0457,
        0.1019,
        0.1038,
        0.0067,
        0.0124,
        0.0152,
        0.0305,
        0.0048,
        0.679
      )
    ),
    class = c("tbl_df", "tbl", "data.frame"),
    row.names = c(NA, -9L)
  )
  
  N_countries <- length(countries$country)
  
  
  # map each row (Id) to a country by sampling from the list of countries
  # following the specified probability of each
  # output is a tibble with the Id and country assigned to said ID
  countryxwalk <- tibble(
    division_country = 1:N,
    country = sample(
      x = countries$country,
      size = N,
      replace = TRUE,
      prob = countries$p
    )
  )
  
  # recreate fake_data tibble
  # adding columns for country and country random effect (And removing assignment prob)
  # merge with time grid as well, to make it a long panel data (row is unit x time)
  # so far, all the vars are time constant and are thus repeated
  fake_data <-
    inner_join(fake_data, countryxwalk, by = "division_country") %>%
    select(-p) %>%
    inner_join(countries, by = "country") %>%
    select(-p) %>%
    inner_join(time_grid, by = "division_country")
  
  
  ## Time
  # create a time random effect table -- each time period
  # is assigned a random shock from the normal dist below
  time_re <-
    tibble(
      time = 1:N_time,
      time_re = rnorm(n = N_time, mean = (0), sd = (0.141) * (1 + log(T_adjustment) / 30))
    )
  # merge the random effect over time.
  # each individual in a given time period will face the same shock.
  fake_data <- inner_join(fake_data, time_re, by = "time")
  
  
  
  # time country interaction
  # redundant
  N_countries <- length(countries$country)
  
  # create a tibble that, in each row, has a unique combo of time and country
  time_country_grid <-
    tidyr::expand_grid(country = 1:N_countries, time = 1:N_time)
  # tibble for each country, in time period 1, with a country specific normal shock
  r_one <-
    tibble(
      country = 1:N_countries,
      time = 1,
      r = rnorm(n = N_countries, mean = 0, sd = 1 / sqrt(T_adjustment))
    )
  
  # merge in the period 1 shock (country specific) into the country x time tibble
  r <-
    full_join(time_country_grid, r_one, by = c("country", "time")) %>%
    group_by(country)
  # now, fill out the rest of the country x time shocks
  # for each time period and country combination,
  # take the country specific shock from the time period before and add normal noise
  # For a given country, we thus have shock2=shock1+norm_shock
  tic("assign country x time shocks")
  for (t in 2:N_time) {
    for (c in 1:N_countries) {
      r$r[r$time == t &
            r$country == c] <- r$r[r$time == (t - 1) &
                                     r$country == c] + rnorm(n = 1, mean = (0), sd = (0.1))
    }
  }
  toc()
  # after filling out the grid for each country/time pair, merge into fake data
  fake_data <- inner_join(fake_data, r, by = c("country", "time"))
  
  ## pitches
  # add to our fake data tibble two columns: the time of first treatment (unit specific)
  # and the relation of the current row to treat time (-1 would be period before treat)
  fake_data <-
    inner_join(fake_data, first_period, by = "division_country") %>%
    mutate(m = time - first)
  # find the longest post treat period of pitch 1 for any unit in data
  tau_one <- rep(NA, max(fake_data$m + 1))
  # adjust the decay so that it decays by the proper amount each month
  # need tau_decay_adj**T_adjustment=tau_decay
  tau_decay_adj <- tau_decay**(1 / T_adjustment)
  # define a treatment decay function: for period one, tau_one is an argument
  # then it decays by tau_decay*prior for each following period until the end of time
  for (m in 1:max(fake_data$m + 1)) {
    if (m == 1) {
      tau_one[m] <- tau_one_zero / T_adjustment
    } else {
      tau_one[m] <- tau_one[m - 1] * tau_decay_adj
    }
  }
  # create a tibble that assings the treatment effect to each post treat period
  tau <- tibble(
    m = 0:max(fake_data$m),
    tau_one = tau_one
  )
  
  fake_data <- fake_data %>%
    mutate(
      new_pitches = 0,
      total_pitches = 0
    )
  
  # FIX: DOES NOT SCALE WELL
  # for each time and unit specific combination, we will create a variable for the number of
  # new pitches in that period (and a cumulative sum tracker)
  # no picthes happen until treatment, but a random amount happen every period after treatment
  tic("assign pitches")
  for (t in 1:N_time) {
    for (i in 1:N) {
      # if we are in the pre treatment periods (m<0), set pitches to 0
      if (fake_data$m[fake_data$division_country == i &
                      fake_data$time == t] < 0) {
        fake_data$new_pitches[fake_data$division_country == i &
                                fake_data$time == t] <- 0
      }
      else if (fake_data$m[fake_data$division_country == i &
                           fake_data$time == t] == 0) {
        # if we are in the treatment period
        # assign a random number of pitches (number of heads out of 10 coin flips) + minimum of 1
        fake_data$new_pitches[fake_data$division_country == i &
                                fake_data$time == t] <- 1 +
          ifelse(T_adjustment == 1, rbinom(n = 1, size = 10, prob = 0.5),
                 ifelse(T_adjustment == 4, rbinom(n = 1, size = 3, prob = 0.4),
                        rbernoulli(1, 2 / 9)
                 )
          )
        # if we are in the treatment period, total (cumulative pitches) is just the new pitched
        fake_data$total_pitches[fake_data$division_country == i &
                                  fake_data$time == t] <-
          fake_data$new_pitches[fake_data$division_country == i &
                                  fake_data$time == t]
      }
      else {
        # if in the post treat period, sample from the binomial without a min
        fake_data$new_pitches[fake_data$division_country == i &
                                fake_data$time == t] <- ifelse(T_adjustment == 1, rbinom(n = 1, size = 10, prob = 0.5),
                                                               ifelse(T_adjustment == 4, rbinom(n = 1, size = 3, prob = 0.4),
                                                                      rbernoulli(1, 2 / 9)
                                                               )
                                )
        # if we are in the proper post treat periods, we add the new pitches to the past cumulative sum
        fake_data$total_pitches[fake_data$division_country == i &
                                  fake_data$time == t] <-
          fake_data$total_pitches[fake_data$division_country == i &
                                    fake_data$time == t - 1] +
          fake_data$new_pitches[fake_data$division_country == i &
                                  fake_data$time == t]
      }
    }
  }
  toc()
  fake_data <- fake_data %>%
    mutate(total_pitches = case_when(
      treated == 1 ~ total_pitches,
      TRUE ~ 0
    ))
  
  
  
  # beta_i
  # define a random variable, specific to each unit
  beta <- tibble(
    division_country = 1:N,
    beta_i = rnorm(
      n = N, mean = 0,
      sd = 0.07789
    )
  )
  # add this random effect to the data
  fake_data <-
    inner_join(fake_data, beta, by = c("division_country"))
  
  
  
  # Initial y
  sd_e <- 0.5929 * (1 / T_adjustment**2)
  sd_0 <- 0.5929 / (sqrt(1 - rho_y^2))
  
  # Generate period 1 y from a normal distribution
  # centered around the sum of individual random effects (ai and beta_i),
  # vertical random effects (vertical_re), and time fixed effects (r)
  # y0 stores the unit and initial y value
  y0 <- fake_data %>%
    filter(time == 1) %>%
    mutate(
      mu_i_zero = log(exp(a_i + vertical_re + r + beta_i) / T_adjustment),
      y_0 = rnorm(n = N, mean = mu_i_zero, sd = sd_0)
    ) %>%
    select(division_country, y_0)
  
  # create an NxTime matrix with the values indicating time relative to first treat (m)
  m_matrix <- fake_data %>%
    select(division_country, time, m) %>%
    tidyr::pivot_wider(names_from = time, values_from = m) %>%
    select(-division_country) %>%
    as.matrix()
  
  # y matrix without treatment
  
  # create a NxTime matrix with values being the ith unit's mean in time j
  # the mean is the sum of random effects, which is where the normal y will be centered
  u <- fake_data %>%
    mutate(mu_i_t = log(exp(a_i + vertical_re + r + beta_i * (time / log(T_adjustment))) / T_adjustment)) %>%
    select(division_country, time, mu_i_t) %>%
    tidyr::pivot_wider(names_from = time, values_from = mu_i_t) %>%
    select(-division_country) %>%
    as.matrix()
  
  # initialize the non-treated outcome matrix, and assign the
  # first col (time 0, unit 1:N) from y0
  y_t0 <- matrix(nrow = N, ncol = N_time)
  y_t0[, 1] <- y0$y_0
  
  # generate the treated y from the same starting point
  y_t1 <- y_t0
  
  # generate a NxTime matrix, where value in i j is
  # unit i's total pitches in period j
  total_pitches_matrix <- fake_data %>%
    select(division_country, time, total_pitches) %>%
    tidyr::pivot_wider(names_from = time, values_from = total_pitches) %>%
    select(-division_country) %>%
    as.matrix()
  
  tic("generating target and counterfactual outcomes")
  # for each unit, by each time period
  for (i in 1:N) {
    for (t in 2:N_time) {
      # define the period t mean as the difference between
      # i;s last period outcome and their last period mean
      # (ie y_t0 has mean u, but is drawn from normal centered around u)
      # then scale that difference by rho, and add it to the current mean,
      # and draw a new outcome centered around this AR mean
      mu_it <- rho_y * (y_t0[i, (t - 1)] - u[i, (t - 1)]) +
        u[i, t]
      y_t0[i, t] <- rnorm(n = 1, mean = mu_it, sd = sd_e)
      
      # if the unit receives treatment in this period -- ie they are in the post period
      # and received treatment
      if (m_matrix[[i, t]] >= 0 & total_pitches_matrix[[i, t]] > 0) {
        # stop(total_pitches_matrix[i,t]<=0)
        # update only the treated matrix -- so that y_t0 is our counterfactual
        # update by adding the decayed impact along with
        # gamma (param) * the log of total pitches to date
        # Ask ignacio -- should we care about the marginal increase in pitches?
        y_t1[i, t] <- y_t0[i, t] +
          tau_one[m_matrix[i, t] + 1] +
          gamma * log(total_pitches_matrix[i, t])
      } else {
        # if untreated, use counterfactual outcome
        y_t1[i, t] <- y_t0[i, t]
      }
    }
  }
  toc()
  
  colnames(y_t1) <- colnames(y_t0) <- 1:N_time
  
  # transform the counterfactual and treated matrices of y into tibbles
  # and pivot long. This yields a tibble with each row having a unique combo of
  # time, unit and the affiliated y0 or y1
  y_t0_long <- as_tibble(y_t0) %>%
    mutate(division_country = 1:N) %>%
    tidyr::pivot_longer(names_to = "time", cols = -division_country) %>%
    mutate(time = as.integer(time)) %>%
    rename(y0 = value)
  
  y_t1_long <- as_tibble(y_t1) %>%
    mutate(division_country = 1:N) %>%
    tidyr::pivot_longer(names_to = "time", cols = -division_country) %>%
    mutate(time = as.integer(time)) %>%
    rename(y1 = value)
  
  
  # merge the treated and counterfactual
  y_long <-
    inner_join(y_t0_long, y_t1_long, by = c("division_country", "time"))
  
  # merge the treated/counterfactual outcomes into fake_data
  # exponentiate all outcomes, and multiply by 1000
  # define y as the "observed" outcome in potential outcome framework
  fake_data <-
    inner_join(fake_data, y_long, by = c("division_country", "time")) %>%
    mutate(y = case_when(
      treated == 1 ~ y1,
      treated == 0 ~ y0
    )) %>%
    mutate(
      y0 = exp(y0) * 1000,
      y1 = exp(y1) * 1000,
      y = exp(y) * 1000
    )
  
  # True lift
  # q1 is the first quarter after treatment -- post treat period 0,1,2
  q1 <- fake_data %>%
    filter(
      m %in% unique(c(0, 1, 2, (T_adjustment > 1) * seq(3, 3 * T_adjustment))),
      treated == 1
    )
  
  # aggregate Y over that period for all indiv, time (and the counterfactual version)
  sum_y0 <- sum(q1$y0)
  sum_y1 <- sum(q1$y1)
  
  # aggregate percent lift in q1
  true_lift <- (sum_y1 / sum_y0) - 1
  
  ## Output data
  # just the x variables, removing v25 (the randomly effected based on verticals)
  xs <- xs %>%
    select(-v25)
  # Tibble of the relevant output data
  # synthetic_data <- fake_data %>%
  #   dplyr::select(division_country,
  #                 time,
  #                 treated,
  #                 country,
  #                 vertical,
  #                 total_pitches,
  #                 y)  %>%
  #   inner_join(xs, by = "division_country")
  
  # theoretically possible to receive 0 pitches on the first potential treat period
  # Fix: adjust First
  synthetic_data_full <- fake_data %>%
    dplyr::select(-y1) %>%
    rename(
      period = time, post_treat_t = m, Treatment_Period = first,
      counterfactual = y0, target = y
    ) %>%
    inner_join(xs, by = "division_country") %>%
    rename(entry = division_country) %>%
    arrange(period, entry)
  
  
  
  # my_fake_data <- synthetic_data_full %>% mutate(
  #   true_lift = true_lift,
  #   N = N,
  #   tau_one_zero = tau_one_zero,
  #   gamma = gamma,
  #   seed = seed,
  #   prop_treated = mean(synthetic_data_full$treated),
  #   rho_y = rho_y
  # ) %>%
  #   tidyr::nest(data = -c(true_lift, tau_one_zero, gamma, N, seed, prop_treated, rho_y))
}


aa_daily_test <- gen_data_lubr(T_freq = "daily", N = 50)
aa_monthly_test <- gen_data_lubr(N=50)

#Need to recover or recreate the aggregate time var
aa_monthly_test2 <- aggregate_data_by_date(aa_daily_test)
aa_monthly_test2 %>%
  group_by(period) %>%
  summarize(mean_cf = mean(y0), agg_cf = sum(y0))
aa_monthly_test %>%
  group_by(period) %>%
  summarize(mean_cf = mean(y0), agg_cf = sum(y0))








#Old benchmarking

####################################################3
#SELECTION VS NO SELECTION
####################################################3
seeds=c(1,2)
tic("Starting DGP")
AA_data_no_sel_base=factor_synthetic_dgp(date_start="2010-01-01",
                                         first_treat="2017-07-01",
                                         date_end="2020-01-01",
                                         num_entries=200,
                                         prop_treated=0.25,
                                         treat_impact_sd = 0, 
                                         treat_impact_mean = 0, 
                                         rho=0.9,
                                         rescale_y_mean = 2.5e3,
                                         cov_overlap_scale = 0,
                                         seed=42)

AA_data_no_sel_unformatted=furrr::future_map(.x=seeds, 
                                             .f=~noisify_draw(data_inp=
                                                                AA_data_no_sel_base,
                                                              seed=.x))

AA_data_no_sel=future_map(AA_data_no_sel_unformatted, format_for_est)
toc()

tic("Starting DGP, Selection")
AA_data_sel_base=factor_synthetic_dgp(date_start="2010-01-01",
                                      first_treat="2017-07-01",
                                      date_end="2020-01-01",
                                      num_entries=200,
                                      prop_treated=0.25,
                                      treat_impact_sd = 0, 
                                      treat_impact_mean = 0, 
                                      rho=0.9,
                                      rescale_y_mean = 2.5e3,
                                      cov_overlap_scale = 0,
                                      intercept_scale = 0.75,
                                      seed=42)


AA_data_sel_unformatted=furrr::future_map(.x=seeds, 
                                          .f=~noisify_draw(data_inp=
                                                             AA_data_sel_base,
                                                           seed=.x))

AA_data_sel=future_map(AA_data_sel_unformatted, format_for_est)
toc()



tic("Estimating Gsynth, Placebo, No Selection")
gsynth_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se=F)
gsynth_tot_AA=future_map(gsynth_AA, compute_tot_se_jackknife, stat_in="mean")

gsynth_AA_bias=compute_jackknife_bias(gsynth_tot_AA)
gsynth_bias_plot=create_gap_ci_plot(gsynth_AA_bias, 
                                    plot_title="Gsynth Rho=0.1, Bias", 
                                    plot_x_lab="Post-Treat Time",
                                    plot_y_lab="ATT Bias", pct_flag = F)
gsynth_overall_metrics=compute_jackknife_metrics(gsynth_AA)
gsynth_AA_tot_var=compute_tot_variance(gsynth_tot_AA)
gsynth_AA_cov=compute_tot_coverage(gsynth_tot_AA)
toc()



tic("Estimating Gsynth, Placebo, Selection")
gsynth_AA_sel=future_map(AA_data_sel, estimate_gsynth_series, se_est=F)
gsynth_tot_AA_sel=future_map(gsynth_AA_sel, compute_tot_se_jackknife,
                             stat_in="mean")

gsynth_AA_sel_bias=compute_jackknife_bias(gsynth_tot_AA_sel)
gsynth_sel_bias_plot=create_gap_ci_plot(gsynth_AA_sel_bias, 
                                        plot_title="Gsynth w/Selection - Rho=0.1, Bias", 
                                        plot_x_lab="Post-Treat Time",
                                        plot_y_lab="ATT Bias", pct_flag = F)
gsynth_sel_overall_metrics=compute_jackknife_metrics(gsynth_AA_sel)
gsynth_sel_AA_tot_var=compute_tot_variance(gsynth_tot_AA_sel)
gsynth_AA_sel_cov=compute_tot_coverage(gsynth_tot_AA_sel )
toc()




tic("Estimating SCDID")
scdid_AA=future_map(AA_data_no_sel, estimate_scdid_series)
scdid_tot_AA=future_map(scdid_AA, compute_tot_se_jackknife, stat_in="mean")

scdid_AA_bias=compute_jackknife_bias(scdid_tot_AA)
scdid_AA_bias_plot=create_gap_ci_plot(scdid_AA_bias, 
                                      plot_title="SCDID Bias", 
                                      plot_x_lab="Post-Treat Time",
                                      plot_y_lab="ATT Bias", pct_flag = F)
scdid_AA_overall_metrics=compute_jackknife_metrics(scdid_AA)
scdid_AA_tot_var=compute_tot_variance(scdid_tot_AA)
scdid_AA_cov=compute_tot_coverage(scdid_tot_AA)
toc()

tic("Estimating SCDID Sel")
scdid_AA_sel=future_map(AA_data_sel, estimate_scdid_series)
scdid_tot_AA_sel=future_map(scdid_AA_sel, compute_tot_se_jackknife, stat_in="mean")

scdid_AA_sel_bias=compute_jackknife_bias(scdid_tot_AA_sel)
scdid_AA_sel_bias_plot=create_gap_ci_plot(scdid_AA_sel_bias, 
                                          plot_title="SCDID w/Selection - Bias", 
                                          plot_x_lab="Post-Treat Time",
                                          plot_y_lab="ATT Bias", pct_flag = F)
scdid_AA_sel_overall_metrics=compute_jackknife_metrics(scdid_AA_sel)
scdid_AA_sel_tot_var=compute_tot_variance(scdid_tot_AA_sel)
scdid_AA_sel_cov=compute_tot_coverage(scdid_tot_AA_sel)
toc()

tic("Estimating MC")

mc_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se_est=F,
                 estimator_type="mc")
mc_tot_AA=future_map(mc_AA, compute_tot_se_jackknife, stat_in="mean")

mc_AA_bias=compute_jackknife_bias(mc_tot_AA)
mc_AA_bias_plot=create_gap_ci_plot(mc_AA_bias, 
                                   plot_title="MC - Bias", 
                                   plot_x_lab="Post-Treat Time",
                                   plot_y_lab="ATT Bias",  pct_flag = F)
mc_AA_overall_metrics=compute_jackknife_metrics(mc_AA)
mc_AA_tot_var=compute_tot_variance(mc_tot_AA)
mc_AA_cov=compute_tot_coverage(mc_tot_AA)
toc()




tic("Estimating MC sel")

mc_AA_sel=future_map(AA_data_sel, estimate_gsynth_series, se_est=F,
                     estimator_type="mc")
mc_tot_AA_sel=future_map(mc_AA_sel, compute_tot_se_jackknife, stat_in="mean")

mc_AA_sel_bias=compute_jackknife_bias(mc_tot_AA_sel)
mc_AA_sel_bias_plot=create_gap_ci_plot(mc_AA_sel_bias, 
                                       plot_title="MC w/Selection - Bias", 
                                       plot_x_lab="Post-Treat Time",
                                       plot_y_lab="ATT Bias",  pct_flag = F)
mc_AA_sel_overall_metrics=compute_jackknife_metrics(mc_AA_sel)
mc_AA_sel_tot_var=compute_tot_variance(mc_tot_AA_sel)
mc_AA_sel_cov=compute_tot_coverage(mc_tot_AA_sel)
toc()

save.image(here("Data/Cloud_AA_50.RData"))








####################################################3
#SELECTION VS NO SELECTION LOW RHO
####################################################3
seeds=c(1,2)
tic("Starting DGP")
AA_data_no_sel_base=factor_synthetic_dgp(date_start="2010-01-01",
                                         first_treat="2017-07-01",
                                         date_end="2020-01-01",
                                         num_entries=200,
                                         prop_treated=0.25,
                                         treat_impact_sd = 0, 
                                         treat_impact_mean = 0, 
                                         rho=0.1,
                                         rescale_y_mean = 2.5e3,
                                         cov_overlap_scale = 0,
                                         seed=42)

AA_data_no_sel_unformatted=furrr::future_map(.x=seeds, 
                                             .f=~noisify_draw(data_inp=
                                                                AA_data_no_sel_base,
                                                              seed=.x))

AA_data_no_sel=future_map(AA_data_no_sel_unformatted, format_for_est)
toc()

tic("Starting DGP, Selection")
AA_data_sel_base=factor_synthetic_dgp(date_start="2010-01-01",
                                      first_treat="2017-07-01",
                                      date_end="2020-01-01",
                                      num_entries=200,
                                      prop_treated=0.25,
                                      treat_impact_sd = 0, 
                                      treat_impact_mean = 0, 
                                      rho=0.1,
                                      rescale_y_mean = 2.5e3,
                                      cov_overlap_scale = 0,
                                      intercept_scale = 0.75,
                                      seed=42)


AA_data_sel_unformatted=furrr::future_map(.x=seeds, 
                                          .f=~noisify_draw(data_inp=
                                                             AA_data_sel_base,
                                                           seed=.x))

AA_data_sel=future_map(AA_data_sel_unformatted, format_for_est)
toc()



tic("Estimating Gsynth, Placebo, No Selection")
gsynth_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se=F, criterion="pc")
gsynth_tot_AA=future_map(gsynth_AA, compute_tot_se_jackknife, stat_in="mean")

gsynth_AA_bias=compute_jackknife_bias(gsynth_tot_AA)
gsynth_bias_plot=create_gap_ci_plot(gsynth_AA_bias, 
                                    plot_title="Gsynth Rho=0.1, Bias", 
                                    plot_x_lab="Post-Treat Time",
                                    plot_y_lab="ATT Bias", pct_flag = F)
gsynth_overall_metrics=compute_jackknife_metrics(gsynth_AA)
gsynth_AA_tot_var=compute_tot_variance(gsynth_tot_AA)
gsynth_AA_cov=compute_tot_coverage(gsynth_tot_AA)
toc()



tic("Estimating Gsynth, Placebo, Selection")
gsynth_AA_sel=future_map(AA_data_sel, estimate_gsynth_series, se_est=F)
gsynth_tot_AA_sel=future_map(gsynth_AA_sel, compute_tot_se_jackknife,
                             stat_in="mean")

gsynth_AA_sel_bias=compute_jackknife_bias(gsynth_tot_AA_sel)
gsynth_sel_bias_plot=create_gap_ci_plot(gsynth_AA_sel_bias, 
                                        plot_title="Gsynth w/Selection - Rho=0.1, Bias", 
                                        plot_x_lab="Post-Treat Time",
                                        plot_y_lab="ATT Bias", pct_flag = F)
gsynth_sel_overall_metrics=compute_jackknife_metrics(gsynth_AA_sel)
gsynth_sel_AA_tot_var=compute_tot_variance(gsynth_tot_AA_sel)
gsynth_AA_sel_cov=compute_tot_coverage(gsynth_tot_AA_sel )
toc()




tic("Estimating SCDID")
scdid_AA=future_map(AA_data_no_sel, estimate_scdid_series)
scdid_tot_AA=future_map(scdid_AA, compute_tot_se_jackknife, stat_in="mean")

scdid_AA_bias=compute_jackknife_bias(scdid_tot_AA)
scdid_AA_bias_plot=create_gap_ci_plot(scdid_AA_bias, 
                                      plot_title="SCDID Rho=0.1, Bias", 
                                      plot_x_lab="Post-Treat Time",
                                      plot_y_lab="ATT Bias", pct_flag = F)
scdid_AA_overall_metrics=compute_jackknife_metrics(scdid_AA)
scdid_AA_tot_var=compute_tot_variance(scdid_tot_AA)
scdid_AA_cov=compute_tot_coverage(scdid_tot_AA)
toc()

tic("Estimating SCDID Sel")
scdid_AA_sel=future_map(AA_data_sel, estimate_scdid_series)
scdid_tot_AA_sel=future_map(scdid_AA_sel, compute_tot_se_jackknife, stat_in="mean")

scdid_AA_sel_bias=compute_jackknife_bias(scdid_tot_AA_sel)
scdid_AA_sel_bias_plot=create_gap_ci_plot(scdid_AA_sel_bias, 
                                          plot_title="SCDID w/Selection - Rho=0.1, Bias", 
                                          plot_x_lab="Post-Treat Time",
                                          plot_y_lab="ATT Bias", pct_flag = F)
scdid_AA_sel_overall_metrics=compute_jackknife_metrics(scdid_AA_sel)
scdid_AA_sel_tot_var=compute_tot_variance(scdid_tot_AA_sel)
scdid_AA_sel_cov=compute_tot_coverage(scdid_tot_AA_sel)
toc()

tic("Estimating MC")

mc_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se_est=F,
                 estimator_type="mc")
mc_tot_AA=future_map(mc_AA, compute_tot_se_jackknife, stat_in="mean")

mc_AA_bias=compute_jackknife_bias(mc_tot_AA)
mc_AA_bias_plot=create_gap_ci_plot(mc_AA_bias, 
                                   plot_title="MC - Rho=0.1, Bias", 
                                   plot_x_lab="Post-Treat Time",
                                   plot_y_lab="ATT Bias",  pct_flag = F)
mc_AA_overall_metrics=compute_jackknife_metrics(mc_AA)
mc_AA_tot_var=compute_tot_variance(mc_tot_AA)
mc_AA_cov=compute_tot_coverage(mc_tot_AA)
toc()




tic("Estimating MC sel")

mc_AA_sel=future_map(AA_data_sel, estimate_gsynth_series, se_est=F,
                     estimator_type="mc")
mc_tot_AA_sel=future_map(mc_AA_sel, compute_tot_se_jackknife, stat_in="mean")

mc_AA_sel_bias=compute_jackknife_bias(mc_tot_AA_sel)
mc_AA_sel_bias_plot=create_gap_ci_plot(mc_AA_sel_bias, 
                                       plot_title="MC w/Selection - Rho=0.1, Bias", 
                                       plot_x_lab="Post-Treat Time",
                                       plot_y_lab="ATT Bias",  pct_flag = F)
mc_AA_sel_overall_metrics=compute_jackknife_metrics(mc_AA_sel)
mc_AA_sel_tot_var=compute_tot_variance(mc_tot_AA_sel)
mc_AA_sel_cov=compute_tot_coverage(mc_tot_AA_sel)
toc()

save.image(here("Data/Cloud_AA_ACF_50.RData"))




####################################################3
#SELECTION VS NO SELECTION (ON LOADINGS)
####################################################3
tic("Starting DGP")
AA_data_no_sel_base=factor_synthetic_dgp(date_start="2010-01-01",
                                         first_treat="2017-07-01",
                                         date_end="2020-01-01",
                                         num_entries=200,
                                         prop_treated=0.25,
                                         treat_impact_sd = 0, 
                                         treat_impact_mean = 0, 
                                         rho=0.9,
                                         rescale_y_mean = 2.5e3,
                                         cov_overlap_scale = 0,
                                         seed=42)

AA_data_no_sel_unformatted=furrr::future_map(.x=seeds, 
                                             .f=~noisify_draw(data_inp=
                                                                AA_data_no_sel_base,
                                                              seed=.x))

AA_data_no_sel=future_map(AA_data_no_sel_unformatted, format_for_est)
toc()

tic("Starting DGP, Selection")
AA_data_sel_base=factor_synthetic_dgp(date_start="2010-01-01",
                                      first_treat="2017-07-01",
                                      date_end="2020-01-01",
                                      num_entries=200,
                                      prop_treated=0.25,
                                      treat_impact_sd = 0, 
                                      treat_impact_mean = 0, 
                                      rho=0.9,
                                      rescale_y_mean = 2.5e3,
                                      cov_overlap_scale = 0,
                                      loading_scale = 0.75,
                                      seed=42)


AA_data_sel_unformatted=furrr::future_map(.x=seeds, 
                                          .f=~noisify_draw(data_inp=
                                                             AA_data_sel_base,
                                                           seed=.x))

AA_data_sel=future_map(AA_data_sel_unformatted, format_for_est)
toc()



tic("Estimating Gsynth, Placebo, No Selection")
gsynth_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se_est=F)
gsynth_tot_AA=future_map(gsynth_AA, compute_tot_se_jackknife, stat_in="mean")

gsynth_AA_bias=compute_jackknife_bias(gsynth_tot_AA)
gsynth_bias_plot=create_gap_ci_plot(gsynth_AA_bias, 
                                    plot_title="Gsynth, Bias", 
                                    plot_x_lab="Post-Treat Time",
                                    plot_y_lab="ATT Bias", pct_flag = F)
gsynth_overall_metrics=compute_jackknife_metrics(gsynth_AA)
gsynth_AA_tot_var=compute_tot_variance(gsynth_tot_AA)
gsynth_AA_cov=compute_tot_coverage(gsynth_tot_AA)
toc()



tic("Estimating Gsynth, Placebo, Selection")
gsynth_AA_sel=future_map(AA_data_sel, estimate_gsynth_series, se_est=F)
gsynth_tot_AA_sel=future_map(gsynth_AA_sel, compute_tot_se_jackknife,
                             stat_in="mean")

gsynth_AA_sel_bias=compute_jackknife_bias(gsynth_tot_AA_sel)
gsynth_sel_bias_plot=create_gap_ci_plot(gsynth_AA_sel_bias, 
                                        plot_title="Gsynth w/Selection on Loadings, Bias", 
                                        plot_x_lab="Post-Treat Time",
                                        plot_y_lab="ATT Bias", pct_flag = F)
gsynth_sel_overall_metrics=compute_jackknife_metrics(gsynth_AA_sel)
gsynth_sel_AA_tot_var=compute_tot_variance(gsynth_tot_AA_sel)
gsynth_AA_sel_cov=compute_tot_coverage(gsynth_tot_AA_sel )
toc()




tic("Estimating SCDID")
scdid_AA=future_map(AA_data_no_sel, estimate_scdid_series)
scdid_tot_AA=future_map(scdid_AA, compute_tot_se_jackknife, stat_in="mean")

scdid_AA_bias=compute_jackknife_bias(scdid_tot_AA)
scdid_AA_bias_plot=create_gap_ci_plot(scdid_AA_bias, 
                                      plot_title="SCDID Bias", 
                                      plot_x_lab="Post-Treat Time",
                                      plot_y_lab="ATT Bias", pct_flag = F)
scdid_AA_overall_metrics=compute_jackknife_metrics(scdid_AA)
scdid_AA_tot_var=compute_tot_variance(scdid_tot_AA)
scdid_AA_cov=compute_tot_coverage(scdid_tot_AA)
toc()

tic("Estimating SCDID Sel")
scdid_AA_sel=future_map(AA_data_sel, estimate_scdid_series)
scdid_tot_AA_sel=future_map(scdid_AA_sel, compute_tot_se_jackknife, stat_in="mean")

scdid_AA_sel_bias=compute_jackknife_bias(scdid_tot_AA_sel)
scdid_AA_sel_bias_plot=create_gap_ci_plot(scdid_AA_sel_bias, 
                                          plot_title="SCDID w/Selection  on Loadings, Bias", 
                                          plot_x_lab="Post-Treat Time",
                                          plot_y_lab="ATT Bias", pct_flag = F)
scdid_AA_sel_overall_metrics=compute_jackknife_metrics(scdid_AA_sel)
scdid_AA_sel_tot_var=compute_tot_variance(scdid_tot_AA_sel)
scdid_AA_sel_cov=compute_tot_coverage(scdid_tot_AA_sel)
toc()

tic("Estimating MC")

mc_AA=future_map(AA_data_no_sel, estimate_gsynth_series, se_est=F,
                 estimator_type="mc")
mc_tot_AA=future_map(mc_AA, compute_tot_se_jackknife, stat_in="mean")

mc_AA_bias=compute_jackknife_bias(mc_tot_AA)
mc_AA_bias_plot=create_gap_ci_plot(mc_AA_bias, 
                                   plot_title="MC Bias", 
                                   plot_x_lab="Post-Treat Time",
                                   plot_y_lab="ATT Bias",  pct_flag = F)
mc_AA_overall_metrics=compute_jackknife_metrics(mc_AA)
mc_AA_tot_var=compute_tot_variance(mc_tot_AA)
mc_AA_cov=compute_tot_coverage(mc_tot_AA)
toc()




tic("Estimating MC sel")

mc_AA_sel=future_map(AA_data_sel, estimate_gsynth_series, se_est=F,
                     estimator_type="mc")
mc_tot_AA_sel=future_map(mc_AA_sel, compute_tot_se_jackknife, stat_in="mean")

mc_AA_sel_bias=compute_jackknife_bias(mc_tot_AA_sel)
mc_AA_sel_bias_plot=create_gap_ci_plot(mc_AA_sel_bias, 
                                       plot_title="MC w/Selection on Loadings, Bias", 
                                       plot_x_lab="Post-Treat Time",
                                       plot_y_lab="ATT Bias",  pct_flag = F)
mc_AA_sel_overall_metrics=compute_jackknife_metrics(mc_AA_sel)
mc_AA_sel_tot_var=compute_tot_variance(mc_tot_AA_sel)
mc_AA_sel_cov=compute_tot_coverage(mc_tot_AA_sel)
toc()

save.image(here("Data/Cloud_AA_Loadings_50.RData"))



#Load Data 1 (hi ACF), run CI
load(here("Data/Cloud_AA_50.RData"))
tic("Estimating CI 1")

causalimpact_AA=future_map(AA_data_no_sel, estimate_causalimpact_series)
causalimpact_tot_AA=future_map(causalimpact_AA, compute_tot_se_jackknife, stat_in="mean")

causalimpact_AA_bias=compute_jackknife_bias(causalimpact_tot_AA)
causalimpact_AA_bias_plot=create_gap_ci_plot(causalimpact_AA_bias, 
                                             plot_title="Causal Impact Bias", 
                                             plot_x_lab="Post-Treat Time",
                                             plot_y_lab="ATT Bias",  pct_flag = F)
causalimpact_AA_overall_metrics=compute_jackknife_metrics(causalimpact_AA)
causalimpact_AA_tot_var=compute_tot_variance(causalimpact_tot_AA)
causalimpact_AA_cov=compute_tot_coverage(causalimpact_tot_AA)
toc()

tic("Estimating CI 1, sel")

causalimpact_AA_sel=future_map(AA_data_sel, estimate_causalimpact_series)
causalimpact_tot_AA_sel=future_map(causalimpact_AA_sel, compute_tot_se_jackknife, stat_in="mean")

causalimpact_AA_sel_bias=compute_jackknife_bias(causalimpact_tot_AA_sel)
causalimpact_AA_sel_bias_plot=create_gap_ci_plot(causalimpact_AA_sel_bias, 
                                                 plot_title="Causal Impact Bias", 
                                                 plot_x_lab="Post-Treat Time",
                                                 plot_y_lab="ATT Bias",  pct_flag = F)
causalimpact_AA_sel_overall_metrics=compute_jackknife_metrics(causalimpact_AA_sel)
causalimpact_AA_sel_tot_var=compute_tot_variance(causalimpact_tot_AA_sel)
causalimpact_AA_sel_cov=compute_tot_coverage(causalimpact_tot_AA_sel)
toc()
save.image(here("Data/Cloud_AA_50_v2.RData"))


#Load Data 2 (low acf)
load(here("Data/Cloud_AA_ACF_50.RData"))
tic("Estimating CI 2")

causalimpact_AA=future_map(AA_data_no_sel, estimate_causalimpact_series)
causalimpact_tot_AA=future_map(causalimpact_AA, compute_tot_se_jackknife, stat_in="mean")

causalimpact_AA_bias=compute_jackknife_bias(causalimpact_tot_AA)
causalimpact_AA_bias_plot=create_gap_ci_plot(causalimpact_AA_bias, 
                                             plot_title="Causal Impact Bias", 
                                             plot_x_lab="Post-Treat Time",
                                             plot_y_lab="ATT Bias",  pct_flag = F)
causalimpact_AA_overall_metrics=compute_jackknife_metrics(causalimpact_AA)
causalimpact_AA_tot_var=compute_tot_variance(causalimpact_tot_AA)
causalimpact_AA_cov=compute_tot_coverage(causalimpact_tot_AA)
toc()

tic("Estimating CI 2, sel")

causalimpact_AA_sel=future_map(AA_data_sel, estimate_causalimpact_series)
causalimpact_tot_AA_sel=future_map(causalimpact_AA_sel, compute_tot_se_jackknife, stat_in="mean")

causalimpact_AA_sel_bias=compute_jackknife_bias(causalimpact_tot_AA_sel)
causalimpact_AA_sel_bias_plot=create_gap_ci_plot(causalimpact_AA_sel_bias, 
                                                 plot_title="Causal Impact Bias", 
                                                 plot_x_lab="Post-Treat Time",
                                                 plot_y_lab="ATT Bias",  pct_flag = F)
causalimpact_AA_sel_overall_metrics=compute_jackknife_metrics(causalimpact_AA_sel)
causalimpact_AA_sel_tot_var=compute_tot_variance(causalimpact_tot_AA_sel)
causalimpact_AA_sel_cov=compute_tot_coverage(causalimpact_tot_AA_sel)
toc()
save.image(here("Data/Cloud_AA_ACF_50_v2.RData"))


#Load Data 3 (loading selection)
load(here("Data/Cloud_AA_Loadings_50.RData"))
tic("Estimating CI 3")

causalimpact_AA=future_map(AA_data_no_sel, estimate_causalimpact_series)
causalimpact_tot_AA=future_map(causalimpact_AA, compute_tot_se_jackknife, stat_in="mean")

causalimpact_AA_bias=compute_jackknife_bias(causalimpact_tot_AA)
causalimpact_AA_bias_plot=create_gap_ci_plot(causalimpact_AA_bias, 
                                             plot_title="Causal Impact Bias", 
                                             plot_x_lab="Post-Treat Time",
                                             plot_y_lab="ATT Bias",  pct_flag = F)
causalimpact_AA_overall_metrics=compute_jackknife_metrics(causalimpact_AA)
causalimpact_AA_tot_var=compute_tot_variance(causalimpact_tot_AA)
causalimpact_AA_cov=compute_tot_coverage(causalimpact_tot_AA)
toc()

tic("Estimating CI 3, sel")

causalimpact_AA_sel=future_map(AA_data_sel, estimate_causalimpact_series)
causalimpact_tot_AA_sel=future_map(causalimpact_AA_sel, compute_tot_se_jackknife, stat_in="mean")

causalimpact_AA_sel_bias=compute_jackknife_bias(causalimpact_tot_AA_sel)
causalimpact_AA_sel_bias_plot=create_gap_ci_plot(causalimpact_AA_sel_bias, 
                                                 plot_title="Causal Impact Bias", 
                                                 plot_x_lab="Post-Treat Time",
                                                 plot_y_lab="ATT Bias",  pct_flag = F)
causalimpact_AA_sel_overall_metrics=compute_jackknife_metrics(causalimpact_AA_sel)
causalimpact_AA_sel_tot_var=compute_tot_variance(causalimpact_tot_AA_sel)
causalimpact_AA_sel_cov=compute_tot_coverage(causalimpact_tot_AA_sel)
toc()
save.image(here("Data/Cloud_AA_Loadings_50_v2.RData"))
