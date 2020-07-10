# Try to write functions for everything:

# ingest data can be a function

# Run the TE can be a function of
# AA or AB (which type of synthetic analysis is the goal)
# How many types of data do you want to run?
# is there a way to do this in parallel??
# maybe want to create helper functions so that we estimate (and return) Gsynth separately from
# CausalImpact, and SCDID (once I have the proper one)
# That way, I can call the furr r package and map each dataset to the gsynth_estimation etc.

# can then create print plot methods
# to print the Density of ATT, the Gap plots, as well as printing tables


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


############################################
# Helper functions for AA testing: selecting the placebo treated

# Potential Improvements/Options:
# First -- find K possible matches for each treated, without replacement.
# If there are ties, break them based on how good the second best is,

# Second -- threshold cutoff? If the closest placebo-treated is far,
# just throw out that treated.
###########################################


nearest_ts_euclidean <- function(ts_tomatch, ts_rest) {
  # Helper function to matching_without_replacement
  # Measures the euclidean distance between ts_tomatch
  # and each of the rows in the dataframe ts_rest.

  # Args
  # ts_tomatch: treated time series dataframe, rows identify ID and columns Time
  # ts_rest: donor time series dataframe, rows identify ID and columns Time.

  # Output
  # Vector of indices indicating the row in ts_rest that are best match
  # (min L2 norm distance) to the row in ts_tomatch.
  # (Vector length equal to number of rows in ts_tomatch, ie treated entries)
  # If ts_tomatch has more than 1 entry, the matching is done with replacement
  # (multiple treated can have the same match).

  apply(ts_tomatch, 1, function(ts_tomatch) {
    which.min(
      apply(
        ts_rest, 1, function(ts_rest, ts_tomatch) {
          dist(rbind(ts_rest, ts_tomatch))
        },
        ts_tomatch
      )
    )
  })
}


matching_without_replacement <- function(treated_block, control_block, id_var, treat_period) {
  # finds the nearest match for each unit in the treated subsample (treated_block)
  # among the donor pool (control_block) without replacement via individual calls
  # to nearest_ts_euclidean.

  # Args
  # treated_block: treated time series dataframe, rows identify ID and columns Time
  # control_block: donor time series dataframe, rows identify ID and columns Time.

  # Output
  # df_toreturn, a dataframe containing a column for the placebo-treated unit ID numbers,
  # the treated unit it was the nearest match to, and the time that treated unit was actually
  # treated (num_rows of the dataframe equal to num_rows of treated_block)

  # store an empty vector for the donor IDs that match
  already_matched <- c()
  # Store the time of treatment and treated ID for the true treated units
  placebo_treat_period <- treated_block %>% pull(!!as.name(treat_period))
  treatment_unit <- treated_block %>% pull(!!as.name(id_var))
  for (i in seq_len(nrow(treated_block))) {
    if (i == 1) {
      # if we are searching for the match of our first treated unit, we can search across all donors
      temp_match <- nearest_ts_euclidean(
        treated_block %>%
          slice(i) %>%
          dplyr::select(-c(
            !!as.name(id_var),
            !!as.name(treat_period)
          )),
        control_block %>%
          dplyr::select(-!!as.name(id_var))
      )
      already_matched[i] <- control_block %>%
        slice(temp_match) %>%
        pull(!!as.name(id_var))
    }

    if (i != 1) {
      # If we have already found a match, restrict the search for future matches to the subset of currently unmatched donors
      temp_match <- nearest_ts_euclidean(treated_block %>%
                                           slice(i) %>%
                                           dplyr::select(-c(!!as.name(id_var), !!as.name(treat_period))),
                                         control_block %>%
        filter(!!as.name(id_var) %in% setdiff(!!as.name(id_var), already_matched)) %>% 
          dplyr::select(-!!as.name(id_var)))

      already_matched[i] <- control_block %>%
        filter(!!as.name(id_var) %in% setdiff(!!as.name(id_var), already_matched)) %>%
        slice(temp_match) %>%
        pull(!!as.name(id_var))
    }
  }
  # Store the resulting vectors in a dataframe for output
  df_toreturn <- data.frame(temp_id = already_matched, temp_treattime = placebo_treat_period, Treatment_Unit = treatment_unit)
  df_toreturn <- df_toreturn %>% rename(!!as.name(id_var) := temp_id, !!as.name(treat_period) := temp_treattime)
  return(df_toreturn)
}



# To-Do: expand for the case when we have covariates? This could be useful in 1) matching the placebo-treated,
# 2) improving performance and range of tests available to use

# To-Do: if covariates are allowed, must adapt functions (pivot wide/long, specifically)
create_placebo_df <- function(data_full, id_var = "entry",
                              time_var = "period", treat_indicator = "treatperiod_0",
                              outcome_var = "target", counterfac_var = "counter_factual") {
  # Generates a placebo only dataframe, using matching methods to select
  # placebo-treated entries as those most similar to truly-treated.

  # Args
  # data_full: long-form dataframe with both treated and control entries.
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_indicator:column name of binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)

  # Higher level description of data_full:
  # rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  # each row should have a treatment indicator (treat_indicator), a period number (time_var),
  # an individual ID (id_var), and an outcome (outcome_var)
  # for associated with that period-ID combination


  # Output
  # placebo_df_long, a dataframe of the same format as data_full, but now entirely
  # consists of donor units, some of which are placebo-treated (based on matching).



  # Split the dataset based on whether they are ever treated
  tr_entries <- data_full %>%
    filter(!!as.name(treat_indicator) > 0) %>%
    distinct(!!as.name(id_var)) %>%
    pull()
  ct_entries <- setdiff(data_full %>%
    distinct(!!as.name(id_var)) %>%
    pull(), tr_entries)
  # Create a dataframe of the subset of control units
  cd <- data_full %>% filter(!!as.name(id_var) %in% ct_entries)
  cd <- merge(
    cd,
    data.frame(
      entry = ct_entries,
      rank = seq_along(ct_entries)
    )
  )

  # Pivot the long cd dataframe to wide..now, each row represents an id_var, with columns for the outcome_var at each time period
  cd_for_match <- pivot_wider(
    data = cd %>%
      arrange(!!as.name(time_var), !!as.name(id_var)),
    names_from = !!as.name(time_var),
    id_cols = c(!!as.name(id_var)),
    values_from = c(!!as.name(outcome_var))
  )

  # After identifyin the treated observations (for now, one by one), want to find the placebo match
  treated_to_match <- data_full %>%
    filter(!!as.name(id_var) %in% tr_entries)

  # We will want to make the data wide???
  # want a treatment data indicator before doing this as well
  treated_to_match <- treated_to_match %>%
    group_by(!!as.name(id_var)) %>%
    mutate(
      Treatment_Period =
        length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1
    ) %>%
    ungroup()

  # for pivoting, potential issues arise if we have several time varying covariates
  # we'd have to take the values_from each of them, and for any constant args we'd presumably have to add them to id_cols
  data_wide_m <- pivot_wider(
    data = treated_to_match,
    names_from = !!as.name(time_var),
    id_cols = c(!!as.name(id_var), Treatment_Period),
    values_from = c(!!as.name(outcome_var))
  ) # %>% as.matrix()

  # Supply the treated units for matching (in wide form) and the donor units (in wide form) to the match function
  matched_placebo_df_temp <- matching_without_replacement(data_wide_m, cd_for_match, id_var, "Treatment_Period")

  # We now have the set of control units that form the placebo set, along with their placebo Treat Period and their correspondning Treat Unit
  # Merge this DF into the cd_for_match using id_var as the key

  placebo_df_wide <- cd_for_match %>%
    left_join(matched_placebo_df_temp, by = id_var)

  # Pivot_long it and recreate the treat_indicator indicator based on Treatment_Period, so that the ouput matches the input
  # If there are covariates or other variables which we do not want to be pivoting, append them to non_time_vars vectors
  # may have to add values_to as well if we have multiple time varying variables
  non_time_vars <- c(id_var, "Treatment_Period", "Treatment_Unit")
  placebo_df_long <- placebo_df_wide %>%
    pivot_longer(-all_of(non_time_vars),
      names_to = time_var, values_to = outcome_var
    ) %>%
    mutate(!!as.name(time_var) := as.numeric(!!as.name(time_var))) %>%
    arrange(!!as.name(time_var), !!as.name(id_var)) %>%
    mutate(!!as.name(treat_indicator) := case_when(
      is.na(Treatment_Period) ~ 0,
      !!as.name(time_var) < Treatment_Period ~ 0,
      !!as.name(time_var) >= Treatment_Period ~ 1
    )) %>%
    left_join(
      data_full %>%
        dplyr::select(id_var, time_var, counterfac_var),
      by = c(id_var, time_var)
    )

  return(placebo_df_long)
}

############################################
# Functions for preproccessing the data
###########################################

tsfeature_pc_by_treatment_plot <- function(data_full,
                                           id_var = "entry",
                                           outcome_var = "target",
                                           treat_indicator = "treatperiod_0",
                                           ts_feature_names = NULL) {
  # Generates a scatterplot with the First and Second Principal Components of a number of time series features on each axis,
  # and the individual time series from data_full plotted in that space, colored by treatment status, with centroids.

  # Args
  # data_full: long-form dataframe with both treated and control entries.
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_indicator:column name of binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)
  # ts_feature_names: null for auto features from tsfeatures package (Fix: when formalizing package, pass args for which features you want)

  # Higher level description of data_full:
  # rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  # each row should have a treatment indicator (treat_indicator), a period number (time_var),
  # an individual ID (id_var), and an outcome (outcome_var)
  # for associated with that period-ID combination


  # Output
  # a ggplot object, with the time series by treated and control projected onto the first 2 PCs.

  # identify the treated entries, and assign them a dummy for Treated (in all periods)
  treated_entries <- (data_full %>%
    filter(!!as.name(treat_indicator) > 0) %>%
      distinct(!!as.name(id_var)) %>% 
      pull())
  data_full <- data_full %>%
    mutate(
      Treated =
        as.numeric(!!as.name(id_var) %in% treated_entries)
    )


  list_of_ts <- data_full %>%
    dplyr::select(!!as.name(id_var), !!as.name(outcome_var)) %>%
    split(.[[id_var]]) %>%
    furrr::future_map(~ .[[outcome_var]]) %>%
    furrr::future_map(~ ts(.))

  # Compute the df of Time Series features using tsfeatures package, and compute the Principal Components
  # Fix: need a more systematic way to select out the constant vars that are problematic
  df1_feat <- tsfeatures(list_of_ts)
  df1_feat_pc <- df1_feat %>%
    dplyr::select(-c(frequency, nperiods, seasonal_period)) %>%
    prcomp(scale = TRUE)

  # first, explore overlap in a plot, colored by treatment status
  df1_feat_pc$x %>%
    as_tibble() %>%
    bind_cols(Treated = factor(data_full %>%
      distinct(!!as.name(id_var), .keep_all = T) %>%
      pull(Treated))) %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(col = Treated))

  # Compute the mean of the centroid of PC1 and PC2 for each treatment group
  pc_centroids <- df1_feat_pc$x %>%
    as_tibble() %>%
    bind_cols(Treated = factor(data_full %>%
      distinct(!!as.name(id_var), .keep_all = T) %>%
      pull(Treated))) %>%
    group_by(Treated) %>%
    summarise(
      mean_pc1 = mean(PC1),
      mean_pc2 = mean(PC2)
    )

  # Compute the difference between the two centroids, to get a (naive) sense of overlap
  centroid_dist <- pc_centroids %>%
    pivot_wider(names_from = Treated, values_from = c(mean_pc1, mean_pc2)) %>%
    mutate(
      dist =
        (mean_pc1_0 - mean_pc1_1)^2 + (mean_pc2_0 - mean_pc2_1)^2
    ) %>%
    pull(dist)

  feature_overlap_plot <- df1_feat_pc$x %>%
    as_tibble() %>%
    bind_cols(Treated = factor(data_full %>%
      distinct(entry, .keep_all = T) %>%
      pull(Treated))) %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(col = Treated)) +
    geom_point(aes(x = mean_pc1, y = mean_pc2, col = Treated),
      data = pc_centroids, size = 12, shape = 18
    ) +
    ggtitle("Scatter Plot of First 2 PC by Treatment",
      subtitle = paste(
        "Centroids have L2 dist:",
        round(centroid_dist, 4)
      )
    )

  return(feature_overlap_plot)
}


tsfeature_by_treat_df <- function(data_full,
                                  id_var = "entry",
                                  outcome_var = "target",
                                  treat_indicator = "treatperiod_0",
                                  ts_feature_names = NULL) {
  # Generates a table of the mean of time series features by treatment status

  # Args
  # data_full: long-form dataframe with both treated and control entries.
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_indicator:column name of binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)
  # ts_feature_names: null for auto features from tsfeatures package (Fix: when formalizing package, pass args for which features you want)

  # Higher level description of data_full:
  # rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  # each row should have a treatment indicator (treat_indicator), a period number (time_var),
  # an individual ID (id_var), and an outcome (outcome_var)
  # for associated with that period-ID combination


  # Output
  # a summary statistics table indicating which features are different by treatment status

  # identify the treated entries, and assign them a dummy for Treated (in all periods)
  treated_entries <- (data_full %>%
    filter(!!as.name(treat_indicator) > 0) %>%
    distinct(!!as.name(id_var)) %>%
    pull())
  data_full <- data_full %>%
    mutate(
      Treated =
        as.numeric(!!as.name(id_var) %in% treated_entries)
    )


  list_of_ts <- data_full %>%
    dplyr::select(!!as.name(id_var), !!as.name(outcome_var)) %>%
    split(.[[id_var]]) %>%
    furrr::future_map(~ .[[outcome_var]]) %>%
    furrr::future_map(~ ts(.))
  # Compute the df of Time Series features using tsfeatures package, and t.test by group
  # Fix: need a more systematic way to select out the constant vars that are problematic
  # or try catch those and move on to the rest?
  # Also, clean up the output
  df1_feat <- tsfeatures(list_of_ts) %>%
    as_tibble() %>%
    bind_cols(
      Treated =
        factor(data_full %>%
          distinct(!!as.name(id_var), .keep_all = T) %>%
          pull(Treated))
    ) %>%
    dplyr::select(-c(frequency, nperiods, seasonal_period)) %>%
    pivot_longer(-Treated,
      names_to = "vars",
      values_to = "val"
    ) %>%
    group_by(vars) %>%
    t_test(val ~ Treated) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  return(df1_feat)
}


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
  return(CausalImpact(
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
    filter(!!as.name(treat_indicator) > 0) %>%
    distinct(!!as.name(id_var)) %>%
    pull() %>%
    sort()
  ct_entries <- setdiff(data_full %>%
    distinct(!!as.name(id_var)) %>%
    pull(), tr_entries)
  # Create a dataframe of the subset of control units
  cd <- data_full %>% filter(!!as.name(id_var) %in% ct_entries)
  cd <- merge(cd, data.frame(entry = ct_entries, rank = seq_along(ct_entries)))



  # construct control data matrix
  # For all the control data: Row of matrix (index i) indicates the time, Column (j) indicates observation
  # matrix entry represents the Target
  donor_outcome_matrix <- as.matrix(sparseMatrix(x = cd[[outcome_var]],
                                                 i = cd[[time_var]],
                                                 j = cd$rank))



  # To output Parallelized process, we gather the inputs as a series of lists
  # first, we want the full treated dataframe to manipulate
  td_all <- data_full %>%
    filter(!!as.name(id_var) %in% tr_entries) %>%
    arrange(!!as.name(time_var), !!as.name(id_var)) # NEW, MIGHT CAUSE ERROR

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
  tic()
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
  toc()

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
    inner_join(td_all %>%
      group_by(!!as.name(id_var)) %>%
      mutate(Treatment_Period = length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1) %>%
      ungroup() %>%
      dplyr::select(id_var, time_var, Treatment_Period), by = c(id_var, time_var)) %>%
    arrange(!!as.name(time_var), !!as.name(id_var))


  if (!is.null(counterfac_var)) {
    causalimpact_series_output <- causalimpact_series_output %>%
      inner_join(
        td_all %>%
          dplyr::select(id_var, time_var, counterfac_var),
        by = c(id_var, time_var)
      ) %>%
      mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = (response / !!as.name(counterfac_var)) - 1
      )
  }

  # add a column with relative (pct) effect
  causalimpact_series_output <- causalimpact_series_output %>% mutate(
    pct.effect = (response / point.pred) - 1
  )

  rownames(causalimpact_series_output) <- NULL




  return(causalimpact_series_output)
}





# Potential Improvement: can I use ... as an argument so that any other gsynth arguments (ex, k=5) can be passed? right now, no way to do that for user
estimate_gsynth_series <- function(data_full, id_var = "entry", time_var = "period", treat_indicator = "treatperiod_0", outcome_var = "target", x_in = NULL,
                                   counterfac_var = "counter_factual", se_est = TRUE, num_boots = 1000, inference_type = "parametric",
                                   factor_range = c(0, 5), force_fe = "unit", cross_val = TRUE,
                                   em_flag = FALSE, estimator_type = "ife",
                                   parallel_boot = FALSE) {
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
  # se_est: boolean as to whether uncertainty estimates are provided/estimated
  # num_boots: number of bootstraps run to obtain SE estimates (only applies if se_est=TRUE)
  # inference_type: string -- can be parametric, non-parametric, or jackknife. Parametric is recommended if treatment units are few (approx 40).
  # factor_range: the sequence of unobservable factors to estimate, selected by cross-validation if cross_val=TRUE
  # force_fe: string (unit, time, two-way, none) indicating the type of fixed effects to estimate
  # cross_val: indicate whether to use cross validation to select the optimal number of factors (or the hyperparameter if matrix completion)
  # em_flag: boolean indicating whether EM aglorithm from (Gobillon and Magnac 2016) is to be used for estimating factors
  # estimator_type: string controlling the estimation method -- either Interactive Fixed Effects "ife" or Matrix Completion ("mc")
  # parallel_boot: boolean for whether parallel computing is to be used for bootstrap/se. (FAILS FOR ME....)

  # Higher level description of data_full:
  # rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  # each row should have a treatment indicator (treat_indicator), a period number (time_var),
  # an individual ID (id_var), and an outcome (outcome_var)
  # for associated with that period-ID combination


  # Output

  # estimate the panel SCM
  gsynth_agg_te_all_t <- gsynth(
    Y = outcome_var, D = treat_indicator, data = data_full, index = c(id_var, time_var), X = x_in,
    se = se_est, nboots = num_boots, inference = inference_type, r = factor_range, force = force_fe, CV = cross_val,
    EM = em_flag, estimator = estimator_type, parallel = parallel_boot
  )



  if (se_est & inference_type == "parametric") {
    # gets all indiv effects for all T
    # renames the multi-d array output so that we can more easily pivot it
    gsynth_indiv_te_series <- gsynth_agg_te_all_t$est.ind[, , 1:(gsynth_agg_te_all_t$Ntr)] %>%
      as.data.frame() %>%
      mutate(!!as.name(time_var) := seq_len(nrow(.))) %>%
      clean_names() %>%
      setNames(ifelse(str_count(names(.), "_") > 1,
                      str_replace(names(.), "[_]", ""),
                      names(.)))


    # Pivot the wide dataframe of effects (and CI) per entry with one row per period into a long data set
    # where each period-entry combination has a column for effect, CI.
    # Then, for each period, summarize over all the entries
    gsynth_series_output <- gsynth_indiv_te_series %>%
      pivot_longer(
        -!!as.name(time_var),
        names_to = c(".value", id_var),
        names_sep = "_"
      ) %>%
      mutate(!!as.name(id_var) := as.numeric(!!as.name(id_var))) %>% 
      inner_join(data_full %>%
        group_by(!!as.name(id_var)) %>% 
          mutate(Treatment_Period = length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1) %>%
        ungroup() %>%
        filter(!is.na(Treatment_Period)) %>%
        dplyr::select(id_var, time_var, outcome_var, Treatment_Period), by = c(id_var, time_var)) %>%
      rename(response = outcome_var, point.effect = eff) %>%
      mutate(point.pred = response - point.effect)
  }
  else {
    gsynth_series_output <- (gsynth_agg_te_all_t$Y.tr - gsynth_agg_te_all_t$Y.ct) %>%
      as.data.frame() %>%
      rownames_to_column(var = time_var) %>%
      pivot_longer(
        -!!as.name(time_var),
        names_to = id_var,
        values_to = "point.effect"
      ) %>%
      mutate(
        !!as.name(id_var) := as.numeric(!!as.name(id_var)),
        !!as.name(time_var) := as.numeric(!!as.name(time_var))
      ) %>%
      inner_join(data_full %>%
        group_by(!!as.name(id_var)) %>% 
          mutate(Treatment_Period = length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1) %>%
        ungroup() %>%
        filter(!is.na(Treatment_Period)) %>%
        dplyr::select(id_var, time_var, outcome_var, Treatment_Period), by = c(id_var, time_var)) %>%
      rename(response = outcome_var) %>%
      mutate(point.pred = response - point.effect)
  }

  if (!is.null(counterfac_var)) {
    gsynth_series_output <- gsynth_series_output %>%
      inner_join(
        data_full %>%
          dplyr::select(id_var, time_var, counterfac_var),
        by = c(id_var, time_var)
      ) %>%
      mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = (response / !!as.name(counterfac_var)) - 1
      )
  }

  # add a column with relative (pct) effect
  gsynth_series_output <- gsynth_series_output %>% mutate(
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
  soln <- solve.QP(d_mat, dvec, t(a_t), bvec, meq = meq)
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
  y <- y_df[, 1]
  treatperiod <- y_df[1, 2]

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

  tic("Identified Treat and Control Data")
  # Split the dataset based on whether they are ever treated
  tr_entries <- data_full %>%
    filter(!!as.name(treat_indicator) > 0) %>%
    distinct(!!as.name(id_var)) %>%
    pull() %>%
    sort()
  ct_entries <- setdiff(data_full %>% distinct(!!as.name(id_var)) %>% pull(), tr_entries)


  # create control data frame, with a new id for the sake of ordering observations later
  control_data <- data_full %>% filter(!!as.name(id_var) %in% ct_entries)
  # n0, number of control entries, is just the number of unique entries in cd
  n0 <- control_data %>%
    distinct(!!as.name(id_var)) %>%
    nrow()


  # In the loop, we also want to compute the SCDID Estimates as they require a single observation at a time
  # Because SCDID cannot handle staggered adoption, introduce one treated unit at a time
  # THIS GOES WITHIN Causal Impact For Loop
  treat_data <- data_full %>%
    filter(!!as.name(id_var) %in% tr_entries) %>%
    mutate(new_id = group_indices(., c(!!as.name(id_var))) + n0) %>%
    arrange(!!as.name(time_var), new_id) %>%
    group_by(new_id) %>%
    mutate(Treatment_Period = length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1) %>%
    ungroup()

  toc()

  # create the control matrix once, which is an input to sdid estimator
  control_matrix <- spread(
    control_data %>% dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(outcome_var)),
    !!as.name(time_var), !!as.name(outcome_var)
  ) %>%
    dplyr::select(-!!as.name(id_var)) %>%
    as.matrix() %>%
    t()

  tic("splitting treat data by ID")
  split_treat_data <- treat_data %>% split(.[[id_var]])
  toc()

  # for each treated unit, find when it was treated
  tic("finding list of treat times")
  list_of_treat_times <- split_treat_data %>%
    lapply(., function(x) {
      x %>%
        dplyr::select(Treatment_Period) %>%
        first()
    })
  toc()

  # split the treated units into individual vectors
  tic("creating indiv outcome matrix")
  list_of_treat_data <- split_treat_data %>% lapply(., function(x) {
    x %>%
      dplyr::select(!!as.name(outcome_var)) %>%
      as.matrix()
  })
  toc()

  tic("binidng indiv target to its treat time")
  list_inputs <- Map(cbind, list_of_treat_data, list_of_treat_times)
  toc()

  tic("Mapping to SCDID fun")
  list_of_scdid_series <- furrr::future_map(list_inputs, nsdid_prediction, control_matrix, pre_sdid, post_sdid, nn_sdid, scale_sdid, period_sdid)
  toc()


  # compute the TE by subtracting the matrix of predictions (list_of_scdid_series) from the outcome_var
  # Reformat so that the output is the same as the other functions, namely,
  # each row represents a period-unit combination, with the outcome_var, the prediction, the effect (and the treatment time/indicator)

  df_scdid_series <- list_of_scdid_series %>%
    as.data.frame() %>%
    rownames_to_column(var = time_var) %>%
    mutate(!!as.name(time_var) := as.numeric(!!as.name(time_var))) %>%
    pivot_longer(
      cols = starts_with("X"),
      names_to = "temp_id",
      names_transform = list(temp_id = readr::parse_number),
      values_to = "point.pred"
    ) %>%
    rename(!!as.name(id_var) := temp_id) %>%
    inner_join(
      treat_data %>% 
        dplyr::select(!!as.name(id_var), !!as.name(time_var), !!as.name(outcome_var), Treatment_Period) %>%
        rename(response = outcome_var),
      by = c(id_var, time_var)
    ) %>%
    mutate(point.effect = response - point.pred)


  if (!is.null(counterfac_var)) {
    df_scdid_series <- df_scdid_series %>%
      left_join(
        treat_data %>%
          dplyr::select(id_var, time_var, counterfac_var),
        by = c(id_var, time_var)
      ) %>%
      mutate(
        cf_point.effect = (response - !!as.name(counterfac_var)),
        cf_pct.effect = (response / !!as.name(counterfac_var)) - 1
      )
  }

  # add a column with relative (pct) effect
  df_scdid_series <- df_scdid_series %>% mutate(
    pct.effect = (response / point.pred) - 1
  )

  return(df_scdid_series)
}


# Could add a pure Interactive Fixed Effects Model (interFE) fairly easily






# lots of edits:
# Fix: Include additional weight selection methods? How about adopting the SCDID code and doing double SCDID?
ensemble_placebo_weights <- function(method1_estimated_df, method2_estimated_df, method3_estimated_df, method4_estimated_df,
                                     time_var = "period", id_var = "entry", treat_time_var = "Treatment_Period",
                                     pred_var = "point.pred", counterfac_var = "counter_factual",
                                     pos_coef_constr = F, intercept_allowed = T) {

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

  # Output
  # weights, as a 3x1 (num methods X 1) vector, from an unconstrained linear reg

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

  if (intercept_allowed == T & pos_coef_constr == T) {
    x <- cbind(1, as.matrix(cbind(
      combined_methods_df %>% pull(m1_pred),
      combined_methods_df %>% pull(m2_pred),
      combined_methods_df %>% pull(m3_pred),
      combined_methods_df %>% pull(m4_pred)
    )))

    r_inv <- solve(chol(t(x) %*% x))
    c <- t(cbind(0, rbind(1, diag(4))))

    b <- c(1, rep(0, 4))

    d <- t(combined_methods_df %>% pull(!!as.name(counterfac_var))) %*% x
    nn2 <- sqrt(norm(d, "2"))

    constr_weights_sol <- solve.QP(Dmat = r_inv * nn2, factorized = TRUE, dvec = d / (nn2^2), Amat = c, bvec = b, meq = 1)
    weight_vec <- constr_weights_sol$solution
  }

  if (intercept_allowed == F & pos_coef_constr == T) {
    x <- as.matrix(cbind(
      combined_methods_df %>% pull(m1_pred),
      combined_methods_df %>% pull(m2_pred),
      combined_methods_df %>% pull(m3_pred),
      combined_methods_df %>% pull(m4_pred)
    ))


    r_inv <- solve(chol(t(x) %*% x))
    c <- cbind(rep(1, 4), diag(4))
    b <- c(1, rep(0, 4))

    d <- t(combined_methods_df %>% pull(!!as.name(counterfac_var))) %*% x
    nn2 <- sqrt(norm(d, "2"))

    constr_weights_sol <- solve.QP(Dmat = r_inv * nn2, factorized = TRUE, dvec = d / (nn2^2), Amat = c, bvec = b, meq = 1)
    weight_vec <- constr_weights_sol$solution
  }

  if (intercept_allowed == T & pos_coef_constr == F) {
    lm_formula <- paste(counterfac_var, "~+offset(m1_pred)+I(m2_pred-m1_pred)+I(m3_pred-m1_pred)+I(m4_pred-m1_pred)")
    ols_ensemble_weights <- lm(as.formula(lm_formula), data = post_treat_combined_df)

    coefs_temp <- unname(ols_ensemble_weights$coefficients)[-1]

    # intercept comes first, then rest
    weight_vec <- c(ols_ensemble_weights$coefficients[[1]], 1 - sum(coefs_temp), coefs_temp)
  }
  if (intercept_allowed == F & pos_coef_constr == F) {
    lm_formula <- paste(counterfac_var, "~+offset(m1_pred)+I(m2_pred-m1_pred)+I(m3_pred-m1_pred)+I(m3_pred-m1_pred)+0")
    ols_ensemble_weights <- lm(as.formula(lm_formula), data = post_treat_combined_df)

    coefs_temp <- unname(ols_ensemble_weights$coefficients)
    weight_vec <- c(
      1 - sum(coefs_temp),
      coefs_temp
    )
  }




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



# New fixes: recreating "point effect" and "pct.effect" variables within this code, wnat to allow that as an option
compute_tot_se_jackknife <- function(estimated_series_df, time_var = "period", treat_period_var = "Treatment_Period",
                                     pred_var = "point.pred", outcome_var = "response",
                                     stat_in = "mean", alpha_ci = 0.95,
                                     compute_cf_eff = T, counterfac_var = "counter_factual") {
  # Computes jackknife estimates (from resample package) of the treatment effect on the treated, by period

  # Args
  # estimated_series_df: long-form dataframe with the estimated effects, aka output from one of the methods above.
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_period_var:column name of the period for which the unit's treatment starts (assumes it ends at T max)
  # stat_in: string indicating the statistic which we hope to get a bootstrap mean and se for, either "mean", "aggregate", or "median" typical
  # alpha_ci: number between 0 and 1, indicating the confidence interval desires
  # compute_cf_eff: boolean flag for whether a column of counterfactual effects should be appended to output
  # counterfac_eff: string variable name of the counterfactual effect to estimate, typically "cf_point.effect or "cf_pct.effect" from above

  # Output
  # Tibble containing, by post_treat_period, the sample stat and bootstrapped mean of boot_var, as well as upper and lower bounds
  # that make up the alpha_ci*100% confidence interval

  effect_series_postonly <- estimated_series_df %>%
    mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
    filter(post_period_t >= 0)

  if (compute_cf_eff) {
    cf_tot_df <- estimated_series_df %>%
      mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
      filter(post_period_t >= 0)
    if (tolower(stat_in) == "aggregate") {
      cf_tot_df <- cf_tot_df %>%
        group_by(post_period_t) %>%
        summarise(
          agg_abs_cf_tot = sum(!!as.name(outcome_var) - !!as.name(counterfac_var)),
          agg_pct_cf_tot = (sum(!!as.name(outcome_var)) / sum(!!as.name(counterfac_var))) - 1
        ) %>%
        ungroup()
    }
    if (tolower(stat_in) == "mean") {
      cf_tot_df <- cf_tot_df %>%
        group_by(post_period_t) %>%
        summarise(
          mean_cf_tot = mean(!!as.name(counterfac_var)),
          mean_pct_cf_tot = mean((!!as.name(outcome_var) / !!as.name(counterfac_var)) - 1)
        ) %>%
        ungroup()
    }
    if (tolower(stat_in) == "median") {
      cf_tot_df <- cf_tot_df %>%
        group_by(post_period_t) %>%
        summarise(
          median_abs_cf_tot = median(!!as.name(outcome_var) - !!as.name(counterfac_var)),
          median_pct_cf_tot = median((!!as.name(outcome_var) / !!as.name(counterfac_var)) - 1)
        ) %>%
        ungroup()
    }
  }

  # store the bootstrap metric by post treatment period
  if (tolower(stat_in) == "aggregate") {
    # Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ jackknife(., sum))

    jackknife_comp_pct <- effect_series_postonly %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ jackknife(., function(x) {
        (sum(x[[outcome_var]]) / sum(x[[pred_var]])) - 1
      }))
  }

  if (tolower(stat_in) == "mean") {
    # Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ jackknife(., mean))

    # Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp_pct <- effect_series_postonly %>%
      mutate(pct_eff = (!!as.name(outcome_var) / !!as.name(pred_var)) - 1) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["pct_eff"]]) %>%
      furrr::future_map(~ jackknife(., mean))
  }

  if (tolower(stat_in) == "median") {
    # Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ jackknife(., median))

    # Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp_pct <- effect_series_postonly %>%
      mutate(pct_eff = (!!as.name(outcome_var) / !!as.name(pred_var)) - 1) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["pct_eff"]]) %>%
      furrr::future_map(~ jackknife(., median))
  }


  statname_abs <- paste("jackknife_", stat_in, "_abs_tot", sep = "")
  statname_pct <- paste("jackknife_", stat_in, "_pct_tot", sep = "")
  jackknife_comp_ci_abs <- furrr::future_pmap(list(jackknife_comp_abs, list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = effect_series_postonly %>% 
             distinct(post_period_t) %>%
             pull()) %>%
    mutate(
      !!as.name(statname_abs) := jackknife_comp_abs %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist()
    )
  names(jackknife_comp_ci_abs)[1:2] <- c("jackknife_lb_abs_tot", "jackknife_ub_abs_tot")

  jackknife_comp_ci_pct <- furrr::future_pmap(list(jackknife_comp_pct,
                                            list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = effect_series_postonly %>%
             distinct(post_period_t) %>%
             pull()) %>%
    mutate(
      !!as.name(statname_pct) := jackknife_comp_pct %>%
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist()
    )

  names(jackknife_comp_ci_pct)[1:2] <- c("jackknife_lb_pct_tot", "jackknife_ub_pct_tot")

  jackknife_comp_ci <- jackknife_comp_ci_abs %>% inner_join(jackknife_comp_ci_pct, by = "post_period_t")
  # Manual version of the above code
  # jackknife_mean_bias= jackknife_comp%>% furrr::future_map(~ .$"stats") %>% furrr::future_map(~ .$"Mean") %>%
  #   unlist() %>% tibble("de_biased_mean"=.,  "post_period_t"=effect_series_postonly %>% distinct(post_period_t) %>% pull())
  #
  # jackknife_se=jackknife_comp%>% furrr::future_map(~ .$"stats") %>% furrr::future_map(~ .$"SE") %>%
  #   unlist() %>% tibble("se"=.,"post_period_t"=effect_series_postonly %>% distinct(post_period_t) %>% pull())


  stat_tot_sample <- effect_series_postonly %>%
    group_by(post_period_t) %>%
    summarise(treated_n = n()) %>%
    ungroup()


  if (compute_cf_eff) {
    jackknife_ci_by_postperiod <- stat_tot_sample %>%
      left_join(
        jackknife_comp_ci,
        by = "post_period_t"
      ) %>%
      left_join(
        cf_tot_df,
        by = "post_period_t"
      )
  }
  if (!compute_cf_eff) {
    jackknife_ci_by_postperiod <- stat_tot_sample %>% 
      left_join(
      jackknife_comp_ci,
      by = "post_period_t"
    )
  }


  return(jackknife_ci_by_postperiod)
}




# Uncertainty Bootstrap
compute_ci_bounds_bootstrap <- function(estimated_series_df, time_var = "period", treat_period_var = "Treatment_Period",
                                        pred_var = "point.pred", outcome_var = "response",
                                        stat_in = "mean", alpha_ci = 0.95,
                                        compute_cf_eff = T, counterfac_var = "counter_factual",
                                        nboots = 10000) {
  # Computes bootstrap estimates (from resample package) of the treatment effect on the treated, by period

  # Args
  # estimated_series_df: long-form dataframe with the estimated effects, aka output from one of the methods above.
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_period_var:column name of the period for which the unit's treatment starts (assumes it ends at T max)
  # boot_var: string indicating the effect we wish to bootstrap, either "point.effect" or "pct.effect" works
  # stat_in: string indicating the statistic which we hope to get a bootstrap mean and se for, either "mean", "aggregate", "median" typical
  # alpha_ci: number between 0 and 1, indicating the confidence interval desires
  # compute_cf_eff: boolean flag for whether a column of counterfactual effects should be appended to output
  # counterfac_eff: string variable name of the counterfactual effect to estimate, typically "cf_point.effect or "cf_pct.effect" from above
  # nboots, number of bootstrap samples

  # Output
  # Tibble containing, by post_treat_period, the sample stat and bootstrapped mean of boot_var, as well as upper and lower bounds
  # that make up the alpha_ci*100% confidence interval

  effect_series_postonly <- estimated_series_df %>%
    mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
    filter(post_period_t >= 0)

  if (compute_cf_eff) {
    cf_tot_df <- estimated_series_df %>%
      mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
      filter(post_period_t >= 0)
    if (tolower(stat_in) == "aggregate") {
      cf_tot_df <- cf_tot_df %>%
        group_by(post_period_t) %>%
        summarise(
          agg_abs_cf_tot = sum(!!as.name(outcome_var) - !!as.name(counterfac_var)),
          agg_pct_cf_tot = (sum(!!as.name(outcome_var)) / sum(!!as.name(counterfac_var))) - 1
        ) %>%
        ungroup()
    }
    if (tolower(stat_in) == "mean") {
      cf_tot_df <- cf_tot_df %>%
        group_by(post_period_t) %>%
        summarise(
          mean_cf_tot = mean(!!as.name(counterfac_var)),
          mean_pct_cf_tot = mean((!!as.name(outcome_var) / !!as.name(counterfac_var)) - 1)
        ) %>%
        ungroup()
    }
    if (tolower(stat_in) == "median") {
      cf_tot_df <- cf_tot_df %>%
        group_by(post_period_t) %>%
        summarise(
          median_abs_cf_tot = median(!!as.name(outcome_var) - !!as.name(counterfac_var)),
          median_pct_cf_tot = median((!!as.name(outcome_var) / !!as.name(counterfac_var)) - 1)
        ) %>%
        ungroup()
    }
  }

  # is the double future-map here the slow thing? How does the speed compare to lapply?
  if (tolower(stat_in) == "aggregate") {
    # Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ bootstrap(., mean, R = nboots))

    bootstrap_comp_pct <- effect_series_postonly %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ bootstrap(., function(x) {
        (sum(x[[outcome_var]]) / sum(x[[pred_var]])) - 1
      }, R = nboots))
  }
  if (tolower(stat_in) == "mean") {
    # Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ bootstrap(., mean, R = nboots))

    # Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp_pct <- effect_series_postonly %>%
      mutate(pct_eff = (!!as.name(outcome_var) / !!as.name(pred_var)) - 1) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["pct_eff"]]) %>%
      furrr::future_map(~ bootstrap(., mean, R = nboots))
  }

  if (tolower(stat_in) == "median") {
    # Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp_abs <- effect_series_postonly %>%
      mutate(abs_eff = !!as.name(outcome_var) - !!as.name(pred_var)) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["abs_eff"]]) %>%
      furrr::future_map(~ bootstrap(., median, R = nboots))

    # Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp_pct <- effect_series_postonly %>%
      mutate(pct_eff = (!!as.name(outcome_var) / !!as.name(pred_var)) - 1) %>%
      split(.[["post_period_t"]]) %>%
      furrr::future_map(~ .[["pct_eff"]]) %>%
      furrr::future_map(~ bootstrap(., median, R = nboots))
  }

  statname_abs <- paste("bootstrap_", stat_in, "_abs_tot", sep = "")
  statname_pct <- paste("bootstrap_", stat_in, "_pct_tot", sep = "")
  bootstrap_comp_ci_abs <- furrr::future_pmap(list(bootstrap_comp_abs, list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = effect_series_postonly %>%
             distinct(post_period_t) %>% 
             pull()) %>%
    mutate(
      !!as.name(statname_abs) := bootstrap_comp_abs %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist()
    )
  names(bootstrap_comp_ci_abs)[1:2] <- c("bootstrap_lb_abs_tot", "bootstrap_ub_abs_tot")

  bootstrap_comp_ci_pct <- furrr::future_pmap(list(bootstrap_comp_pct, list(probs = c(0.5 - alpha_ci / 2, 0.5 + alpha_ci / 2))), CI.t) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(post_period_t = effect_series_postonly %>%
             distinct(post_period_t) %>%
             pull()) %>%
    mutate(
      !!as.name(statname_pct) := bootstrap_comp_pct %>% 
        furrr::future_map(~ .[["stats"]]) %>% 
        furrr::future_map(~ .[["Mean"]]) %>%
        unlist()
    )
  names(bootstrap_comp_ci_pct)[1:2] <- c("bootstrap_lb_pct_tot", "bootstrap_ub_pct_tot")

  bootstrap_comp_ci <- bootstrap_comp_ci_abs %>% inner_join(bootstrap_comp_ci_pct, by = "post_period_t")


  stat_tot_sample <- effect_series_postonly %>%
    group_by(post_period_t) %>%
    summarise(treated_n = n()) %>%
    ungroup()



  if (compute_cf_eff) {
    bootstrap_ci_by_postperiod <- stat_tot_sample %>%
      left_join(
        bootstrap_comp_ci,
        by = "post_period_t"
      ) %>%
      left_join(
        cf_tot_df,
        by = "post_period_t"
      )
  }
  if (!compute_cf_eff) {
    bootstrap_ci_by_postperiod <- stat_tot_sample %>% left_join(
      bootstrap_comp_ci,
      by = "post_period_t"
    )
  }



  return(bootstrap_ci_by_postperiod)
}




#################################################################
# Plots
#################################################################



# print PDF Likely does not work with future_map???
create_gap_ci_plot <- function(bootstrapped_effects_df, time_var = "post_period_t", effect_var = "jackknife_mean_pct_tot",
                               upper_ci = "jackknife_ub_pct_tot", lower_ci = "jackknife_lb_pct_tot",
                               cf_plot = T, cf_var = "mean_pct_cf_tot",
                               print_to_pdf = NULL, plot_title = NULL, plot_x_lab = NULL, plot_y_lab = NULL) {

  # Computes bootstrap estimates (from resample package) of the treatment effect on the treated, by period

  # Args
  # bootstrapped_effects_df: data on the bootstrapped (and potentially counterfactual) effects for each post treat period,
  # from compute_ci_bounds_bootstrap or jackknife version
  # time_var:column name of numeric period number indicating the post treatment time period, in increasing order (eg 0 is the first time)
  # effect_var: string indicating the effect we wish to plot, either "jackknife_mean_tot" or "jackknife_median_tot" works, depending on input
  # upper_ci: string indicating thename of the upper CI variable in bootstrapped_effects_df
  # lower_ci: string indicating thename of the lower CI variable in bootstrapped_effects_df
  # cf_plot: boolean flag for whether a counterfactual effect should be plotted
  # cf_var: string variable name of the counterfactual effect to plot, typically "mean_cf_tot or "median_cf_tot" from bootstrapped_effects_df
  # print_to_pdf: string, file path where the pdfs should be printed.
  # plot_title: title of the plot to be printed
  # plot_x_lab: title of the x axis label
  # plot_y_lab: title of the y axis label

  # Output
  # gap plot, potentially with counterfactual predictions

  if (!is.null(print_to_pdf)) {
    pdf(print_to_pdf)
  }

  if (cf_plot) {
    plot_out <- bootstrapped_effects_df %>%
      ggplot(aes(x = !!as.name(time_var), y = !!as.name(effect_var), color = "Estimate")) +
      geom_line() +
      geom_ribbon(aes(ymin = !!as.name(lower_ci), ymax = !!as.name(upper_ci)), alpha = 0.3, color = NA) +
      geom_line(aes(y = !!as.name(cf_var), color = "True")) +
      ggtitle(ifelse(is.null(plot_title), paste("Gap Plot of", effect_var), plot_title)) +
      labs(
        x = ifelse(is.null(plot_x_lab), time_var, plot_x_lab),
        y = ifelse(is.null(plot_y_lab), effect_var, plot_y_lab)
      ) +
      scale_colour_manual(
        name = "",
        values = c("Estimate" = "black", "True" = "red"), labels = c("Estimate", "True")
      )
  }

  if (!cf_plot) {
    plot_out <- bootstrapped_effects_df %>%
      ggplot(aes(x = !!as.name(time_var), y = !!as.name(effect_var))) +
      geom_line() +
      geom_ribbon(aes(ymin = !!as.name(lower_ci), ymax = !!as.name(upper_ci)), alpha = 0.3) +
      ggtitle(ifelse(is.null(plot_title), paste("Gap Plot of", effect_var), plot_title)) +
      labs(
        x = ifelse(is.null(plot_x_lab), time_var, plot_x_lab),
        y = ifelse(is.null(plot_y_lab), effect_var, plot_y_lab)
      )
  }


  if (!is.null(print_to_pdf)) {
    print(plot_out)
    dev.off()
  }
  return(plot_out)
}



#################################################################
# Metrics
#################################################################

compute_avg_metric_per_t <- function(estimated_series_df, time_var = "period", outcome_var = "response", prediction_var = "point.pred",
                                     counterfac_var = "counter_factual", treat_period_var = "Treatment_Period",
                                     metric_str = "mae", pct_eff_flag = F) {

  # Create the gap plot for each of the post treatment periods

  # Args
  # estimated_series_df: long-form dataframe with predicted counterfactual, actual counterfactual, period, and Treatment Period by ID.
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # outcome_var: the true y var for the time series
  # prediction_var: predicted counterfactual by time and id
  # counterfac_var: true counterfactual by time and id
  # treat_period_var: str name of the column which indicates the treatment period for each id_var (at each time_var)
  # metric_str: string indicating desired metric. can be mae, mse, or both
  # pct_eff_flag: binary flag for whether the percent error should be computed

  # Output
  # Dataframe containing number of rows equal to the longest post-treat period
  # and for each post treat period (from 0 -- time of treat, to max), an average of the metric for all observations
  # that experienced that particular post treat period (eg treated in t=5, total T=15 means 10 post treat periods)

  # absolute error computation
  if (pct_eff_flag == FALSE) {
    if (tolower(metric_str) == "mse") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_error_sq = (!!as.name(prediction_var) - !!as.name(counterfac_var))^2) %>%
        group_by(post_period_t) %>%
        summarise(MSE_TE_byT = mean(indiv_error_sq))
    }

    if (tolower(metric_str) == "mae") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error = abs(!!as.name(prediction_var) - !!as.name(counterfac_var))) %>%
        group_by(post_period_t) %>%
        summarise(MAE_TE_byT = mean(indiv_t_error))
    }


    if (tolower(metric_str) == "both") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(
          indiv_t_error = abs(!!as.name(prediction_var) - !!as.name(counterfac_var)),
          indiv_error_sq = (!!as.name(prediction_var) - !!as.name(counterfac_var))^2
        ) %>%
        group_by(post_period_t) %>%
        summarise(
          MAE_TE_byT = mean(indiv_t_error),
          MSE_TE_byT = mean(indiv_error_sq)
        )
    }
  }

  # percent error computation
  if (pct_eff_flag) {
    if (tolower(metric_str) == "mse") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueT- EstTE=(Target/counterFac)-1-(Target/Pred-1)
        mutate(indiv_error_sq = (!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))^2) %>%
        group_by(post_period_t) %>%
        summarise(MSE_TE_byT = mean(indiv_error_sq))
    }

    if (tolower(metric_str) == "mae") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error = abs(!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))) %>%
        group_by(post_period_t) %>%
        summarise(MAE_TE_byT = mean(indiv_t_error))
    }


    if (tolower(metric_str) == "both") {
      avg_metric_series <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(
          indiv_t_error = abs(!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var))),
          indiv_error_sq = (!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))^2
        ) %>%
        group_by(post_period_t) %>%
        summarise(
          MAE_TE_byT = mean(indiv_t_error),
          MSE_TE_byT = mean(indiv_error_sq)
        )
    }
  }




  return(avg_metric_series)
}



compute_avg_metric <- function(estimated_series_df, time_var = "period", outcome_var = "response", prediction_var = "point.pred",
                               counterfac_var = "counter_factual", treat_period_var = "Treatment_Period",
                               metric_str = "mae", pct_eff_flag = F) {
  # Estimates the relevant metric for the predicted Treatment Effect by comparing the predicted and true counterfactuals

  # Args
  # estimated_series_df: long-form dataframe with predicted counterfactual, actual counterfactual, period, and Treatment Period by ID.
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # outcome_var: the true y var for the time series
  # prediction_var: predicted counterfactual by time and id
  # counterfac_var: true counterfactual by time and id
  # treat_period_var: str name of the column which indicates the treatment period for each id_var (at each time_var)
  # metric_str: string indicating desired metric. can be mae, mse, or both
  # pct_eff_flag: binary flag for whether the percent error should be computed

  # Output
  # tibble containing an average of the metric for all observations, over all T

  # absolute error computation
  if (pct_eff_flag == FALSE) {
    if (tolower(metric_str) == "mse") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_error_sq = (!!as.name(prediction_var) - !!as.name(counterfac_var))^2) %>%
        summarise(MSE_TE = mean(indiv_error_sq))
    }

    if (tolower(metric_str) == "mae") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error = abs(!!as.name(prediction_var) - !!as.name(counterfac_var))) %>%
        summarise(MAE_TE = mean(indiv_t_error))
    }


    if (tolower(metric_str) == "both") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(
          indiv_t_error = abs(!!as.name(prediction_var) - !!as.name(counterfac_var)),
          indiv_error_sq = (!!as.name(prediction_var) - !!as.name(counterfac_var))^2
        ) %>%
        summarise(
          MAE_TE = mean(indiv_t_error),
          MSE_TE = mean(indiv_error_sq)
        )
    }
  }

  # percent error computation
  if (pct_eff_flag) {
    if (tolower(metric_str) == "mse") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueT- EstTE=(Target/counterFac)-1-(Target/Pred-1)
        mutate(indiv_error_sq = (!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))^2) %>%
        summarise(MSE_TE = mean(indiv_error_sq))
    }

    if (tolower(metric_str) == "mae") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error = abs(!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))) %>%
        summarise(MAE_TE = mean(indiv_t_error))
    }


    if (tolower(metric_str) == "both") {
      avg_metric <- estimated_series_df %>%
        mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
        filter(post_period_t >= 0) %>% # note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(
          indiv_t_error = abs(!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var))),
          indiv_error_sq = (!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))^2
        ) %>%
        summarise(
          MAE_TE = mean(indiv_t_error),
          MSE_TE = mean(indiv_error_sq)
        )
    }
  }




  return(avg_metric)
}


compute_tot_coverage <- function(bootstrapped_series_df, time_var = "post_period_t", ub_var = "jackknife_ub_pct_tot",
                                 lb_var = "jackknife_lb_pct_tot", counterfac_var = "mean_pct_cf_tot") {
  coverage_out <- bootstrapped_series_df %>%
    mutate(
      contains_truth = as.numeric(!!as.name(counterfac_var) <= !!as.name(ub_var) &
        !!as.name(counterfac_var) >= !!as.name(lb_var))
    ) %>%
    pull(contains_truth) %>%
    mean()

  return(coverage_out)
}




plot_tot_bias_per_t <- function(estimated_series_df, time_var = "period", id_var = "entry", outcome_var = "response", prediction_var = "point.pred",
                                counterfac_var = "counter_factual", treat_period_var = "Treatment_Period",
                                max_post_t = 12, pct_eff_flag = F, plot_title = NULL) {

  # Create the gap plot for each of the post treatment periods

  # Args
  # estimated_series_df: long-form dataframe with predicted counterfactual, actual counterfactual, period, and Treatment Period by ID.
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # id_var: column name of the unit id
  # outcome_var: the true y var for the time series
  # prediction_var: predicted counterfactual by time and id
  # counterfac_var: true counterfactual by time and id
  # treat_period_var: str name of the column which indicates the treatment period for each id_var (at each time_var)
  # max_post_t: numeric variable, the number of time periods over which we are averaging per unit
  # pct_eff_flag: binary flag for whether the percent error should be computed
  # plot_title: string title of the plot

  # Output
  # Dataframe containing number of rows equal to the longest post-treat period
  # and for each post treat period (from 0 -- time of treat, to max), an average of the metric for all observations
  # that experienced that particular post treat period (eg treated in t=5, total T=15 means 10 post treat periods)

  # absolute error computation
  if (pct_eff_flag == FALSE) {
    avg_bias_series <- estimated_series_df %>%
      mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
      filter(post_period_t >= 0, post_period_t <= max_post_t) %>%
      group_by(!!as.name(id_var)) %>%
      summarise(bias_per_id = mean(!!as.name(prediction_var) - !!as.name(counterfac_var))) %>%
      ggplot(aes(x = bias_per_id)) +
      geom_density(fill = "blue", alpha = 0.4) +
      theme_bw() +
      theme(
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      scale_x_continuous(name = "ATT by ID") +
      ggtitle(ifelse(is.null(plot_title), paste("ToT Bias by Unit, averaged over", max_post_t, "post-treat Periods"), plot_title))
  }

  # percent error computation
  if (pct_eff_flag) {
    avg_bias_series <- estimated_series_df %>%
      mutate(post_period_t = !!as.name(time_var) - !!as.name(treat_period_var)) %>%
      filter(post_period_t >= 0) %>%
      group_by(!!as.name(id_var)) %>%
      summarise(bias_per_id = mean(!!as.name(outcome_var) / !!as.name(counterfac_var) - (!!as.name(outcome_var) / !!as.name(prediction_var)))) %>%
      ggplot(aes(x = bias_per_id)) +
      geom_density(fill = "blue", alpha = 0.4) +
      theme_bw() +
      theme(
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      scale_x_continuous(name = "ATT by ID") +
      ggtitle(ifelse(is.null(plot_title), paste("ToT Pct Bias by Unit, averaged over", max_post_t, "post-treat Periods"), plot_title))
  }




  return(avg_bias_series)
}
