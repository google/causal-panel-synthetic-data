pacman::p_load(dplyr, tidyr, stats, tibble)
############################################
# Functions to create a placebo dataset,
# assigning treatment based on matching.
###########################################

#' Helper function for .MatchingWithoutReplacement.
#' Measures the euclidean distance between ts_tomatch
#' and each of the rows in the tibble ts_rest.
#' 
#' @param ts_tomatch Treated time series tibble, rows are ID and columns Time.
#' @param ts_rest Donor time series tibble, rows identify ID and columns Time.
#' 
#' @return Vector of indices each indicating the row in ts_rest that is best 
#'      matched (min L2 norm distance) to the row in ts_tomatch.
#'      
#' @noRd
.NearestTSEuclidean <- function(ts_tomatch, ts_rest) {
  nearest_ts <- apply(ts_tomatch, 1, function(ts_tomatch) {
    which.min(
      apply(
        ts_rest, 1, function(ts_rest, ts_tomatch) {
          stats::dist(rbind(ts_rest, ts_tomatch))
        },
        ts_tomatch
      )
    )
  })
  return(nearest_ts)
}

#' Finds the nearest match for each treated unit (treated_block)
#' among the donor pool (control_block) without replacement by calling helper
#' .NearestTSEuclidean.
#'
#' @param treated_block Treated time series tibble, dimension ID by Time.
#' @param control_block Donor time series tibble, dimension ID by Time.
#' 
#' @return A tibble containing a column for the placebo-treated unit ID,
#'     the treated unit it was the nearest match to, and the treatment time of 
#'     the true treated unit.
#' @noRd
.MatchingWithoutReplacement <- function(treated_block,
                                         control_block,
                                         id_var,
                                         treat_period) {
  # Initialize an empty vector for the donor IDs that match.
  already_matched <- c()
  # Randomize the order of matching for this greedy approach.
  treated_block <- treated_block %>% dplyr::sample_n(dplyr::n())
  # Store the time of treatment and treated ID for the true treated units.
  placebo_treat_period <- treated_block %>% 
    dplyr::pull(tidyselect::all_of(treat_period))
  treatment_unit <- treated_block %>% dplyr::pull(!!as.name(id_var))
  for (i in seq_len(nrow(treated_block))) {
    # Find nearest euclidean match among unmatched controls for 
    # each treated observation.
    temp_match <- .NearestTSEuclidean(
      treated_block %>%
        dplyr::slice(i) %>%
        dplyr::select(-tidyselect::all_of(c(id_var,treat_period))),
      control_block %>%
        dplyr::filter(!!as.name(id_var) %in% setdiff(!!as.name(id_var), 
                                              already_matched)) %>% 
        dplyr::select(-tidyselect::all_of(id_var)))
    # Update the vector of already matched donors.
    already_matched[i] <- control_block %>%
      dplyr::filter(!!as.name(id_var) %in% 
               setdiff(!!as.name(id_var), already_matched)) %>%
      dplyr::slice(temp_match) %>%
      dplyr::pull(!!as.name(id_var))
  }
  # Store the resulting vectors in a tibble for output.
  df_toreturn <- tibble::tibble(temp_id = already_matched, 
                        temp_treattime = placebo_treat_period, 
                        Treatment_Unit = treatment_unit) %>%
    dplyr::rename(!!as.name(id_var) := temp_id,
                  !!as.name(treat_period) := temp_treattime)
  
  return(df_toreturn)
}
#' Generates a placebo tibble, using matching methods to select
#' placebo-treated entries as those most similar to truly-treated.
#'
#' @param data_full Long-form tibble with both treated and control entries.
#'    rows represent period-entry combinations  (e.g. N (total num of entry) 
#'    rows for period t). Each row should have a treatment indicator 
#'    (treat_indicator),  a period number (time_var), an individual ID 
#'    (id_var), and an outcome (outcome_var) for associated with that 
#'    period-ID combination.
#' @param id_var Column name of numeric, unique ID representing the entry.
#' @param time_var Column name of numeric period number (the time period),
#'    in increasing order (e.g. 0 is the first time, 120 is the last).
#' @param treat_indicator Column name of binary (0, 1) indicator for whether 
#'     the unit in that time period was treated. Should never decrease.
#'     
#' @return A tibble of the same format as data_full, 
#'    entirely consisting of donor units, some of which are placebo-treated 
#'    (based on matching).
# TODO(alexdkellogg):  Revisit whether pivoting will cause issues in the 
#    presence of several time varying covariates.
# TODO(alexdkellogg): Account for multiple treatments, treatment end dates.
# TODO(alexdkellogg): Explore alternative methods of matching TS, e.g. estimated
#    probability of having been treated given lagged outcomes (up to treatment).
CreatePlaceboData <- function(data_full, id_var = "entry",
                              time_var = "period", 
                              treat_indicator = "treatperiod_0",
                              outcome_var = "target", 
                              counterfac_var = "counter_factual") {
  # Split the dataset based on whether the unit is ever treated.
  tr_entries <- data_full %>%
    dplyr::filter(!!as.name(treat_indicator) > 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>%
    dplyr::pull()
  
  ct_entries <- setdiff(data_full %>%
                          dplyr::distinct(!!as.name(id_var)) %>%
                          dplyr::pull(), tr_entries)
  # Create a dataframe of the subset of control units.
  cd <- data_full %>% dplyr::filter(!!as.name(id_var) %in% ct_entries)
  cd <- merge(
    cd,
    data.frame(
      entry = ct_entries,
      rank = seq_along(ct_entries)
    )
  )
  # Pivot the long cd dataframe to wide. Each row will represent an id_var, 
  # with columns for the outcome_var at each time period.
  cd_for_match <- tidyr::pivot_wider(
    data = cd %>%
      dplyr::arrange(!!as.name(time_var), !!as.name(id_var)),
    names_from = !!as.name(time_var),
    id_cols = c(!!as.name(id_var)),
    values_from = c(!!as.name(outcome_var))
  )
  # Store treated data in a tibble.
  treated_to_match <- data_full %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries)
  # Create a variable indicating the time treatment is assigned by unit.
  treated_to_match <- treated_to_match %>%
    dplyr::group_by(!!as.name(id_var)) %>%
    dplyr::mutate(
      Treatment_Period =
        length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1
    ) %>%
    dplyr::ungroup()
  # Pivot the data wide, so each row has all the time series data for 
  # a given unit.
  data_wide_m <- tidyr::pivot_wider(
    data = treated_to_match,
    names_from = !!as.name(time_var),
    id_cols = c(!!as.name(id_var), Treatment_Period),
    values_from = c(!!as.name(outcome_var))
  )
  
  # Determine control units to assign to placebo treatment group via matching.
  matched_placebo_df_temp <- .MatchingWithoutReplacement(data_wide_m, 
                                                          cd_for_match, 
                                                          id_var,
                                                          "Treatment_Period")
  # We now have the set of control units that form the placebo set, along with 
  # their placebo Treat Period and their corresponding Treat Unit.
  # Assign treatment to our donor only tibble based on matches.
  placebo_df_wide <- cd_for_match %>%
    dplyr::left_join(matched_placebo_df_temp, by = id_var)
  
  # Pivot_long the tibble and recreate the treat_indicator indicator
  # equal to 1 if unit i in time t is treated, 
  # so that the output matches the input.
  non_time_vars <- c(id_var, "Treatment_Period", "Treatment_Unit")
  placebo_df_long <- placebo_df_wide %>%
    tidyr::pivot_longer(-all_of(non_time_vars),
                 names_to = time_var, values_to = outcome_var
    ) %>%
    dplyr::mutate(!!as.name(time_var) := as.numeric(!!as.name(time_var))) %>%
    dplyr::arrange(!!as.name(time_var), !!as.name(id_var)) %>%
    dplyr::mutate(!!as.name(treat_indicator) := dplyr::case_when(
      is.na(Treatment_Period) ~ 0,
      !!as.name(time_var) < Treatment_Period ~ 0,
      !!as.name(time_var) >= Treatment_Period ~ 1
    )) %>%
    dplyr::left_join(
      data_full %>%
        dplyr::select(tidyselect::all_of(c(id_var, time_var, counterfac_var))),
      by = c(id_var, time_var)
    )
  # If no counterfactual variable, create one in the placebo from the outcome.
  if(is.null(counterfac_var)) {
    placebo_df_long <- placebo_df_long %>%
      dplyr::mutate("counter_factual"=!!as.name(outcome_var))
  }
  return(placebo_df_long)
}
