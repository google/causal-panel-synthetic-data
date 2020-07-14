library(tidyverse)
library(dplyr)
library(furrr)
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



#TODO(alexdkellogg): Expand to matching with and storing covariates
#                   (pivot wide/long, specifically).
#TODO(alexdkellogg): account for multiple treatments, treatment end dates
create_placebo_df <- function(data_full, id_var = "entry",
                              time_var = "period", 
                              treat_indicator = "treatperiod_0",
                              outcome_var = "target", 
                              counterfac_var = "counter_factual") {
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
  
  # Determine control units to assign to placebo treatment group via matching.
  matched_placebo_df_temp <- matching_without_replacement(data_wide_m, 
                                                          cd_for_match, 
                                                          id_var,
                                                          "Treatment_Period")
  
  # We now have the set of control units that form the placebo set, along with 
  #their placebo Treat Period and their corresponding Treat Unit
  # Assign treatment to our donor only tibble based on matches.
  placebo_df_wide <- cd_for_match %>%
    dplyr::left_join(matched_placebo_df_temp, by = id_var)
  
  # Pivot_long the tibble and recreate the treat_indicator indicator
  # equal to 1 if unit i in time t is treated, 
  # so that the output matches the input.
  # If there are covariates or other variables which we do not want to be 
  # pivoting, append them to non_time_vars vectors
  # may have to add values_to as well if we have multiple time varying 
  # variables (added as TODO)
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
  
  return(placebo_df_long)
}
