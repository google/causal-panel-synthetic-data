library(tidyverse)
library(dplyr)
library(furrr)
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