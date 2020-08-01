pacman::p_load(dplyr, ggplot2, quadprog,tsfeatures,tibble, tidyr,rstatix)


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
                        dplyr::filter(!!as.name(treat_indicator) > 0) %>%
                        dplyr::distinct(!!as.name(id_var)) %>% 
                        dplyr::pull())
  data_full <- data_full %>%
    dplyr::mutate(Treated = as.numeric(!!as.name(id_var) %in% treated_entries))
  
  
  list_of_ts <- data_full %>%
    dplyr::select(!!as.name(id_var), !!as.name(outcome_var)) %>%
    split(.[[id_var]]) %>%
    furrr::future_map(~ .[[outcome_var]]) %>%
    furrr::future_map(~ stats::ts(.))
  
  # Compute the df of Time Series features using tsfeatures package, and compute the Principal Components
  # Fix: need a more systematic way to select out the constant vars that are problematic
  df1_feat <- tsfeatures::tsfeatures(list_of_ts)
  df1_feat_pc <- df1_feat %>%
    dplyr::select(-c(frequency, nperiods, seasonal_period,
                     e_acf10, x_acf10, diff2_acf10, diff1_acf10)) %>%
    stats::prcomp(scale = TRUE)
  
  # first, explore overlap in a plot, colored by treatment status
  df1_feat_pc$x %>%
    tibble::as_tibble() %>%
    dplyr::bind_cols(
      Treated = factor(data_full %>%
                         dplyr::distinct(!!as.name(id_var), .keep_all = T) %>%
                         dplyr::pull(Treated))) %>%
    ggplot2::ggplot(aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(aes(col = Treated))
  
  # Compute the mean of the centroid of PC1 and PC2 for each treatment group
  pc_centroids <- df1_feat_pc$x %>%
    tibble::as_tibble() %>%
    dplyr::bind_cols(
      Treated = factor(data_full %>%
                         dplyr::distinct(!!as.name(id_var), .keep_all = T) %>%
                         dplyr::pull(Treated))) %>%
    dplyr::group_by(Treated) %>%
    dplyr::summarise(
      mean_pc1 = mean(PC1),
      mean_pc2 = mean(PC2)
    )
  
  # Compute the difference between the two centroids, to get a (naive) sense of overlap
  centroid_dist <- pc_centroids %>%
    tidyr::pivot_wider(names_from = Treated, 
                       values_from = c(mean_pc1, mean_pc2)) %>%
    dplyr::mutate(
      dist =
        (mean_pc1_0 - mean_pc1_1)^2 + (mean_pc2_0 - mean_pc2_1)^2
    ) %>%
    dplyr::pull(dist)
  
  feature_overlap_plot <- df1_feat_pc$x %>%
    tibble::as_tibble() %>%
    dplyr::bind_cols(
      Treated = factor(data_full %>%
                         dplyr::distinct(entry, .keep_all = T) %>%
                         dplyr::pull(Treated))) %>%
    ggplot2::ggplot(aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(aes(col = Treated)) +
    ggplot2::geom_point(aes(x = mean_pc1, y = mean_pc2, col = Treated),
               data = pc_centroids, size = 12, shape = 18 ) +
    ggplot2::ggtitle("Scatter Plot of First 2 PC by Treatment",
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
                        dplyr::filter(!!as.name(treat_indicator) > 0) %>%
                        dplyr::distinct(!!as.name(id_var)) %>%
                        dplyr::pull())
  data_full <- data_full %>%
    dplyr::mutate( Treated = as.numeric(!!as.name(id_var) %in% treated_entries))
  
  
  list_of_ts <- data_full %>%
    dplyr::select(!!as.name(id_var), !!as.name(outcome_var)) %>%
    split(.[[id_var]]) %>%
    furrr::future_map(~ .[[outcome_var]]) %>%
    furrr::future_map(~ ts(.))
  # Compute the df of Time Series features using tsfeatures package, and t.test by group
  # Fix: need a more systematic way to select out the constant vars that are problematic
  # or try catch those and move on to the rest?
  # Also, clean up the output
  df1_feat <- tsfeatures::tsfeatures(list_of_ts) %>%
    tibble::as_tibble() %>%
    dplyr::bind_cols(
      Treated =
        factor(data_full %>%
                 dplyr::distinct(!!as.name(id_var), .keep_all = T) %>%
                 dplyr::pull(Treated))
    ) %>%
    dplyr::select(-c(frequency, nperiods, seasonal_period,
                     e_acf10, x_acf10, diff2_acf10, diff1_acf10)) %>%
    tidyr::pivot_longer(-Treated,
                 names_to = "vars",
                 values_to = "val"
    ) %>%
    dplyr::group_by(vars) %>%
    rstatix::t_test(val ~ Treated) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance()
  return(df1_feat)
}
