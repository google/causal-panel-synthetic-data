pacman::p_load(dplyr, ggplot2, janitor,tsfeatures,tibble, tidyr,rstatix)
############################################
# Functions for visualizing and analyzing time series (TS) features.
###########################################

#' Visualize time series characteristics by treatment status.
#'
#' @param data_full Raw data in long form, one row per unit time combination.
#' @param id_var String name of the unit identifier.
#' @param outcome_var String name of the time series outcome identifier.
#' @param treat_indicator String name of the treatment status indicator.
#' @param type String name for the desired plot type.
#' @param ... Additional arguments to be passed to tsfeatures function from the 
#'    tsfeatures package. Examples include the name of the features and whether
#'    to trim the outliers.
#'
#' @return A scatterplot of the first two principal components of the time 
#'    series features, by treatment status.
TSFeaturesPlot <- function(data_full,
                                    id_var = "entry",
                                    outcome_var = "target",
                                    treat_indicator = "treatperiod_0",
                                    type=c("rug", "centroid", "contour"),
                                    ...) {
  type <- match.arg(type)
  # Identify treated observations.
  data_full <- .IdentifyTreated(data_full=data_full, id_var=id_var,
                                treat_indicator=treat_indicator)
  
  # Compute the features for each TS, and the principal components of these.
  feature_pc_df <- .ComputeTSFeatures(data_full,id_var = "entry",
                                      outcome_var = "target", ...) %>%
    stats::prcomp(scale = TRUE)
  # Combine PC information with treatment status of the TS.
  pc_entry_df <- feature_pc_df$x %>% tibble::as_tibble() %>%
    dplyr::bind_cols(
      Treated =
        factor(data_full %>%
                 dplyr::distinct(!!as.name(id_var), .keep_all = T) %>%
                 dplyr::pull(Treated))
    )
  # Create a base plot with the PCs by treatment.
  feature_overlap_plot <- pc_entry_df %>%
    ggplot2::ggplot(aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(aes(col = Treated)) +
    ggplot2::ggtitle("First 2 PCs of TS Features by Treatment")+ 
    ggplot2::theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"))
  # Add to the plot based on the input type.
  if (type=="rug"){
    feature_overlap_plot=feature_overlap_plot+ 
      ggplot2::geom_rug(aes(color=Treated)) 
  } else if (type=="centroid"){
    feature_overlap_plot <- .CentroidHelper(feature_overlap_plot, pc_entry_df)
  } else if (type=="contour"){
    feature_overlap_plot <- feature_overlap_plot +
      ggplot2::geom_density_2d(aes(color = Treated,alpha=..level..), size=2)
  }
  # Change the color for slightly improved visibility.
  feature_overlap_plot=feature_overlap_plot +
    ggplot2::scale_color_manual(values=c("#E69F00", "#56B4E9"))
  return(feature_overlap_plot)
}

#' Formats and estimates the TS features for each series.
#'
#' @inheritParams TSFeaturesPCScatterplot
#'
#' @return Data frame of the TS features per series, without constant columns.
.ComputeTSFeatures <- function(data_full, id_var, outcome_var,...){
  
  # Split the data into a list of individual time series.
  list_of_ts <- data_full %>%
    dplyr::select(!!as.name(id_var), !!as.name(outcome_var)) %>%
    split(.[[id_var]]) %>%
    furrr::future_map(~ .[[outcome_var]]) %>%
    furrr::future_map(~ ts(.))
  # Compute the TS features for each series, remove constant columns.
  df1_feat <- tsfeatures::tsfeatures(list_of_ts) %>%
    janitor::remove_constant()
  
  return(df1_feat)
}

#' Add PC centroids to plot.
#'
#' @param base_plot Existing scatter plot from TSFeaturesPCScatterplot.
#' @param pc_df Dataframe of PC of TS features by unit, with treatment status.
#'
#' @return Updated plot with included centroids of the PCs by treatment and 
#'    customized title for the distance between the two.
.CentroidHelper <- function(base_plot, pc_df) {
  # Compute the mean of the centroid of PC1 and PC2 for each treatment group.
  pc_centroids <- pc_df %>%
    dplyr::group_by(Treated) %>%
    dplyr::summarise(
      mean_pc1 = mean(PC1),
      mean_pc2 = mean(PC2)
    )
  # Compute the difference between the two centroids.
  centroid_dist <- pc_centroids %>%
    tidyr::pivot_wider(names_from = Treated, 
                       values_from = c(mean_pc1, mean_pc2)) %>%
    dplyr::mutate(
      dist =
        (mean_pc1_0 - mean_pc1_1)^2 + (mean_pc2_0 - mean_pc2_1)^2
    ) %>%
    dplyr::pull(dist)
  # Add centroids to the plot and add relevant subtitle.
  base_plot <- base_plot  +
    ggplot2::geom_point(aes(x = mean_pc1, y = mean_pc2, col = Treated),
                        data = pc_centroids, size = 12, shape = 18 ) +
    ggplot2::labs(subtitle = paste(
      "Centroids have L2 dist:",
      round(centroid_dist, 4)))
  return(base_plot)
}

#' Test whether TS features are significantly different by treatment status.
#'
#' @inheritParams TSFeaturesPCScatterplot
#'
#' @return A tibble displaying the results of hypothesis tests for differences
#'    in ts features by treatment status.
TSFeatureTest <- function(data_full,
                          id_var = "entry",
                          outcome_var = "target",
                          treat_indicator = "treatperiod_0",
                          ...){
  # Identify the treated units.
  data_full <- .IdentifyTreated(data_full=data_full, id_var=id_var,
                                treat_indicator=treat_indicator)
  # Add treatment status to the tibble of TS features per entry.
  feature_treat_df <- .ComputeTSFeatures(data_full,id_var = "entry",
                                         outcome_var = "target", ...) %>%
    tibble::as_tibble() %>%
    dplyr::bind_cols(
      Treated =
        factor(data_full %>%
                 dplyr::distinct(!!as.name(id_var), .keep_all = T) %>%
                 dplyr::pull(Treated))
    ) 
  # Pivot tibble into long form to compute t-tests.
  test_tib=feature_treat_df%>%
    tidyr::pivot_longer(-Treated,
                        names_to = "vars",
                        values_to = "val"
    ) %>%
    dplyr::group_by(vars) %>%
    rstatix::t_test(val ~ Treated) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance()
  # Remove extraneous column.
  test_tib= test_tib %>%
    dplyr::select(-".y.")
  return(test_tib)
}

#' Plot the density of a particular feature by treatment.
#'
#' @inheritParams TSFeaturesPCScatterplot
#' @param feature String name of the TS feature to plot. This should be the
#'    specific name of the feature outputs from the tsfeature package (e.g. 
#'    call "x_acf1" rather than "acf_features").
#'
#' @return 
FeatureDensity <- function(data_full, id_var = "entry",
                           outcome_var = "target",
                           treat_indicator = "treatperiod_0",
                           feature = "x_acf1",
                           ...) {
  # Identify treated entries.
  data_full <- .IdentifyTreated(data_full=data_full, id_var=id_var,
                                treat_indicator=treat_indicator)
  # Obtain TS features for each series, and identify treatment status.
  feature_treat_df <- .ComputeTSFeatures(data_full,id_var = "entry",
                                         outcome_var = "target", ...) %>%
    tibble::as_tibble() %>%
    dplyr::bind_cols(
      Treated = factor(data_full %>%
                         dplyr::distinct(!!as.name(id_var), .keep_all = T) %>%
                         dplyr::pull(Treated)))
  
  feature_density <- ggplot2::ggplot(feature_treat_df) +
    ggplot2::geom_density(aes(x=!!as.name(feature), color=Treated),alpha=0.4) +
    ggplot2::theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black")) +
    ggplot2::scale_color_manual(values=c("#E69F00", "#56B4E9")) +
    ggplot2::labs(title = paste(feature, "Density by Treatment"))
  
  return(feature_density)
  
}

#' Identify and tag the treated units by ID.
#'
#' @inheritParams TSFeaturesPCScatterplot
#'
#' @return The original tibble of data, with a new Treated indicator column. 
.IdentifyTreated <- function(data_full, id_var, treat_indicator) {
  # Identify treated entries.
  treated_entries <- data_full %>%
    dplyr::filter(!!as.name(treat_indicator) > 0) %>%
    dplyr::distinct(!!as.name(id_var)) %>% 
    dplyr::pull()
  data_full <- data_full %>%
    dplyr::mutate(Treated = as.numeric(!!as.name(id_var) %in% treated_entries))
  return(data_full)
}