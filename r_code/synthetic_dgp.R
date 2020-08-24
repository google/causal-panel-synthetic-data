pacman::p_load(
  dplyr, furrr, tibble, tidyr, stats, lubridate,
  glue, truncnorm, MASS
)

#' Generates tibble of long form panel data using a factor-augmented AR 1.
#'
#' @param num_entries Number of units to be generated.
#' @param num_periods Number of time periods to be generated.
#' @param first_treat Integer time of the first possible treatment period. Time
#'    will be drawn from a geometric distribution for each unit, starting here.
#' @param freq String indicating time unit, one of "daily", "weekly", "monthly".
#' @param prop_treated Proportion of entries that should receive treatment.
#' @param treat_impact_mean Initial period treatment impact mean, drawn from
#'    truncated normal distribution centered here. Bounds are [0,0.25] by
#'    default, but shift to [a,b] where b=mean+0.25 if the mean is larger than
#'    0.25 (max of 1) or a=mean-0.25 if mean is below 0.
#' @param treat_impact_sd Standard deviation of the truncated normal impact.
#' @param treat_decay_mean Initial period treatment decay mean, drawn from
#'    truncated normal distribution centered here. Bounds are [0,0.9] by
#'    default, but shift to [0,1+eps] if mean>0.9 -- allowing units to have no
#'    decay (value of 1). Propagates as mean**(post_treat_period).
#' @param treat_decay_sd Standard deviation of the truncated normal decay.
#' @param selection String indicating treatment assignment mechanism; one of
#'    "random", "observables", "unobservables".
#' @param rho Mean of truncated normal distribution for the autocorrelation of
#'    outcome process, with bounds [0, 0.995] by default. If input is negative
#'    (must be >-0.25), lower bound shifts to mean-0.25.
#' @param rho_scale Standard deviation of truncated normal for autocorrelation.
#' @param rho_shift Multiplier on the mean rho for control units. Overall
#'    mean must stay below 1 and above -0.25 (-0.25<rho*rho_shift<1).
#' @param cov_overlap_scale Number between (-1,1) shifting distribution of
#'    covariates for a random fraction (prop_treated) of generated x variables.
#'    A value of 1 shifts the distribution for treatment up on all x's,
#'    whereas -1 shifts up the distribution for donors.
#' @param num_factors Integer number of time-varying, unobserved factors to
#'    simulate. For freq="monthly", "weekly", "daily", number of factors must be
#'     at least 3, 4, and 5 respectively.
#' @param loading_scale Number in (-1,1) shifting loading distribution by
#'    treatment. Negative values shift loadings up for control units, positive
#'     values shift loadings distribution up for treated units.
#' @param intercept_scale Number in (-1,1) that shifts the mean of a truncated
#'    normal distribution intercept by treatment assignment. Positive input
#'    shifts treatment to be larger in size and shifts controls to be smaller.
#' @param conditional_impact_het Constant added to the treatment impact. Applies
#'    only to units in top 25% (added) and bottom 25% (subtracted) of y in time
#'    1. Positive input increases impact for large units.
#' @param rescale_y_mean Number representing the target mean of monthly
#'    exp(counterfactual). Output will be shifted towards this.
#' @param seed Starting seed for random generation processes, for replicability.
#'
#' @return Long form tibble with columns for observed, potential treated, and
#'    potential untreated outcomes, treatment time, and the relevant x variables
#'    (loadings, intercept, observables).
# TODO(alexdkellogg): Allow dgp to be covariates only or factors only?
# TODO(alexdkellogg): Extend to multiple treatment and continuous treatment.
SyntheticDGP <- function(num_entries = 500,
                         num_periods = 120,
                         first_treat = 70,
                         freq = c("monthly", "weekly", "daily"),
                         prop_treated = 0.4,
                         treat_impact_mean = 0.1,
                         treat_impact_sd = 0.1,
                         treat_decay_mean = 0.7,
                         treat_decay_sd = 0.1,
                         selection = c(
                           "random", "observables",
                           "unobservables"
                         ),
                         rho = 0.9,
                         rho_scale = 0.2,
                         rho_shift = 1,
                         cov_overlap_scale = 0,
                         num_factors = 3,
                         loading_scale = 0,
                         intercept_scale = 0,
                         conditional_impact_het = -0,
                         rescale_y_mean = 2.5e4,
                         seed = NULL) {

  # Match the arguments to valid entries.
  if (!is.null(seed)) set.seed(seed)
  freq <- match.arg(freq)
  
  # Require 5% min in both treat and control.
  stopifnot(0.05 <= prop_treated & prop_treated <= 0.95)
  # Require cov_overlap_scale to be between -1 (shift treat down) and 1.
  stopifnot(-1 <= cov_overlap_scale & cov_overlap_scale <= 1)
  # Require 3 factors of more for monthly, 4 for weekly, and 5 daily.
  if (freq == "monthly") stopifnot(num_factors >= 3)
  if (freq == "weekly") stopifnot(num_factors >= 4)
  if (freq == "daily") stopifnot(num_factors >= 5)
  # Require less than 100% TE.
  stopifnot(treat_impact_mean < 1)
  # Require ACF to less than one, and not too negative.
  stopifnot(-0.25 < rho & rho < 1)
  stopifnot(-0.25 < rho * rho_shift & rho * rho_shift < 1)
  # Require the intercept scale to be between -1 and 1.
  stopifnot(-1 <= intercept_scale & intercept_scale <= 1)
  # Require the loading scale to be between -1 and 1.
  stopifnot(-1 <= loading_scale & loading_scale <= 1)
  # Stop if there are too few pre treat periods, 20 is a rule of thumb.
  stopifnot(20 <= first_treat | first_treat < 0.8 * num_periods)

  # Assign covariates and treatment given selection and overlap.
  synth_data_unit <- .UnitLevelSimulation(
    n_inp = num_entries,
    type_inp = selection,
    cov_overlap_inp = cov_overlap_scale,
    loading_scale_inp = loading_scale,
    num_factors_inp = num_factors,
    int_scale_inp = intercept_scale,
    rho_inp = rho,
    rho_scale_inp = rho_scale,
    rho_shift_inp = rho_shift,
    prop_treated_inp = prop_treated,
    treat_start = first_treat,
    num_periods_inp = num_periods,
    impact_mean_inp = treat_impact_mean,
    impact_sd_inp = treat_impact_sd,
    decay_mean_inp = treat_decay_mean,
    decay_sd_inp = treat_decay_sd
  )

  # Work on time varying characteristics.
  synth_data_factors <- .GenerateFactors(
    num_factors_inp = num_factors,
    num_periods_inp = num_periods,
    freq_inp = freq
  )
  # Combine time and unit characteristics into a grid.
  unit_time_grid <- tidyr::expand_grid(
    entry = seq_len(num_entries),
    time = seq_len(num_periods)
  )

  synth_data_full <- unit_time_grid %>%
    dplyr::left_join(synth_data_unit, by = c("entry")) %>%
    dplyr::left_join(synth_data_factors, by = c("time"))

  rescale_y <- dplyr::case_when(
    freq == "monthly" ~ log(rescale_y_mean),
    freq == "weekly" ~ log(rescale_y_mean / 4),
    freq == "daily" ~ log(rescale_y_mean / 30)
  )
  # Generate the counterfactual outcomes for each unit time combo.
  synth_data_full <- .GenerateCounterfactual(synth_data_full,
    num_periods_inp = num_periods,
    rescale_y = rescale_y,
    freq_inp = freq
  )
  # Generate the per period impact (taking account of decay).
  synth_data_full <- .GenerateTreatImpact(
    data_inp = synth_data_full,
    cond_impact_inp =
      conditional_impact_het
  )
  synth_data_full <- synth_data_full %>%
    dplyr::arrange(time, entry)
  return(synth_data_full)
}

#' Draw IID noise and add to the outcome to finalize a draw from DGP.
#'
#' @param data_inp Long form Tibble, output of SyntheticDGP.
#' @param seed Seed from which we draw noise, for reproducibility.
#' @param log_output Boolean indicator for whether output should be log(y).
#' @param sig_y Standard deviation of the iid normal (mean 0) noise.
#'
#' @return Noised long form tibble.
# TODO(alexkellogg): Add non-iid noise process functionality.
NoisifyDraw <- function(data_inp, seed=NULL, log_output = T, sig_y = 0.2) {
  # Draws a vector or random normal (mean zero) epsilon values.
  # This noise is added to the outcome, representing a new iid draw.
  if (!is.null(seed)) set.seed(seed)
  eps <- stats::rnorm(nrow(data_inp), sd = sig_y)
  data_inp <- data_inp %>%
    dplyr::mutate(
      y = exp(y + eps) * (1 - log_output) + log_output * (y + eps),
      y0 = exp(y0 + eps) * (1 - log_output) + log_output * (y0 + eps),
      y1 = exp(y1 + eps) * (1 - log_output) + log_output * (y1 + eps)
    )
  return(data_inp)
}

#' Compute individual-time specific factor loading product.
#'
#' @param data_inp Long form tibble created by private functions containing
#' individual loadings parameters as well as time specific factor.
#'
#' @return Tibble of computed factor-loadings by unit and time.
.ComputeFactorLoadings <- function(data_inp) {
  # Create a list of factor loadings split by individual.
  loadings_vec_list <- data_inp %>%
    dplyr::select(entry, tidyselect::contains("loading")) %>%
    dplyr::distinct(entry, .keep_all = T) %>%
    dplyr::group_split(entry, .keep = F) %>%
    lapply(as.matrix)
  # Create a list of factor matrices split by individual.
  factor_mat_list <- data_inp %>%
    dplyr::select(entry, tidyselect::contains("factor")) %>%
    dplyr::group_split(entry, .keep = F) %>%
    lapply(as.matrix)
  # Compute the matrix product of loadings and factors for each individual.
  factor_loading <- furrr::future_map2(
    .x = loadings_vec_list, .y = factor_mat_list,
    .f = ~ (.y) %*% t(.x)
  ) %>%
    unlist()
  return(factor_loading)
}

#' Generate the counterfactual outcome process given innovations and AR model.
#'
#' @param data_inp Long form tibble containing individual autocorrelation
#'    parameters as well as factor*loadings by unit and time.
#' @param num_periods_inp Number of time periods in the tibble.
#'
#' @return Individual specific counterfactual outcome series.
.ComputeOutcomeProcess <- function(data_inp, num_periods_inp) {
  # Create a list of AR models, with unit specific noise and autocorrelation.
  ar_param_inp <- data_inp %>%
    dplyr::select(entry, autocorr) %>%
    dplyr::distinct(entry, .keep_all = T) %>%
    dplyr::group_split(entry, .keep = F) %>%
    lapply(function(x) {
      list(order = c(1, 0, 0), ar = x[[1]])
    })
  # Create a list of the unit/time specific innovations.
  innov_list <- data_inp %>%
    dplyr::select(entry, tidyselect::contains("loading")) %>%
    dplyr::group_split(entry, .keep = F) %>%
    lapply(as.matrix)
  # Compute the unit series in parallel.
  outcome_series <- furrr::future_map2(
    .x = ar_param_inp, .y = innov_list,
    .f = ~ stats::arima.sim(
      model = .x,
      n = num_periods_inp,
      n.start = 500,
      innov = .y
    )
  ) %>%
    unlist()
  return(outcome_series)
}

#' Combine variables from the model to form counterfactual outcomes.
#'
#' @param data_inp Long form tibble containing individual autocorrelation
#'    parameters as well as factors, loadings, and Xs by unit and time.
#' @param num_periods_inp Number of time periods in the tibble.
#' @param rescale_y Number representing the adjusted (by freq) target mean of
#'    exp(counterfactual). Output will be shifted towards this.
#' @param freq_inp String indicating time unit, one of "daily", "weekly",
#'    "monthly".
#'
#' @return Long form tibble with counterfactual (rescaled) outcome appended.
# TODO(alexdkellogg): Scaling should probably become a function of rho.
.GenerateCounterfactual <- function(data_inp, num_periods_inp, rescale_y,
                                    freq_inp = c("monthly", "weekly", "daily")
                                    ){
  freq_inp <- match.arg(freq_inp)
  # Compute each of the component pieces for generating counterfactual y.
  outcome_series <- data_inp %>%
    dplyr::select(
      time, entry, intercept, autocorr,
      tidyselect::matches("loading|factor")
    )
  # For each unit and time, compute the factor-loading.
  computed_factor_vec <- .ComputeFactorLoadings(outcome_series)
  outcome_series <- outcome_series %>%
    dplyr::select(-tidyselect::matches("loading|factor")) %>%
    dplyr::mutate(factor_loadings = computed_factor_vec)
  # Compute the AR process for y0 given individual rho, factor-loadings.
  outcome_ar <- .ComputeOutcomeProcess(outcome_series, num_periods_inp)

  # Compute the covariate contribution to outcome.
  xvars <- data_inp %>%
    dplyr::distinct(entry, .keep_all = T) %>%
    dplyr::select(tidyselect::matches("obs")) %>%
    as.matrix()
  # Multiply the covariates by a random coefficient to aggregate.
  x_beta <- xvars %*% stats::rnorm(n = ncol(xvars), sd = 0.5) %>%
    as.vector()
  # Map the aggregated covariates to the units.
  x_beta_x_walk <- tibble::tibble(
    entry = seq_len(length(x_beta)),
    xbeta = x_beta
  )
  # Append these computed components to our tibble.
  outcome_series <- outcome_series %>%
    dplyr::left_join(x_beta_x_walk, by = "entry") %>%
    dplyr::mutate(
      factor_loadings = computed_factor_vec,
      indiv_component = intercept + xbeta,
      y0 = outcome_ar + indiv_component
    )
  # Adjust the tightness of the rescaled distribution by frequency.
  freq_scaler <- dplyr::case_when(
    freq_inp == "monthly" ~ 0.7,
    freq_inp == "weekly" ~ 0.9,
    freq_inp == "daily" ~ 1.1
  )

  # Rescale the outcome series to avoid massive outliers.
  outcome_series_rescaled <- outcome_series %>%
    dplyr::mutate(
      de_sd_y0 = y0 / (freq_scaler * sd(y0)),
      rescaled_y = (de_sd_y0 - mean(de_sd_y0)) + rescale_y
    )
  # Keep outcomes positive for easy aggregation, if desired.
  min_rescaled <- min(outcome_series_rescaled[["rescaled_y"]])
  if (min_rescaled < 0) {
    outcome_series_rescaled <- outcome_series_rescaled %>%
      dplyr::mutate(rescaled_y = rescaled_y - min_rescaled)
  }
  outcome_series_rescaled <- outcome_series_rescaled %>%
    dplyr::mutate(y0 = rescaled_y) %>%
    dplyr::select(time, entry, factor_loadings, y0)
  data_inp <- data_inp %>%
    dplyr::inner_join(outcome_series_rescaled, by = c("time", "entry")) %>%
    dplyr::select(
      time, entry, treated, treatment_period, y0,
      factor_loadings, dplyr::everything()
    )
  return(data_inp)
}

#' Compute and append treatment impact by unit and time.
#'
#' @param data_inp Long form tibble already countaining counterfactual y.
#' @param cond_impact_inp Constant added to the treatment impact. Applies
#'    only to units in top 25% (added) and bottom 25% (subtracted) of y in time
#'    1. Positive input increases impact for large units.
#'
#' @return Long form tibble with treatment impact, potential treated outcome,
#'    and observed outcome appended.
.GenerateTreatImpact <- function(data_inp,
                                 cond_impact_inp) {
  # Determine which observations get a conditional treatment boost.
  cond_x_walk <- data_inp %>%
    dplyr::filter(time == 1) %>%
    dplyr::select(entry, y0) %>%
    dplyr::mutate(
      cond_treat_impact = dplyr::case_when(
        y0 > quantile(y0, 0.75) ~ cond_impact_inp,
        y0 < quantile(y0, 0.25) ~ -cond_impact_inp,
        TRUE ~ 0
      )
    ) %>%
    dplyr::select(-y0)
  
  # Define the treatment propagation for each unit and time combo.
  data_inp <- data_inp %>%
    dplyr::left_join(cond_x_walk, by = "entry") %>%
    dplyr::group_by(entry) %>%
    dplyr::mutate(
      post_treat_t = time - treatment_period,
      decay_t = dplyr::case_when(
        post_treat_t < 0 ~ 0,
        post_treat_t >= 0 ~ treat_decay**post_treat_t
      ),
      impact_t = decay_t * (treat_impact + cond_treat_impact),
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(treat_decay, treat_decay, decay_t, cond_treat_impact))
  
  # Add the treatment impact to create potential treated outcome (y1).
  # Also note the empirically observed y.
  data_inp <- data_inp %>%
    dplyr::mutate(
      y1 = impact_t + y0,
      y = treated * y1 + (1 - treated) * y0
    ) %>%
    dplyr::select(
      time, entry, treated,
      post_treat_t, treatment_period, impact_t,
      y, y0, y1,
      dplyr::everything()
    )
  return(data_inp)
}

#' Create a grid of date times.
#'
#' @param num_periods_inp Number of time periods to be generated.
#' @param freq_inp String indicating time unit; "daily", "weekly", "monthly".
#'
#' @return Grid in freq_inp units spanning num_period_inp periods, with a
#'    default start date of 2002-01-01.
.GenerateTimeGrid <- function(num_periods_inp,
                              freq_inp = c("monthly", "weekly", "daily")) {
  freq_inp <- match.arg(freq_inp)
  # Creates a tibble that spans the input dates by the input frequency.
  # Stores the relevant info about that date.
  date_start_temp <- "2002-01-01"
  period_dates <- switch(freq_inp,
    "daily" = format(lubridate::ymd(date_start_temp) +
      lubridate::days(0:(num_periods_inp - 1))),
    "weekly" = format(lubridate::ymd(date_start_temp) +
      lubridate::weeks(0:(num_periods_inp - 1))),
    "monthly" = format(lubridate::ymd(date_start_temp) +
      months(0:(num_periods_inp - 1)))
  )

  # Identify the relevant components of the date (day/week/etc).
  date_info_tib <- tibble::tibble(
    time = seq_len(num_periods_inp),
    date_t = period_dates,
    day_num = lubridate::day(period_dates),
    yday_num = lubridate::yday(period_dates),
    week_num = lubridate::week(period_dates),
    month_num = lubridate::month(period_dates),
    quarter_num = lubridate::quarter(period_dates),
    year_num = lubridate::year(period_dates)
  )

  return(date_info_tib)
}

#' Creates the requisite number of factors beyond trend, given frequency.
#'
#' @param date_tib Tibble with columns for the date time, as well as for each
#'    component of the date (e.g. day) spanning the desired number of periods.
#' @param shock_name String indicating the type of shock, one of "quarter_num",
#'    "month_num", "week_num", "day_num".
#' @param factor_name String indicating the desired factor name.
#' @param ar_model AR model list object for AR(1) process with rho=0.2.
#' @param ar_sd AR(1) parameter noise term for the factors.
#'
#' @return Tibble with the shocks by shock_name, as well as the factors.
.BaseFactorHelper <- function(date_tib, 
                              shock_name= c("quarter_num", "month_num", 
                                            "week_num", "day_num"), 
                              factor_name,
                              ar_model, ar_sd) {
  shock_name <- match.arg(shock_name)
  # Generate the bounds for the random uniform noise.
  lim <- switch(shock_name,
    "day_num" = 0.1,
    "week_num" = 0.5,
    1
  )
  # Create a map of the distinct shock frequency and the associated shock.
  f_shocks <- paste(shock_name[1], "shock", sep = "_")
  shock_map <- date_tib %>%
    dplyr::select(tidyselect::all_of(shock_name)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(!!as.name(f_shocks) := stats::runif(dplyr::n(), -lim, lim))
  # Map the shocks into the factors following the AR model.
  tib_out <- date_tib %>%
    dplyr::inner_join(shock_map, by = shock_name) %>%
    dplyr::group_by(!!as.name(shock_name)) %>%
    dplyr::mutate(
      !!as.name(factor_name) := stats::arima.sim(
        model = ar_model, n = dplyr::n(),
        innov = !!as.name(f_shocks) + stats::rnorm(dplyr::n(),
          sd = ar_sd
        ),
        n.start = 500
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(tidyselect::all_of(factor_name))

  return(tib_out)
}

#' Creates extra noisy factors, beyond the base requisiste number.
#'
#' @param date_tib Grid of date time, spanning the desired number of periods.
#' @param factor_name String indicating the desired factor name.
#' @param ar_model AR model list object for AR(1) process with rho=0.2.
#' @param ar_sd AR(1) parameter noise term for the factors.
#'
#' @return Tibble with the shocks by shock_name, as well as the factors.
.ExtraFactorHelper <- function(date_tib, factor_name, ar_model, ar_sd) {
  # Generate a set of unique dates which will be mapped to a shock.
  date_tib <- date_tib %>%
    dplyr::group_by(year_num) %>%
    dplyr::mutate(extra_t = dplyr::row_number()) %>%
    dplyr::ungroup()
  unique_dates <- date_tib %>%
    dplyr::select(time, extra_t) %>%
    dplyr::distinct(extra_t)

  effective_t <- max(unique_dates)
  num_shocks <- sample(3:min(25, effective_t), 1) - 1

  shock_map <- unique_dates %>%
    sample_n(num_shocks) %>%
    dplyr::mutate(e_shock = stats::runif(dplyr::n(), -1, 1)) %>%
    dplyr::full_join(unique_dates, by = "extra_t") %>%
    dplyr::mutate(e_shock = tidyr::replace_na(e_shock, 0)) %>%
    dplyr::arrange(extra_t)

  # Map the shocks into factors, and structure the output.
  tib_out <- date_tib %>%
    dplyr::inner_join(shock_map, by = "extra_t") %>%
    dplyr::group_by(e_shock) %>%
    dplyr::mutate(
      !!as.name(factor_name) := stats::arima.sim(
        model = ar_model, n = dplyr::n(),
        innov = e_shock + stats::rnorm(dplyr::n(),
          sd = ar_sd
        ),
        n.start = 500
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(tidyselect::all_of(factor_name))

  return(tib_out)
}

#' Generate time-varying factors.
#'
#' @param num_factors_inp Integer number of time-varying, unobserved factors to
#'    simulate. For freq="monthly", "weekly", "daily", number of factors must be
#'     at least 3, 4, and 5 respectively.
#' @param num_periods_inp Number of time periods to be generated.
#' @param freq_inp String indicating time unit: "daily", "weekly", "monthly".
#' @param ar_rho AR(1) parameter for the factors, fixed to rho=0.2.
#' @param ar_sd AR(1) parameter noise term for the factors, 0.75
#'
#' @return Time only tibble with factors appended.
# TODO(alexdkellogg): for some reason, get very different results when I
#    specify the sd=sqrt(ar_sd) into the ar_model -- should be same though!
.GenerateFactors <- function(num_factors_inp, num_periods_inp, 
                             freq_inp = c("monthly", "weekly", "daily"),
                             ar_rho = 0.2, ar_sd = 0.75) {
  freq_inp <- match.arg(freq_inp)
  # AR 1 with auto correlation and sd inputs.
  ar_model <- list(order = c(1, 0, 0), ar = ar_rho)

  # Create the date indicator grid along with the first factor (trend).
  factor_tib <- .GenerateTimeGrid(
    num_periods_inp = num_periods_inp,
    freq_inp = freq_inp
  ) %>%
    dplyr::mutate(factor1 = time / dplyr::n() +
      stats::rnorm(dplyr::n(), mean = 0, sd = 0.1))

  # Set up the input for the requisite base number of factors.
  base_vec <- switch(freq_inp,
    "daily" = c("quarter_num", "month_num", "week_num", "day_num"),
    "weekly" = c("quarter_num", "month_num", "week_num"),
    "monthly" = c("quarter_num", "month_num")
  )
  # Generate the base factors.
  shocks_tib <- purrr::map2_dfc(
    .x = base_vec,
    .y = c(glue::glue("factor{2:(length(base_vec)+1)}")),
    .f = ~ .BaseFactorHelper(
      date_tib = factor_tib, shock_name = .x,
      factor_name = .y, ar_model = ar_model,
      ar_sd = ar_sd
    )
  )

  factor_tib <- factor_tib %>% dplyr::bind_cols(shocks_tib)
  
  if (num_factors_inp > length(base_vec) + 1) {
    # Add additional factors beyond the number implied by freq_inp.
    extra_factors_names <-
      c(glue::glue("factor{(length(base_vec)+2):num_factors_inp}"))

    extra_factor_tib <-
      purrr::map_dfc(
        .x = extra_factors_names,
        .f = ~ .ExtraFactorHelper(
          date_tib = factor_tib,
          factor_name = .x, ar_model = ar_model,
          ar_sd = ar_sd
        )
      )

    factor_tib <- factor_tib %>%
      dplyr::bind_cols(extra_factor_tib)
  }
  # Keep just the factors and the date, for later merging.
  factor_tib <- factor_tib %>% dplyr::select(-tidyselect::contains("shock"))
  return(factor_tib)
}

#' Create unit specific (time constant) variables.
#'
#' @param n_inp Number of units to be generated.
#' @param type_inp String indicating treatment assignment mechanism; one of
#'    "random", "observables", "unobservables".
#' @param cov_overlap_inp Number between (-1,1) shifting distribution of
#'    covariates for a random fraction (prop_treated) of generated x variables.
#'    A value of 1 shifts the distribution for treatment up on all x's,
#'    whereas -1 shifts up the distribution for donors.
#' @param loading_scale_inp Number in (-1,1) shifting loading distribution by
#'    treatment. Negative values shift loadings up for control units, positive
#'    values shift loadings distribution up for treated units.
#' @param num_factors_inp Integer number of time-varying, unobserved factors to
#'    simulate. For freq="monthly", "weekly", "daily", number of factors must be
#'    at least 3, 4, and 5 respectively.
#' @param int_scale_inp Number in (-1,1) that shifts the mean of a truncated
#'    normal distribution intercept by treatment assignment. Positive input
#'    shifts treatment to be larger in size and shifts controls to be smaller.
#' @param rho_inp Mean of truncated normal distribution for the autocorrelation
#'    of outcome process, with default bounds [0, 0.995]. If input is negative
#'    (must be >-0.25). lower bound shifts to mean-0.25.
#' @param rho_scale_inp Standard deviation of truncated normal for
#'    autocorrelation.
#' @param rho_shift_inp Multiplier on the mean rho for control units. Overall
#'    mean must stay below 1 and above -0.25 (-0.25<rho*rho_shift<1).
#' @param prop_treated_inp Proportion of entries that should receive treatment.
#' @param treat_start Integer time of the first possible treatment period. Time
#'    will be drawn from a geometric distribution for each unit, starting here.
#' @param num_periods_inp Number of time periods to be generated.
#' @param impact_mean_inp Initial period treatment impact mean, drawn from
#'    truncated normal distribution centered here. Bounds are [0,0.25] by
#'    default, but shift to [a,b] where b=mean+0.25 if the mean is larger than
#'    0.25 (max of 1)  or a=mean-0.25 if mean is below 0.
#' @param impact_sd_inp Standard deviation of the truncated normal impact.
#' @param decay_mean_inp Initial period treatment decay mean, drawn from
#'    truncated normal distribution centered here. Bounds are [0,0.9] by
#'    default, but shift to [0,1+eps] if mean>0.9 -- allowing units to have no
#'    decay (value of 1). Propagates as mean**(post_treat_period).
#' @param decay_sd_inp Standard deviation of the truncated normal decay.
#'
#' @return Tibble with all individual level work completed.
.UnitLevelSimulation <- function(n_inp,
                                 type_inp = c(
                                   "random", "observables",
                                   "unobservables"
                                 ),
                                 cov_overlap_inp,
                                 loading_scale_inp,
                                 num_factors_inp,
                                 int_scale_inp,
                                 rho_inp,
                                 rho_scale_inp,
                                 rho_shift_inp,
                                 prop_treated_inp,
                                 treat_start,
                                 num_periods_inp,
                                 impact_mean_inp,
                                 impact_sd_inp,
                                 decay_mean_inp,
                                 decay_sd_inp) {
  type_inp <- match.arg(type_inp)
  # Gather data that is time invariant.
  # First, call to assign treatment, which generates covariates.
  unit_level_tib <- .AssignTreat(
    n_inp = n_inp,
    type_inp = type_inp,
    cov_overlap_inp = cov_overlap_inp,
    loading_scale_inp = loading_scale_inp,
    num_factors_inp = num_factors_inp,
    rho_inp = rho_inp,
    rho_scale_inp = rho_scale_inp,
    rho_shift_inp = rho_shift_inp,
    int_scale_inp = int_scale_inp,
    prop_treated_inp = prop_treated_inp
  )

  # Merge in the treatment period assignment.
  unit_level_tib <- .AssignTreatTime(
    unit_level_tib, treat_start,
    num_periods_inp
  )
  # Allow wiggle room for the TE without getting NA.
  impact_ub <- ifelse(impact_mean_inp < 0.25, 0.25,
    min(1, impact_mean_inp + (0.25))
  )
  impact_lb <- ifelse(impact_mean_inp <= 0, impact_mean_inp - 0.25, 0)
  decay_ub <- ifelse(decay_mean_inp > 0.9, 1.000000001, 0.9)
  # Assign treatment effect impact and decay parameters per unit.
  unit_level_tib <- unit_level_tib %>%
    dplyr::mutate(
      treat_impact = truncnorm::rtruncnorm(
        n = dplyr::n(), a = impact_lb, b = impact_ub,
        mean = impact_mean_inp, sd = impact_sd_inp
      ),
      treat_decay = truncnorm::rtruncnorm(
        n = dplyr::n(), a = 0, b = decay_ub,
        mean = decay_mean_inp, sd = decay_sd_inp
      )
    )
  # Require that no treat decay is larger than 1.
  unit_level_tib <- unit_level_tib %>%
    dplyr::mutate(treat_decay = dplyr::case_when(
      treat_decay > 1 ~ 1,
      treat_decay <= 1 ~ treat_decay
    ))

  return(unit_level_tib)
}

#' Treatment assignment mechanism, with shifted covariates appended.
#'
#' @inheritParams .UnitLevelSimulation
#'
#' @return Tibble with treatment assignment and individual specific covariates
#'    shifted conditional on treatment assignment and inputs.
.AssignTreat <- function(n_inp, 
                         type_inp = c( "random",
                                       "observables", 
                                       "unobservables"),
                         cov_overlap_inp,
                         loading_scale_inp, rho_inp,
                         rho_scale_inp, rho_shift_inp, num_factors_inp,
                         int_scale_inp, prop_treated_inp) {
  type_inp <- match.arg(type_inp)
  # Generate covariates with varying levels of overlap for shifted subset.
  covariate_tib <- .GenCovariates(
    n_inp = n_inp, cov_overlap_inp = cov_overlap_inp,
    frac_shifted = prop_treated_inp
  )

  if (type_inp == "random") {
    # Randomly assign treatment to prop_treated_inp.
    treat_covar_tib <- tibble::tibble(
      entry = seq_len(n_inp),
      treated = as.numeric(dplyr::percent_rank(
        runif(n = n_inp, min = 0, max = 1)
      ) > 1 - prop_treated_inp)
    ) %>%
      inner_join(covariate_tib, by = "entry")
  }
  else {
    # Assign treatment based on either observable or unobservable Xs.
    relevant_xs <- ifelse(type_inp == "observables", "obs", "unobs")
    # Select the subset of relevant covariates given selection type.
    z <- covariate_tib %>%
      dplyr::select(tidyr::starts_with(relevant_xs)) %>%
      as.matrix()

    # Combine the variables for each observation to form a predicted score.
    # Rescale with logistic function to map into the probability space (0,1).
    prob_treat <- tibble::tibble(
      score =
        z %*% stats::rnorm(n = ncol(z), mean = 0, sd = 1)
    ) %>%
      dplyr::mutate(prop_score = 1 / (1 + exp(-score)))

    # Add random noise to the score and take the top fraction as treated.
    treat_covar_tib <- prob_treat %>%
      dplyr::mutate(
        u = runif(n = n_inp, min = 0, max = 1),
        p_new = prop_score + u,
        treated = as.numeric(dplyr::percent_rank(p_new) >
          1 - prop_treated_inp),
        entry = seq_len(n_inp)
      ) %>%
      dplyr::select(entry, treated) %>%
      inner_join(covariate_tib, by = "entry")
  }
  # Generate loadings and intercept given treatment assignment.
  # Allows these to be correlated with treatment, for overlap assumptions.
  unit_tib <- .GenerateLoadings(
    treat_tib_inp = treat_covar_tib,
    loading_scale_inp = loading_scale_inp,
    num_factors_inp = num_factors_inp,
    int_scale_inp = int_scale_inp,
    rho_inp = rho_inp,
    rho_scale_inp = rho_scale_inp,
    rho_shift_inp = rho_shift_inp
  )
  return(unit_tib)
}

#' Generates unit specific loadings, intercepts, and autocorrelation.
#'
#' @param treat_tib_inp Unit level tibble with treatment assignment column.
#' @inheritParams .UnitLevelSimulation
#'
#' @return Tibble with loadings, intercepts, and autocorrelation assigned based
#'    on selection mechanism.
.GenerateLoadings <- function(treat_tib_inp, loading_scale_inp, num_factors_inp,
                              int_scale_inp, rho_inp, rho_scale_inp,
                              rho_shift_inp) {
  # Generate loadings for each unit and each factor.
  loadings_mat <- matrix(0, nrow = nrow(treat_tib_inp), ncol = num_factors_inp)
  colnames(loadings_mat) <- glue::glue("loading{1:num_factors_inp}")
  # In case rho mean is negative, assign a flexible lower bound to autocorr.
  rho_lb <- ifelse(rho_inp <= 0, rho_inp - 0.25, 0)

  # Assign potentially shifted rho and intercept for each individual.
  tib_pre_loadings <- treat_tib_inp %>%
    dplyr::group_by(treated) %>%
    dplyr::mutate(
      intercept = truncnorm::rtruncnorm(
        n = dplyr::n(), a = 4, b = 11,
        mean = 7.5 - 2 * (1 - treated) * int_scale_inp +
          2 * treated * (int_scale_inp),
        sd = 0.7
      ),
      autocorr = truncnorm::rtruncnorm(
        n = dplyr::n(), a = rho_lb, b = 0.995,
        mean = ifelse(treated, rho_inp,
          rho_inp * rho_shift_inp
        ),
        sd = rho_scale_inp
      )
    ) %>%
    dplyr::bind_cols(tibble::as_tibble(loadings_mat)) %>%
    dplyr::ungroup()
  # Assign potentially shifted loadings by unit and treatment status.
  tib_pre_loadings <- tib_pre_loadings %>%
    dplyr::group_by(treated) %>%
    dplyr::mutate(
      dplyr::across(tidyselect::starts_with("loading"),
        .fns = list(
          load = ~ stats::rbeta(
            dplyr::n(),
            2 - (1 - treated) * (loading_scale_inp),
            2 - (treated) * (loading_scale_inp)
          )
        ),
        .names = "{col}"
      )
    ) %>%
    dplyr::ungroup()

  return(tib_pre_loadings)
}

#' Generate observable and unobservable Xs for use in selection.
#'
#' @param n_inp Number of units.
#' @param cov_overlap_inp Number between (-1,1) shifting distribution of
#'    covariates for a random fraction (prop_treated) of generated x variables.
#'    A value of 1 shifts the distribution for treatment up on all x's,
#'    whereas -1 shifts up the distribution for donors.
#' @param frac_shifted Fraction of units (random) to receive covariate shift.
#'
#' @return A tibble with observable and unobservable x variables.
# TODO(alexkellogg): Generalize to allow for user to input number of xs,
#    also, to break IFE models, allow for time-varying x's uncorrelated to
#    factors.
.GenCovariates <- function(n_inp, cov_overlap_inp, frac_shifted) {
  # Create a set of correlated variables, mean zero, var/covar sigma.
  sigma <- matrix(c(
    16, 4, -4.8,
    4, 25, 9,
    -4.8, 9, 9
  ), 3, 3)
  corr_xs <- MASS::mvrnorm(n = n_inp, rep(0, 3), sigma)
  colnames(corr_xs) <- c("obs_mvnorm1", "unobs_mvnorm", "obs_mvnorm2")
  # Generate covariates, one observed and unobserved from beta, binom, and norm.
  x_tib <- tibble::tibble(
    entry = seq_len(n_inp),
    to_shift = as.numeric(dplyr::percent_rank(
      runif(n = n_inp, min = 0, max = 1)
    ) > 1 - frac_shifted)
  ) %>%
    cbind(corr_xs) %>%
    dplyr::group_by(to_shift) %>%
    dplyr::mutate(
      obs_beta = stats::rbeta(
        n = dplyr::n(),
        shape1 = 6,
        shape2 = 5 + 3 * (to_shift * -(cov_overlap_inp))
      ),
      obs_binom = stats::rbinom(
        n = dplyr::n(),
        size = 15,
        prob = 0.5 + 0.25 * (to_shift * (cov_overlap_inp))
      ),
      obs_norm = stats::rnorm(
        n = dplyr::n(),
        mean = 0 + 5 * (to_shift * (cov_overlap_inp)),
        sd = 10
      ),
      unobs_beta = stats::rbeta(
        n = dplyr::n(),
        shape1 = 6,
        shape2 = 5 + 3 * (to_shift * -(cov_overlap_inp))
      ),
      ubobs_binom = stats::rbinom(
        n = dplyr::n(),
        size = 15,
        prob = 0.5 + 0.25 * (to_shift * (cov_overlap_inp))
      ),
      ubobs_norm = stats::rnorm(
        n = dplyr::n(),
        mean = 0 + 5 * (to_shift * (cov_overlap_inp)),
        sd = 10
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-to_shift)

  return(x_tib)
}

#' Assign treatment time to each unit.
#'
#' @param treat_tib_inp Unit level tibble with treatment assignment column.
#' @param treat_start Integer time of the first possible treatment period. Time
#'    will be drawn from a geometric distribution for each unit, starting here.
#' @param num_periods_inp Number of time periods total.
#'
#' @return Tibble with treatment time randomly assigned from geometric
#'    distribution, with probability shifting to ensure adequate post periods.
.AssignTreatTime <- function(treat_tib_inp, treat_start, num_periods_inp) {
  # Geometric distribution parameter selected based on time.
  # Goal is to have most treatment assigned before end date.
  geom_prob <- dplyr::case_when(
    num_periods_inp - treat_start <= 20 ~ 0.25,
    num_periods_inp - treat_start <= 30 ~ 0.2,
    num_periods_inp - treat_start <= 40 ~ 0.15,
    TRUE ~ 0.1
  )
  # Draw treatment period as first date plus geometric random variable.
  # If the assignment date is past the final period, change it to final period.
  treat_tib_inp <- treat_tib_inp %>%
    dplyr::mutate(
      treatment_period = treat_start +
        stats::rgeom(dplyr::n(), geom_prob)
    ) %>%
    dplyr::group_by(entry) %>%
    dplyr::mutate(
      treatment_period =
        min(treatment_period, num_periods_inp)
    ) %>%
    dplyr::ungroup()
  return(treat_tib_inp)
}

#' Reformats output of SyntheticDGP for easy estimation using our functions.
#'
#' @param data_output Long form tibble, output of SyntheticDGP after noised.
#'
#' @return A reformatted tibble with default column names.
FormatForEst <- function(data_output) {
  # Reformats the synthetic data for default estimation by my estimators.
  # For estimation, we want a dummy for the treated units in treated times.
  reformatted_data <- data_output %>%
    dplyr::mutate(treatperiod_0 = (treated == 1) * (post_treat_t >= 0)) %>%
    dplyr::rename(
      potential_post_t = post_treat_t,
      potential_treatment_period = treatment_period,
      potential_impact_t = impact_t,
      target = y,
      counter_factual = y0,
      period = time
    ) %>%
    dplyr::select(
      period, entry, treated, treatperiod_0, target,
      counter_factual
    )

  return(reformatted_data)
}

#' Transforms SyntheticDGP output into an encouragement design.
#'
#' @param data_output Long form tibble, output of SyntheticDGP after noised.
#' @param complier_frac Fraction of units assigned to treatment that actually
#'    take-up the treatment. Also, 1-complier_frac of those not assigned to
#'    treatment take-up treatment regardless.
#' @param selection_vars Variables on which take-up is based, as a vector of
#'    column names. 
#' @param seed Random number seed for replicability of who takes-up treatment.
#'
#' @return A tibble that takes the treatment assignment from data_output and
#'    perturbs it into an encouragement design, where entries decide themselves
#'    whether to accept their treatment assignment. New columns are 
#'    encouragement, the assigned treatment group, and treated, which is the 
#'    actual take-up for that unit -- which depends on the fraction of compliers
#'    and the selection mechanism into take-up. The observed outcome, y, is also
#'    altered to reflect the take-up of treatment rather than the assignment.
EncouragementDesignDGP <- function(data_output, complier_frac=0.7,
                                   selection_vars=NULL, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  # Grab the columns relevant to the encouragement design definitions.
  itt_tib <- data_output %>%
    dplyr::distinct(entry, .keep_all=T) %>%
    dplyr::select(tidyselect::all_of(c("entry", "treated", "y0", "y1",
                                       selection_vars))) %>%
    dplyr::mutate(sel_score = treated) %>% 
    dplyr::rename(encouraged="treated")
  # If selection into take up is desired, model it.
  if (!is.null(selection_vars)) {
    selection_mat <- itt_tib %>% 
      dplyr::select(tidyselect::all_of(selection_vars)) %>%
      as.matrix()
    random_proj <- rnorm(n=ncol(selection_mat), 0, 1)
    proj_vec <- as.vector(selection_mat %*% random_proj)
    # Create an N x 1 vector of selection scores to rank each unit.
    itt_tib <- itt_tib %>% dplyr::mutate(sel_score = sel_score + proj_vec)
  }
  # Add additional random noise to the selection score, assign take up.
  itt_tib <- itt_tib %>%
    dplyr::group_by(encouraged) %>%
    dplyr::mutate(
      p_score = 1 / (1 + exp(-sel_score)),
      u = runif(n = dplyr::n(), min = 0, max = 1),
      p_score = p_score + u,
      treated = as.numeric(dplyr::percent_rank(p_score) > 
                             encouraged * (1-complier_frac) + 
                             (1-encouraged) * complier_frac) 
      ) %>%
    dplyr::ungroup() %>%
    dplyr::select(entry, encouraged, treated)
  # Combine all the data together for output.
  data_out <- data_output %>% 
    dplyr::select(-treated) %>%
    dplyr::left_join(itt_tib, by = "entry") %>%
    dplyr::mutate(y = treated*y1 + (1-treated)*y0) %>%
    dplyr::select(time, entry, encouraged, treated, dplyr::everything())
 return(data_out)
}

#' Aggregate per week, month, quarter data to a coarser frequency.
#'
#' @param disagg_data Data tibble output of SyntheticDGP.
#' @param logged_input Boolean indicating whether the outcomes in the input
#'    tibble are logged or not - which is relevant for aggregating over time.
#' @param to String for desired output frequency: "weekly", "monthly", "yearly".
#' @param from String for current frequency: "daily", weekly", "monthly".
#' @param id_var String name of unit identifier -- default "entry".
#' @param time_var String name of time identifier -- default "time".
#'
#' @return Tibble with covariates and output aggregated from frequency "from"
#'    to frequency "to". Output for the outcome variables is return in logs.
AggregateDataByFreq <- function(disagg_data, logged_input = T, to = "monthly",
                                from = "daily", id_var = "entry",
                                time_var = "time") {
  date_levels <- factor(c("daily", "weekly", "monthly", "yearly"),
    ordered = T,
    levels = c("daily", "weekly", "monthly", "yearly")
  )
  # If no aggregation is happening, break.
  stopifnot(match.arg(to, date_levels) > match.arg(from, date_levels))

  agg_ids <- c("day_num", "week_num", "month_num", "year_num")

  # Identify the larger dates.
  to_ind <- match(to, date_levels)
  to_vars <- agg_ids[seq(to_ind, length(agg_ids))]

  disagg_data <- disagg_data %>% dplyr::mutate(period = !!as.name(time_var))

  # Create the aggregated tibble of outcomes.
  new_date_name <- paste0(to, "_date_end", sep = "")
  agg_data <- disagg_data %>%
    dplyr::mutate(treatperiod_0 = (treated == 1) * (post_treat_t >= 0)) %>%
    dplyr::group_by_at(c(id_var, to_vars)) %>%
    summarise(
      y0 = log(logged_input * sum(exp(y0)) + (1 - logged_input) * sum(y0)),
      y1 = log(logged_input * sum(exp(y1)) + (1 - logged_input) * sum(y1)),
      y = log(logged_input * sum(exp(y)) + (1 - logged_input) * sum(y)),
      !!as.name(new_date_name) := max(date_t),
      treatperiod_0 = max(treatperiod_0)
    ) %>%
    dplyr::ungroup()

  agg_data <- agg_data %>%
    dplyr::mutate(period = agg_data %>%
      dplyr::group_by_at(to_vars) %>%
      dplyr::select(tidyselect::all_of(to_vars)) %>%
      dplyr::group_indices())
  # Create the aggregated tibble with relevant covariates.
  agg_data_covariates <- agg_data %>%
    dplyr::inner_join(
      disagg_data %>% dplyr::select(-c(
        period, date_t, treatment_period,
        post_treat_t, time,
        setdiff(agg_ids, to_vars),
        y0, y1, y
      )) %>%
        distinct_at(c("entry", to_vars), .keep_all = T),
      by = c(id_var, to_vars)
    ) %>%
    dplyr::group_by(!!as.name(id_var)) %>%
    dplyr::mutate(
      Treatment_Period = case_when(
        treated == 0 ~ NA_real_,
        treated == 1 ~ (length(treatperiod_0) - sum(treatperiod_0) + 1)
      ),
      post_treat_t = period - Treatment_Period
    ) %>%
    dplyr::ungroup()

  return(agg_data_covariates)
}
