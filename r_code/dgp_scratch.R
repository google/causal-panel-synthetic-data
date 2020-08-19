pacman::p_load(dplyr, furrr, tidyr, stats,
               caret, glmnet, bigstatsr, forecast)


library(dbarts)

#Goal: rewrite factor creations
#if you want weekly, need min 4 factors; daily minimum 5 factos
#factor 4 becomes a weekly effect, factor 5 becomes a daily effect
#create a grid that is daily by default. append to factor_tib as normal
#Use the structure in factor tib to build a helper?
#For factor 4, regardless of how the data is frequent, grab the distinct week_nums
#and assign 13 shocks (the rest are 0).
#for factor 5, grab the distinct day nums and assign 20 shocks at random locations
#note that if the data is monthly, the week num and day nums will have the same
#number of distinct values -- that's not a problem.
#For all factors>5, grab the distinct day num (most granular) and assign 22
# more shocks shocks

.BaseFactorHelper <- function(date_tib, freq_inp, shock_name, factor_name,
                          ar_model){
  lim=switch (shock_name,
    "day_num" = 0.1,
    "week_num" = 0.2,
    1
  )
  f_shocks=paste(shock_name[1], "shock", sep="_")
  shock_map= date_tib %>% 
    dplyr::select(tidyselect::all_of(shock_name)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(!!as.name(f_shocks):=stats::runif(dplyr::n(), -lim, lim ))
  
  tib_out= date_tib %>%
    dplyr::inner_join(shock_map, by=shock_name) %>%
    dplyr::group_by(!!as.name(shock_name)) %>% 
    dplyr::mutate(
      !!as.name(factor_name) := stats::arima.sim(
        model = ar_model, n = dplyr::n(),
        innov = !!as.name(f_shocks) + stats::rnorm(dplyr::n(),
                                       sd = 0.1
        ),
        n.start = 500 ) ) %>%  dplyr::ungroup() %>%
    dplyr::select(tidyselect::all_of(c(factor_name, f_shocks)))
  
  return(tib_out)
}


.ExtraFactorHelper <- function(date_tib,factor_name, ar_model){
 
  date_tib=date_tib %>%
    dplyr::group_by(year_num) %>%
    dplyr::mutate(extra_t=dplyr::row_number()) %>%
    dplyr::ungroup() 
  unique_dates=date_tib %>%
    dplyr::select(time, extra_t) %>% 
    dplyr::distinct(extra_t) 
  
  effective_t=max(unique_dates)
  num_shocks <- sample(1:min(25,effective_t), 1)
  
  #TODO(alexdkellogg): which of these approaches is right??
  #one: only have a shock on the particular randomly selected date
  # shock_map=unique_dates %>%
  #   sample_n(num_shocks) %>% 
  #   dplyr::mutate(e_shock=stats::runif(dplyr::n(), -1, 1 )) %>%
  #   dplyr::full_join(unique_dates, by="extra_t") %>%
  #   dplyr::mutate(e_shock=tidyr::replace_na(e_shock,0)) %>%
  #   dplyr::arrange(extra_t)
  
  #two: shocks propogate until they switch
  # If it's two, then just change the join and uncomment.
  shock_locs <- c(0, sort(sample(1:52, size = num_shocks, replace = F)), 52)
  extra_shocks <- stats::runif(n = num_shocks + 1, min = -1, max = 1)
  shock_seq <- rep(rep(extra_shocks, diff(shock_locs)),
                   length.out = max(date_tib[["time"]]))



  # Map the randomly determined shocks into the dates
  shock_map <- tibble::tibble(
    time = seq_len(max(date_tib[["time"]])),
    e_shock = shock_seq
  )
  
  tib_out= date_tib %>%
    dplyr::inner_join(shock_map, by="time") %>%
    dplyr::group_by(e_shock) %>% 
    dplyr::mutate(
      !!as.name(factor_name) := stats::arima.sim(
        model = ar_model, n = dplyr::n(),
        innov = e_shock + stats::rnorm(dplyr::n(),
                                                   sd = 0.1
        ),
        n.start = 500 ) ) %>%  
    dplyr::ungroup() %>%
    dplyr::select(tidyselect::all_of(c(factor_name, "e_shock")))
  
  return(tib_out)
}

generate_factors <- function(num_factors_inp,
                             num_periods_inp, num_entry_inp,
                             date_start_inp, date_end_inp, freq_inp) {
  # ar model description -- AR 1 with auto correlation and sd inputs
  ar_model <- list(order = c(1, 0, 0), ar = 0.2)
 
  # combine the zero matrix of factors with date indicators
  factor_tib <- generate_time_grid(
    date_start_inp = date_start_inp,
    num_periods_inp = num_periods_inp,
    freq_inp = freq_inp,
    num_entry_inp = num_entry_inp
  ) %>%
    dplyr::mutate(factor1 = time / dplyr::n() +
                    stats::rnorm(dplyr::n(), mean = 0, sd = 0.1))
  
  #Month and Quarter Factors
  base_vec=switch (freq_inp,
    "daily" = c("quarter_num", "month_num", "week_num","day_num"),
    "weekly" = c("quarter_num", "month_num","week_num"),
    "monthly" = c("quarter_num", "month_num")
  )
  
  shocks_tib=purrr::map2_dfc(.x=base_vec,
                          .y=c(glue::glue("factor{2:(length(base_vec)+1)}")),
                          .f=~.BaseFactorHelper(date_tib=factor_tib,
                                            freq_inp=freq_inp, shock_name=.x, 
                                            factor_name=.y, ar_model=ar_model))
  
  factor_tib <- factor_tib %>% dplyr::bind_cols(shocks_tib)
  
  
  if (num_factors_inp > length(base_vec)+1) {
    # Add additional factors beyond the number implied by freq_inp.
    extra_factors_names <- 
      c(glue::glue("factor{(length(base_vec)+2):num_factors_inp}"))
    
    extra_factor_tib <-
      purrr::map_dfc(.x=extra_factors_names,
                      .f=~.ExtraFactorHelper(date_tib=factor_tib,
                                            factor_name=.x, ar_model=ar_model))
    
    factor_tib <- factor_tib %>%
      dplyr::bind_cols(extra_factor_tib)
  }
  browser()
  factor_tib <- factor_tib %>% dplyr::select(-tidyselect::contains("shock"))
  return(factor_tib)
}



# TODO(alexdkellogg): extra factors are essentially weekly noise
#    perhaps make the 4th weekly and the 5th daily (by selecting more number
#    as well as their location among 1-365)
#    Issue is that if this is monthly, then shocks wont be repeated by day...
add_extra_factors <- function(factor_tib, col_in, num_factors_inp,
                              num_periods_inp,
                              ar_model_inp) {
  browser()
  # compute the number shocks and their respective locations
  # Shock number is randomly determined between 1 and 13
  # Shock location is randomly determined as one of the weeks
  num_shocks <- sample(1:13, 1)
  shock_locs <- c(0, sort(sample(1:52, size = num_shocks, replace = F)), 52)
  extra_shocks <- stats::runif(n = num_shocks + 1, min = -1, max = 1)
  shock_seq <- rep(rep(extra_shocks, diff(shock_locs)), 
                   length.out = num_periods_inp)
  
  
  
  # Map the randomly determined shocks into the dates
  shockXwalk <- tibble::tibble(
    time = seq_len(num_periods_inp),
    e_shocks = shock_seq
  )
  factor_tib <- factor_tib %>%
    dplyr::left_join(shockXwalk, by = "time") %>%
    dplyr::mutate(
      !!as.name(col_in) :=
        stats::arima.sim(
          model = ar_model_inp, n = dplyr::n(),
          innov = e_shocks + stats::rnorm(dplyr::n(),
                                          sd = 0.1
          ),
          n.start = 500
        )
    )
  factor_tib <- factor_tib %>% dplyr::select(tidyselect::all_of(col_in))
  return(factor_tib)
}










library(gsynth)
library(panelView)
library(dbarts)
library(ggplot2)
library(dplyr)
ci <- 0.75
data(gsynth)
panelView(Y ~ D, data = simdata,  index = c("id","time"), pre.post = TRUE) 


controls <-
  simdata %>% group_by(id) %>% summarise(D = sum(D)) %>% filter(D == 0) %>% 
  pull(id)

treated <- simdata %>% dplyr::select(id) %>% distinct() %>% 
  filter(!(id %in% controls)) %>% pull(id)

control_df <- simdata %>% dplyr::select(id, time, Y) %>%
  filter(id %in% controls) %>% 
  tidyr::pivot_wider(names_from = id, values_from = Y, names_prefix = "Y_")

#long form df, with each row and id_time combination, and columns Y are the output
#Y_# means the output of control unit #
#Note that for each treated unit we essentially have their 20 periods
# and the 20 periods of every control unit
#IDs only exist for the treated
analisys_df <- simdata %>% dplyr::select(id, time, Y, X1, X2) %>%
  filter(id %in% treated) %>% 
  inner_join(control_df, by = "time") %>% 
  mutate(id = as.factor(id))

#store only the pre periods
pre_data <- analisys_df %>%
  filter(time < 21)
#format for bart
train_data <- dbartsData(pre_data %>% dplyr::select(-Y), pre_data$Y)

nskip <- 20000
nchain <- 8
ndpost <- 4000
keepevery <- (ndpost * nchain)/2000

#fir model to training data (pre period)
bartFit <-
  bart2(
    train_data,
    keepTrees = TRUE,
    combineChains = TRUE,
    n.burn = nskip,
    n.samples = ndpost,
    n.thin = keepevery,
    n.chains = nchain,
    n.threads = nchain,
    verbose = FALSE
  )


df_plot <- tibble(sigma = bartFit$sigma) %>% mutate(draw=1:n()) 

ggplot(data=df_plot, aes(y=sigma,x=draw)) + geom_line()

neff_sigma <- coda::effectiveSize(x = df_plot$sigma)



### Pre intervention

#unique identifier for each entry time combination of TREATED units
idxXwalk <- pre_data %>% 
  dplyr::select(id, time) %>% 
  mutate(idx = 1:n())

#predict seems to output a 100 column df, with column presumably related to id above
#Each row is a draw of this
#so i think one row is one draw, with columns 1:20 being the predicted outcome 
#treated 1, then 21:40 predicted of treated 2, etc
y_hat_draws <- predict(bartFit, pre_data, group.by = pre_data$id) %>%
  as.data.frame() %>%
  mutate(draw = 1:n()) %>%
  tidyr::pivot_longer(names_to = "idx",
                      values_to = "y",
                      -draw) %>%
  mutate(idx = gsub(pattern = "V", replacement = "", x = idx, fixed = T),
         idx = as.numeric(idx)) %>% 
  inner_join(idxXwalk, by = "idx")





y_hat <- y_hat_draws %>%
  group_by(id, time) %>%
  summarise(y_hat = median(y),
            LB = quantile(y, (1-ci)/2),
            UB = quantile(y, 1-(1-ci)/2))

df_plot <- inner_join(pre_data, y_hat, by = c("id", "time"))

ggplot(data = df_plot, aes(x=time, y=Y)) +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2, color = "blue") +
  geom_line() +
  facet_wrap(~id)

### Post intervention

post_data <- analisys_df %>%
  filter(time >= 21) 

idxXwalk <- post_data %>% 
  dplyr::select(id, time) %>% 
  mutate(idx = 1:n())


#now using the post data, same thing
#predict spits out a column for each treated unit*time combo in post t
#this means we have row 1 = draw1, col1=treat1 time 1, col2=treat1 time 2, etc
y_hat_draws_post <- predict(bartFit, post_data) %>%
  as.data.frame() %>%
  mutate(draw = 1:n()) %>%
  tidyr::pivot_longer(names_to = "idx",
                      values_to = "y",
                      -draw) %>%
  mutate(idx = gsub(pattern = "V", replacement = "", x = idx, fixed = T),
         idx = as.numeric(idx)) %>% 
  inner_join(idxXwalk, by = "idx")

y_hat_post <- y_hat_draws_post %>%
  group_by(id, time) %>%
  summarise(y_hat = median(y),
            LB = quantile(y, (1-ci)/2),
            UB = quantile(y, 1-(1-ci)/2))


df_plot_post <- inner_join(post_data, y_hat_post, by = c("id", "time"))
df_plot <- dplyr::bind_rows(df_plot, df_plot_post)

ggplot(data = df_plot, aes(x=time, y=Y)) +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2, color = "blue") +
  geom_line() +
  facet_wrap(~id)

df_ate_pre <- y_hat_draws %>% 
  inner_join(pre_data, by = c("id", "time")) %>% 
  dplyr::select(y_hat = y, Y, draw, id, time) %>% 
  mutate(tau_draw = Y - y_hat) %>% 
  group_by(draw, time) %>% 
  summarise(att_draw = mean(tau_draw)) %>% 
  group_by(time) %>%
  summarise(ATT = median(att_draw),
            LB = quantile(att_draw, (1-ci)/2),
            UB = quantile(att_draw, 1-(1-ci)/2))

df_ate <- y_hat_draws_post %>% 
  inner_join(post_data, by = c("id", "time")) %>% 
  dplyr::select(y_hat = y, Y, draw, id, time) %>% 
  mutate(tau_draw = Y - y_hat) %>% 
  group_by(draw, time) %>% 
  summarise(att_draw = mean(tau_draw)) %>% 
  group_by(time) %>%
  summarise(ATT = median(att_draw),
            LB = quantile(att_draw, (1-ci)/2),
            UB = quantile(att_draw, 1-(1-ci)/2)) %>% 
  bind_rows(df_ate_pre)

ggplot(data = df_ate, aes(x=time, y=ATT)) +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2, color = "blue") +
  geom_line() +
  geom_vline(xintercept = 21, linetype="dashed")




#Ensemble usage
furrr::future_map2(.x=formatted_data[c(1,2)], 
                   .y=preprocess_ensemble_input(list(gs_ens, mc_ens)),
                   .f=~estimate_ensemble(
                     method_names=c("gsynth", "mc"),
                     true_data=.x, pred_list=.y))


estimate_ensemble( method_names=c("gsynth", "mc"),
                   true_data=formatted_data[[1]], pred_list=list(gs_ens[[1]], mc_ens[[1]]),
                   indiv_weights = F)



#check if the selection is occuring on the right components
#specifically, want selection on growth and or size
test_data_raw=do.call(factor_synthetic_dgp,list(date_start="2010-01-01",
                                                 first_treat="2017-07-01",
                                                 date_end="2020-01-01",
                                                 num_entries=1000,
                                                 prop_treated=0.50,
                                                 treat_impact_sd = 0, 
                                                 treat_impact_mean = 0, 
                                                 rho=0.1,
                                                 num_factors=4,
                                                 rescale_y_mean = 2.5e3,
                                                 cov_overlap_scale = 0,
                                                 intercept_scale = 0,
                                                 loading_scale = 0.8,
                                                 seed=42))


#obs: loading scale=0 leads to no difference, roughly, in size or growth
#obs: loading scale=0.95 has large differences in both, particularly growth... however, can  be negative!!
#obs: loading scale=0.75 has a bit of a difference in growth, not size
#

check_selection <- function(input){
  input %>% group_by(treated) %>% filter(time<treatment_period) %>% 
    summarise(avg_load1=mean(loading1),
              avg_load2=mean(loading2),
              avg_load3=mean(loading3),
              mean_inter=mean(intercept)) 
  size=input %>% group_by(treated) %>% filter(time<treatment_period) %>% summarise(avg_size=mean(y0)) 
  growth=input %>% group_by(entry)  %>% filter(time==treatment_period-30 | time==treatment_period) %>% 
    mutate(Growth = (y0 - lag(y0))/lag(y0)) %>% ungroup() %>% group_by(treated) %>%
    summarise(avg_growth=mean(Growth, na.rm=T))
  
  return(list(size, growth))
}
check_selection(test_data_raw)

#Attempt at my own SDID


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
  
  
  
  
  sdid_weights(treat_data = split_treat_data_pre[[1]],
               control_mat = control_matrix,
               treat_time = list_of_treat_times[[1]])
  #gonna need weights from this method
  #Then, need time weights (similar regression except now the outcome var
  # is the first Treatment Period for each control unit, and the we predict using
  #all the past values for that particular unit?)
  list_inputed_y=furrr::future_map2(.x=split_treat_data_pre,
                                    .y=list_of_treat_times,
                                    .f=~sdid_weights(treat_data = .x,
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




sdid_weights <- function(treat_data, control_mat, treat_time){
  #define the control data matrix
  browser()
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
gfoo_match <- function(ts_tomatch, ts_pool, treatperiod, 
                                 num.neighbours = 8, neighbour.weighting = F,
                                 scaling.option = c("centered", "mean", "none")) {
  
  neighbour.method <- match.arg(neighbour.method)
  scaling.option <- match.arg(scaling.option)
  
  # Scale the histories if needed.
  yx <- cbind(ts_tomatch, ts_pool)
  if (scaling.option == "centered") {
    yx <- scale(yx)
  } else if (scaling.option == "mean") {
    yx <- t(t(yx) / colSums(yx)) * nrow(yx)
  }
  yx[, which(is.na(yx[1, ]))] <- 0
  
  # Set the weights for distance calculation.
  if (neighbour.weighting) {
    #increasing weights up to Treatment period, then decreasing again
    weights <- c(seq_len(treatperiod), 
                 seq(treatperiod, by=-1, length.out = nrow(yx)-treatperiod))
  } else {
    weights <- rep(1, nrow(yx))
  }
  distance <- colSums(((yx[, seq(from = 2, to = ncol(yx), by = 1)]
                        - yx[, 1]) * weights)^2)
  
  # id <- order(distance)[seq_len(num.neighbours)]
  # Modification: We choose 1.5 times as many neighbours but only use those
  # with closest volume
  id.init <- order(distance)[seq_len(floor(num.neighbours * 1.5))]
  
  xx <- ts_pool[, id.init]
  yxx <- cbind(ts_tomatch, xx)
  distance <- colSums((yxx[, seq(from = 2, to = ncol(yxx), by = 1)]
                       - yxx[, 1])^2)
  id <- id.init[order(distance)[seq_len(num.neighbours)]]
  return(id[1])
}


gfoo_match_helper <- function(treated_block,
                                         control_block,
                                         id_var, time_var,
                                         treat_period, outcome_var) {
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
  placebo_treat_period <- treated_block %>%
    dplyr::distinct(!!as.name(id_var), .keep_all=T) %>%
    dplyr::select(tidyselect::all_of(c(id_var,treat_period))) %>%
    dplyr::rename(Treatment_Unit=id_var)
  treat_pool=treated_block %>% 
    dplyr::distinct(!!as.name(id_var)) %>% 
    dplyr::pull() %>% 
    sample()
  for (i in seq_along(treat_pool)) {
    #Find nearest euclidean match among unmatched controls for 
    #each treated observation
    control_mat=control_block %>%
      dplyr::filter(!!as.name(id_var) %in% setdiff(!!as.name(id_var), 
                                                   already_matched)) %>%
      tidyr::spread(!!as.name(time_var), !!as.name(outcome_var) ) %>%
      dplyr::select(-!!as.name(id_var)) %>%
      as.matrix() %>%
      t()
    
    treat_unit=treated_block %>%
      dplyr::filter(!!as.name(id_var)==treat_pool[i]) %>%
      dplyr::pull(!!as.name(outcome_var))
    
    treat_time=treated_block %>%
      dplyr::filter(!!as.name(id_var)==treat_pool[i]) %>%
      dplyr::distinct(!!as.name(treat_period)) %>%
      dplyr::pull()
    
    temp_match <- gfoo_match(
      ts_tomatch =treat_unit,
      ts_pool=control_mat,
      treatperiod = treat_time,...)
    #Update the vector of already matched donors
    already_matched[i] <- control_block %>%
      dplyr::filter(!!as.name(id_var) %in% 
                      setdiff(!!as.name(id_var), already_matched)) %>%
      dplyr::distinct(!!as.name(id_var)) %>%
      dplyr::slice(temp_match) %>%
      dplyr::pull()
  }
  # Store the resulting vectors in a tibble for output
  df_toreturn <- tibble::tibble(temp_id = already_matched,
                                Treatment_Unit = treat_pool) %>%
    dplyr::inner_join(placebo_treat_period, by="Treatment_Unit") %>% 
    dplyr::rename(!!as.name(id_var) := temp_id)
  
  return(df_toreturn)
}


#attempting to match based on TS features and compare the performance (seems worse)
#TODO(alexdkellogg): account for multiple treatments, treatment end dates
create_placebo_df <- function(data_full, id_var = "entry",
                              time_var = "period", 
                              treat_indicator = "treatperiod_0",
                              outcome_var = "target", 
                              counterfac_var = "counter_factual",
                              match_type="old") {
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
  cd <- data_full %>% 
    dplyr::filter(!!as.name(id_var) %in% ct_entries) %>%
    dplyr::select(tidyselect::all_of(c(time_var, id_var, outcome_var )))
  
  
  # Store treated data in a tibble
  treated_to_match <- data_full %>%
    dplyr::filter(!!as.name(id_var) %in% tr_entries) %>%
    dplyr::group_by(!!as.name(id_var)) %>%
    dplyr::mutate(
      Treatment_Period =
        length(!!as.name(treat_indicator)) - sum(!!as.name(treat_indicator)) + 1
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(tidyselect::all_of(c(time_var, id_var, outcome_var,
                                       "Treatment_Period" )))
  
  if(match_type=="gfoo"){
    # for pivoting, potential issues arise if we have several time varying 
    # covariates -- added as TODO.
    # we'd have to take the values_from each of them, 
    #and for any constant args we'd presumably have to add them to id_cols
    #Pivot the data wide, so each row has all the time series data for 
    #a given unit.

    # Determine control units to assign to placebo treatment group via matching.
    matched_placebo_df_temp <- 
      matching_without_replacement(treated_to_match, cd, id_var=id_var,
                                   treat_period = "Treatment_Period",
                                   outcome_var = outcome_var,
                                   time_var = time_var)
  }
  if(match_type=="feature"){
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
    
    matched_placebo_df_temp<- matching_without_replacement(treated_ts_features, 
                                                           control_ts_features, 
                                                           id_var,
                                                           "Treatment_Period")
  }
  
  if(match_type=="old"){
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






