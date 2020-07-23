library(dplyr)

#TODO(alexdkellogg): allow dgp to be covariates only or factors only
#TODO(alexdkellogg): allow num_periods directly rather than dates
#TODO(alexdkellogg): introduce coeff vector and have xBeta go into Y
factor_synthetic_dgp<-function(num_entries=2000,
                               date_start="2017-01-01",
                               first_treat="2018-07-01",
                               date_end="2020-01-01",
                               freq=c("daily", "weekly", "monthly"),
                               prop_treated=0.4,
                               treat_impact_mean=0.1,
                               treat_impact_sd=0.1,
                               treat_decay_mean=0.7,
                               treat_decay_sd=1,
                               selection=c("random", "observables", 
                                           "unobservables"),
                               rho=0.9, 
                               rho_scale=0.2, 
                               rho_shift=0, 
                               cov_overlap_scale=0, 
                               num_factors=3,
                               loading_scale=0,
                               intercept_scale=0,
                               conditional_impact_het=-0,
                               rescale_y_mean=2.5e4,
                               seed=19,
                               log_output=T){ 
  #rescale_y_max=5e8
  #=0,
  
  #Generates tibble of long form panel data using a factor-augmented AR 1
  #Each row is time period x unit unique combination.
  
  #Args
  #num_entries: number of units to be generated
  #date_start: string "yyyy-mm-dd" input for the first simulated date
  #first_treat: string "yyyy-mm-dd" input first treatment period.
  #date_end:: string "yyy-mm-dd" input for the last simulated date
  #freq: string indicating time unit, either "daily", "weekly", "monthly"
  #prop_treated: proportion of entries that should receive treatment
  #treat_impact_mean: initial period treatment impact mean, drawn from truncated
  #   normal distribution centered here. The end points of the dist are [0,0.25]
  #   by default, but shift to [a,b] where b=mean+0.25 if the mean is larger
  #   than 0.25 (max of 1)  or a=mean-0.25 if mean is below 0.
  #treat_impact_sd: standard deviation of the truncated normal mean impact. 
  #treat_decay_mean: initial period treatment decay mean, drawn from truncated
  #   normal distribution centered here. The end points of the dist are [0,0.9]
  #   by default, but shift to [0,1+eps] if mean>0.9. This allows for units to
  #   have no decay (value of 1). Propagates as mean**(post_treat_period).
  #treat_impact_sd: standard deviation of the truncated normal decay factor 
  #selection= string in "random", "observables", "unobservables" dictating 
  #   treatment assignment mechanism.
  #rho: mean of truncated normal distribution of the autocorrelation of outcome,
  #   with bounds [0, 0.995].
  #rho_scale: standard deviation of truncated normal for the autocorrelation.
  #rho_shift: multiplier on the mean rho for control units, 
  #   overall mean stay below 1 (rho*rho_shift<1).
  #cov_overlap_scale: (-1,1) shifts distribution of covariates for a fraction 
  #   (prop_treated) of x variables. A value of 1 shifts the distribution for 
  #   treatment up on all x's, whereas -1 shifts up the distribution for donors.
  #num_factors: number (3+) of time-varying, unobserved factors to simulate
  #loading_scale: (-1,1) shift in factor distribution, -1 shifts loadings up
  #   for control units, positive values shift loadings distribution up for 
  #   treated units.
  #intercept_scale: (-1,1) shifts the mean of a truncated normal distribution
  #   for control unit intercepts: >0 shifts the mean down (treatment is larger)
  #   and <0 shift control mean above the treatment mean.
  #conditional_impact_het: constant added to the treatment impact for top 25%
  #   and subtracted from bottom 25% of counterfactual y at time 1.
  #log_output: boolean to determine whether output is log scale or exp
  #seed: random number seed
  
  
  
  #Output
  #Long form tibble with columns for observed, potential treated, and
  #  potential untreated outcomes, treatment time, and
  #  the relevant x variables (loadings, intercept, observables).
  
  
  
  
  #Match the arguments to valid entries
  set.seed(seed)
  if (missing(selection)) {
    selection <- "random"
  } else{
    selection <- match.arg(selection)
  }
  
  if (missing(freq)) {
    freq <- "monthly"
  } else{
    freq <- match.arg(freq)
  }
  #require 10% min in both treat and control
  stopifnot(prop_treated >= 0.1 & prop_treated <= 0.9)
  #require cov_overlap_scale to be between -1 (shift treat down) and 1
  stopifnot(cov_overlap_scale <= 1 & cov_overlap_scale >= -1)
  #require 3 factors of more
  stopifnot(num_factors>2)
  #require less than 100% TE
  stopifnot(treat_impact_mean<1)
  
  #given the dates and frequency, identify total number of periods
  num_periods=time_interval_calculator(start_t=date_start, end_t=date_end,
                                  freq_t=freq)
  #Store how long after start date the treatment begins
  treat_start_int=time_interval_calculator(start_t=date_start, 
                                           end_t=first_treat, freq_t=freq)
  #Stop if there are too few pre treat periods
  stopifnot(treat_start_int < 0.8*num_periods)
  
  #
  #assign covariates and treatment given selection and overlap
  synth_data_unit=unit_level_simulation(n_inp=num_entries,
                               type_inp=selection, 
                               cov_overlap_inp=cov_overlap_scale,
                               loading_scale_inp=loading_scale,
                               num_factors_inp=num_factors,
                               int_scale_inp=intercept_scale,
                               rho_inp=rho,
                               rho_scale_inp=rho_scale,
                               rho_shift_inp=rho_shift,
                               prop_treated_inp=prop_treated,
                               treat_start=treat_start_int,
                               num_periods_inp=num_periods,
                               impact_mean_inp=treat_impact_mean,
                               impact_sd_inp=treat_impact_sd,
                               decay_mean_inp=treat_decay_mean,
                               decay_sd_inp=treat_decay_sd)

  #next, work on time varying characteristics
  synth_data_factors=generate_factors(num_factors_inp=num_factors,
                                      num_periods_inp=num_periods, 
                                      num_entry_inp=num_entries,
                                      date_start_inp=date_start,
                                      date_end_inp=date_end,
                                      freq_inp=freq)
  
  #Next, we generate both potential outcomes
  #Option for factor only (y=factor*loadings), Random Effect only
  # (y=xB), or both.
  #Then, generate the treatment impact
  unit_time_grid <-tidyr::expand_grid(entry = seq_len(num_entries),
                                      time = seq_len(num_periods)) 
  
  synth_data_full=unit_time_grid %>% 
    dplyr::left_join(synth_data_unit,by=c("entry")) %>%
    dplyr::left_join(synth_data_factors,by=c("time"))
  
  #generate the counterfactual outcomes
  synth_data_full=generate_counterfactual(synth_data_full,
                                          num_periods_inp=num_periods,
                                          rescale_y=rescale_y_mean)
  
  
  #generate the per period impact (taking acocunt of decay)
  synth_data_full=generate_treat_impact(data_inp=synth_data_full,
                                        cond_impact_inp=
                                          conditional_impact_het) %>%
    dplyr::mutate(y=exp(y)*(1- log_output)+log_output*y,
                  y0=exp(y0)*(1- log_output)+log_output*y0,
                  y1=exp(y1)*(1- log_output)+log_output*y1) 
  
  return(synth_data_full)
  
  
  
}


compute_factor_loadings <- function(data_inp){
  #create a list of factor loadings split by individual
  loadings_vec_list=data_inp %>%
    dplyr::select(entry, tidyselect::contains("loading")) %>%
    dplyr::distinct(entry, .keep_all=T) %>%
    dplyr::group_split(entry, .keep=F) %>%
    lapply(as.matrix)
  #create a list of factor matrices split by individual
  factor_mat_list=data_inp %>%
    dplyr::select(entry, tidyselect::contains("factor")) %>%
    dplyr::group_split(entry, .keep=F) %>% 
    lapply(as.matrix)
  #Compute the matrix product of loadings and factors for each individual
  return(furrr::future_map2(.x=loadings_vec_list,.y=factor_mat_list,
                                         .f=~(.y)%*%t(.x) ) %>%
    unlist())
}

#TODO(alexdkellogg): check with AP about 0 noise AR for y
compute_outcome_process<-function(data_inp,num_periods_inp){
  #create a list of AR models, with individual specific noise and autocorr
  ar_param_inp=data_inp %>%
    dplyr::select(entry, autocorr) %>%
    dplyr::distinct(entry, .keep_all=T) %>%
    dplyr::group_split(entry, .keep=F) %>%
    lapply(function(x){ list(order=c(1,0,0), ar=x[[1]])} )
  
  
  innov_list=data_inp %>%
    dplyr::select(entry, tidyselect::contains("loading")) %>%
    dplyr::group_split(entry, .keep=F)   %>%
    lapply(as.matrix)
  
  
  return(furrr::future_map2(.x=ar_param_inp, .y=innov_list,
                               .f=~stats::arima.sim(model=.x,
                                                    n=num_periods_inp,
                                                    n.start = 500,
                                                    innov = .y)) %>%
    unlist())
  
}


generate_counterfactual<-function(data_inp,num_periods_inp, rescale_y){
  
  outcome_series=data_inp %>%
    dplyr::select(time, entry, intercept, autocorr, 
                  tidyselect::matches("loading|factor"))
  
  computed_factor_vec=compute_factor_loadings(outcome_series)
  
  outcome_series=outcome_series %>%
    dplyr::select(-tidyselect::matches("loading|factor")) %>%
    dplyr::mutate(factor_loadings=computed_factor_vec)
  
  outcome_ar=compute_outcome_process(outcome_series,num_periods_inp)
  
  xvars=data_inp %>% 
    dplyr::distinct(entry, .keep_all=T) %>%
    dplyr::select(tidyselect::matches("obs")) %>%
    as.matrix() 
  
  x_beta=xvars %*% stats::rnorm(n=ncol(xvars), sd=0.5) %>%
    as.vector()
  
  x_betaXwalk=tibble::tibble(entry=seq_len(length(x_beta)),
                             xbeta=x_beta) 
  outcome_series=outcome_series %>%
    dplyr::left_join(x_betaXwalk, by="entry") %>%
    dplyr::mutate(factor_loadings=computed_factor_vec,
                  indiv_component=intercept+xbeta,
                  y0=outcome_ar+indiv_component+stats::rnorm(dplyr::n(), sd=0.5))
  
  outcome_series_rescaled=outcome_series %>%
    dplyr::mutate(de_sd_y0=(y0/sd(y0)),
                  rescaled_y=de_sd_y0-(mean(de_sd_y0)-log(rescale_y)),
                  y0=1.125*rescaled_y) %>%
    dplyr::select(time, entry, factor_loadings, y0)
  
  return(data_inp %>% 
           dplyr::inner_join(outcome_series_rescaled, by=c("time", "entry")) %>%
           dplyr::select(time, entry, treated, treatment_period,y0,
                         factor_loadings, dplyr::everything())
         )
  

}

#TODO(alexdkellogg): add conditional impact
generate_treat_impact<-function(data_inp=synth_data_full, 
                                cond_impact_inp){
  #determine which observations get a conditional treatment boost
  condXwalk=data_inp %>%
    dplyr::filter(time==1) %>%
    dplyr::select(entry, y0) %>%
    dplyr::mutate(
      cond_treat_impact=dplyr::case_when(
        y0>quantile(y0,0.75)~cond_impact_inp,
        y0<quantile(y0,0.25)~-cond_impact_inp,
        TRUE~0 )) %>% 
    dplyr::select(-y0)

    
  #Define the treatment propogation for each unit and time combo
  data_inp=data_inp %>%
    dplyr::left_join(condXwalk, by="entry") %>%
    dplyr::group_by(entry) %>%
    dplyr::mutate(
      post_treat_t=time-treatment_period,
      decay_t=dplyr::case_when(
        post_treat_t<0~0,
        post_treat_t>=0~treat_decay**post_treat_t),
      impact_t=decay_t*(treat_impact+cond_treat_impact),
      ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(treat_decay, treat_decay, decay_t, cond_treat_impact)) 
  
  #Add the treatment impact to create y1
  return(data_inp %>% 
    dplyr::mutate(y1=impact_t+y0,
                  y=treated*y1+(1-treated)*y0) %>%
      dplyr::select(time, entry, treated, 
                    post_treat_t, treatment_period, impact_t,
                    y, y0,y1,
                    dplyr::everything()))
}


generate_time_grid <- function(date_start_inp,num_periods_inp, 
                               freq_inp, num_entry_inp){
  period_dates=switch(freq_inp,
                      "daily"=format(lubridate::ymd(date_start_inp)+lubridate::days(0:(num_periods_inp-1))),
                      "weekly"=format(lubridate::ymd(date_start_inp)+lubridate::weeks(0:(num_periods_inp-1))),
                      "monthly"=format(lubridate::ymd(date_start_inp)+months(0:(num_periods_inp-1)))
  )
  
  #Identidy the relevant components of the date (day/week/etc)
  date_info_tib=switch(freq_inp,
                        "daily" = tibble::tibble(time=seq_len(num_periods_inp),
                                         date_t=period_dates,
                                         day_num=lubridate::day(period_dates), 
                                         week_num=lubridate::week(period_dates),
                                         month_num=lubridate::month(period_dates),
                                         quarter_num=lubridate::quarter(period_dates),
                                         year_num=lubridate::year(period_dates)),
                        "weekly"=tibble::tibble(time=seq_len(num_periods_inp),
                                        date_t=period_dates,
                                        week_num=lubridate::week(period_dates),
                                        month_num=lubridate::month(period_dates),
                                        quarter_num=lubridate::quarter(period_dates),
                                        year_num=lubridate::year(period_dates)),
                        "monthly"=tibble::tibble(time=seq_len(num_periods_inp),
                                         date_t=period_dates,
                                         month_num=lubridate::month(period_dates),
                                         quarter_num=lubridate::quarter(period_dates),
                                         year_num=lubridate::year(period_dates))
                        )
  
  return( date_info_tib )
}


#plan: factors are Txnum_fac dimensional, and must be multiplied by loadings
#can group by treated and call generate factors with varying rho using
#rho_het_inp, which multiplies rho for the counterfactuals

generate_factors<-function(num_factors_inp, 
                           num_periods_inp, num_entry_inp,
                           date_start_inp, date_end_inp, freq_inp){
  
  #generate factors from an AR 1 process
  factor_mat <-matrix(0, nrow=num_periods_inp,ncol = num_factors_inp )
  colnames(factor_mat) <- glue::glue("factor{1:num_factors_inp}")
  #ar model description -- AR 1 with auto correlation and sd inputs
  ar_model=list(order=c(1,0,0), ar=0.2)
  #TODO(alexdkellogg): check with AP if shocks trend over time, assumed fixed
  quarter_effects=tibble::tibble(q_shock=stats::runif(4, -1, 1), 
                                 quarter_num=seq_len(4))
  month_effects=tibble::tibble(m_shock=stats::runif(12, -1, 1), 
                               month_num=seq_len(12))
    

  #combine the zero matrix of factors with date indicators
  factor_tib=generate_time_grid(date_start_inp=date_start_inp,
                                num_periods_inp=num_periods_inp, 
                                freq_inp=freq_inp,
                                num_entry_inp=num_entry_inp) %>%
    dplyr::bind_cols(tibble::as_tibble(factor_mat)) %>%
    dplyr::inner_join(quarter_effects, by="quarter_num") %>%
    dplyr::inner_join(month_effects, by="month_num")
  
  
  #add first factor, period/total + noise
  factor_tib=factor_tib %>%
    dplyr::mutate(factor1=time/dplyr::n()+
                    stats::rnorm(dplyr::n(),mean=0,sd=0.1)) %>%
    dplyr::group_by(quarter_num) %>%
    dplyr::mutate(
      factor2=stats::arima.sim(model=ar_model, n=dplyr::n(),
                               innov = q_shock+stats::rnorm(dplyr::n(),
                                                            sd=0.1),
                               n.start = 500)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(month_num) %>%
    dplyr::mutate(
      factor3=stats::arima.sim(model=ar_model, n=dplyr::n(),
                               innov = m_shock+stats::rnorm(dplyr::n(),
                                                            sd=0.1),
                               n.start = 500)) %>%
    dplyr::ungroup()
  if(num_factors_inp>3){
    #this process works for one factor. Want to lapply to all columns >4
    extra_factors_names=setdiff(
      names(factor_tib %>% dplyr::select(tidyselect::contains("factor"))),
      c(glue::glue("factor{1:3}"))
      )
    extra_factor_tib=purrr::map(.x=extra_factors_names,
                                       .f=~add_extra_factors(factor_tib=factor_tib,
                                 num_factors_inp=num_factors_inp,
                                 num_periods_inp=num_periods_inp,
                                 ar_model_inp=ar_model, col_in = .x)) %>%
      dplyr::bind_cols() 
    
    factor_tib=factor_tib %>%
      dplyr::select(-tidyselect::all_of(extra_factors_names)) %>%
      dplyr::bind_cols(extra_factor_tib)
  }
  
  return(factor_tib %>% dplyr::select(-tidyselect::contains("shock")))
                  
}


add_extra_factors<-function(factor_tib,col_in,num_factors_inp,
                            num_periods_inp,
                             ar_model_inp){
  
  #compute the number shocks and their respective locations
  num_shocks=sample(1:13,1)
  shock_locs=c(0,sort(sample(1:52, size =num_shocks, replace = F )),52)
  extra_shocks=stats::runif(n=num_shocks+1, min=-1, max=1)
  shock_seq=rep(rep(extra_shocks, diff(shock_locs)), length.out=num_periods_inp)
    
  
  
  
  shockXwalk=tibble::tibble(time=seq_len(num_periods_inp),
                            e_shocks=shock_seq)
  factor_tib=factor_tib %>% 
    dplyr::left_join(shockXwalk, by="time") %>%
    dplyr::mutate(
    !!as.name(col_in):=
      stats::arima.sim(model=ar_model_inp, n=dplyr::n(),
                       innov = e_shocks+stats::rnorm(dplyr::n(),
                                                     sd=0.1), 
                       n.start = 500)) 
  
  return(factor_tib %>% dplyr::select(tidyselect::all_of(col_in)))
  
}


#will need to add individual specific loadings here
#with heterogeneity given by loading_scale
unit_level_simulation <- function(n_inp,
                                  type_inp, 
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
                                  decay_sd_inp){
  #Gather data that is time invariant
  #First, call to assign treatment, which generates covariates
  unit_level_tib=assign_treat(n_inp=n_inp,
               type_inp=type_inp, 
               cov_overlap_inp=cov_overlap_inp,
               loading_scale_inp=loading_scale_inp,
               num_factors_inp=num_factors_inp,
               rho_inp=rho_inp,
               rho_scale_inp=rho_scale_inp,
               rho_shift_inp=rho_shift_inp,
               int_scale_inp= int_scale_inp,
               prop_treated_inp=prop_treated_inp)
  
  #Merge in the treatment period assignment
  unit_level_tib=assign_treat_time(unit_level_tib, treat_start, num_periods_inp)
  #allow wiggle room for the TE without getting NA
  impact_ub=ifelse(impact_mean_inp<0.25,0.25, 
                   min(1, impact_mean_inp+(0.25)) )
  impact_lb=ifelse(impact_mean_inp<=0,impact_mean_inp-0.25, 0)
  decay_ub=ifelse(decay_mean_inp>0.9,1.000000001, 0.9)
  #assign treatment effect impact and decay parameters per unit
  unit_level_tib=unit_level_tib %>% 
    dplyr::mutate(
      treat_impact=truncnorm::rtruncnorm(n=dplyr::n(), a=impact_lb,b=impact_ub, 
                                         mean=impact_mean_inp, sd=impact_sd_inp),
      treat_decay=truncnorm::rtruncnorm(n=dplyr::n(), a=0,b=decay_ub, 
                                        mean=decay_mean_inp, sd=decay_sd_inp))
  
  return(unit_level_tib)
}








time_interval_calculator <- function(start_t, end_t, freq_t){
  #Computes the amount of time (in units defined by fred_t)
  #between two input dates
  
  #Args:
  # start_t: date string (yyyy-mm-dd) indicating start of interval
  # end_t: date string (yyyy-mm-dd) indicating end of interval
  # freq_t: string, one of "daily", "weekly", "monthly", "yearly"
  
  return(switch(freq_t,
                "daily"=ceiling(
                  lubridate::interval(lubridate::ymd(start_t),
                                      lubridate::ymd(end_t))/
                    lubridate::days(1)+1),
                "weekly"=ceiling(
                  lubridate::interval(lubridate::ymd(start_t),
                                      lubridate::ymd(end_t))/
                    lubridate::weeks(1)+1),
                "monthly"=ceiling(
                  lubridate::interval(lubridate::ymd(start_t),
                                      lubridate::ymd(end_t))/
                    months(1)+1),
                "yearly"=ceiling(
                  lubridate::interval(lubridate::ymd(start_t),
                                      lubridate::ymd(end_t))/
                    lubridate::years(1)+1))
  )
}



#TODO(alexdkellogg): discuss with AP difference between 
#   overlap and selection (do we need both?)
assign_treat<-function( n_inp, type_inp, cov_overlap_inp,
                        loading_scale_inp,rho_inp=rho_inp, 
                        rho_scale_inp,rho_shift_inp, num_factors_inp,
                        int_scale_inp, prop_treated_inp){
  #Creates a tibble with treatment assignment and time constant covariates
  #Covariate distribution can differ by treatment depending on cov_overlap_inp,
  #and selection into treatment is modelled.
  
  #generate covariates with varying levels of overlap
  covariate_tib=gen_covariates(n_inp=n_inp, cov_overlap_inp=cov_overlap_inp,
                               frac_shifted=prop_treated_inp)
  
  if (type_inp == "random") {
    #Randomly assign treatment to prop_treated_inp
    treat_covar_tib=tibble::tibble(
      entry = seq_len(n_inp),
      treated = as.numeric(dplyr::percent_rank(
        runif(n = n_inp,min = 0,max = 1)
      )>1-prop_treated_inp)) %>% 
      inner_join(covariate_tib, by="entry")
    
    unit_tib=generate_loadings(treat_tib_inp=treat_covar_tib,
                                     loading_scale_inp=loading_scale_inp,
                               num_factors_inp=num_factors_inp,
                                     int_scale_inp=int_scale_inp,
                               rho_inp=rho_inp,
                                     rho_scale_inp=rho_scale_inp,
                               rho_shift_inp = rho_shift_inp)
    return(unit_tib)
  } 
  else{ 
    #Assign treatment based on either observable or unobservable xs
    relevant_xs=ifelse(type_inp=="observables", "obs", "unobs")
    #First, create covariates with 
    z <- covariate_tib %>% 
      dplyr::select(tidyr::starts_with(relevant_xs)) %>%
      as.matrix()

    # combine the variables for each observation to form a predicted score
    # rescale with logistic function to map into the probability space (0,1)
    prob_treat= tibble::tibble(score=
                                 z %*% stats::rnorm(n = ncol(z),mean = 0,sd = 1)) %>%
      dplyr::mutate(prop_score = 1 / (1 + exp(score)))
    
    #Ad random noise to the score and take the top fraction as treated
    treat_covar_tib= prob_treat %>%
      dplyr::mutate( u = runif(n = n_inp, min = 0, max = 1),
                     p_new = prop_score+u,
                     treated=as.numeric(dplyr::percent_rank(p_new)>
                                          1-prop_treated_inp),
                     entry = seq_len(n_inp)) %>% 
      dplyr::select(entry, treated) %>%
      inner_join(covariate_tib, by="entry")
    
    unit_tib=generate_loadings(treat_tib_inp=treat_covar_tib,
                               loading_scale_inp=loading_scale_inp,
                               num_factors_inp=num_factors_inp,
                               int_scale_inp=int_scale_inp,
                               rho_inp=rho_inp,
                               rho_scale_inp=rho_scale_inp,
                               rho_shift_inp = rho_shift_inp)
    
    return(unit_tib)
  }
}


generate_loadings <- function(treat_tib_inp, loading_scale_inp,num_factors_inp,
                              int_scale_inp,rho_inp, rho_scale_inp, rho_shift_inp){
  
  loadings_mat=matrix(0, nrow=nrow(treat_tib_inp), ncol=num_factors_inp)
  colnames(loadings_mat)=glue::glue("loading{1:num_factors_inp}")
  
  
  tib_pre_loadings=treat_tib_inp %>% 
    dplyr::group_by(treated) %>%
    dplyr::mutate(
      #intercept=stats::rexp(dplyr::n(), 1+(1-treated)*int_scale_inp),
      intercept=truncnorm::rtruncnorm(n=dplyr::n(), a=5,b=10,
                                      mean=7.5-(1-treated)*int_scale_inp,
                                      sd=0.7),
      #intercept=stats::runif(dplyr::n(), min=5,max=8-(1-treated)*int_scale_inp),
      # autocorr=stats::runif(dplyr::n(),ifelse(treated,rho_inp,
      #                                      rho_inp*rho_scale_inp),
      #                       max=0.95),
      autocorr=truncnorm::rtruncnorm(n=dplyr::n(), a=0,b=0.995,
                                     mean=ifelse(treated,rho_inp,
                                                 rho_inp*rho_shift_inp),
                                     sd=rho_scale_inp)) %>%
    dplyr::bind_cols(tibble::as_tibble(loadings_mat)) %>%
    dplyr::ungroup()
  return(
    tib_pre_loadings %>%
      dplyr::group_by(treated) %>%
      dplyr::mutate(
        dplyr::across(tidyselect::starts_with("loading"),
        .fns=list(load=~stats::rbeta(dplyr::n(), 2 -(1-treated)*(loading_scale_inp) ,
                     2 - (treated)*(loading_scale_inp)) ),
        .names ="{col}") 
      ) %>% dplyr::ungroup()
    
  )
    
  
}

#TODO(alexdkellogg): Generalize to allow for user inputted number of xs
gen_covariates<-function(n_inp, cov_overlap_inp, frac_shifted){
  #Generate several covariates, both observed and unobserved.
  #Shift means according to overlap_inp and frac_shifted 
  
  #first, create a set of correlated variables, mean zero, var/covar Sigma
  Sigma <- matrix(c(16 ,  4,  -4.8,
                    4  , 25,    9,
                    -4.8,  9,    9),3,3)
  corr_xs= MASS::mvrnorm(n = n_inp, rep(0, 3), Sigma)
  colnames(corr_xs)=c("obs_mvnorm1", "unobs_mvnorm", "obs_mvnorm2")
  
  return( 
    x_tib=tibble::tibble(
      entry = seq_len(n_inp),
      to_shift = as.numeric(dplyr::percent_rank(
        runif(n = n_inp,min = 0,max = 1)
                         )>1-frac_shifted)) %>%
    cbind(corr_xs) %>% 
    dplyr::group_by(to_shift) %>%
    dplyr::mutate(
      obs_beta=stats::rbeta(n = dplyr::n(),
                         shape1 = 6,
                         shape2 = 5+3*(to_shift*-(cov_overlap_inp))),
      obs_binom=stats::rbinom(n = dplyr::n(),
                         size = 15,
                         prob = 0.5+0.25*(to_shift*(cov_overlap_inp))),
      obs_norm=stats::rnorm(n = dplyr::n(),
                          mean = 0+5*(to_shift*(cov_overlap_inp)),
                          sd = 10),
      unobs_beta=stats::rbeta(n = dplyr::n(),
                         shape1 = 6,
                         shape2 = 5+3*(to_shift*-(cov_overlap_inp))),
      ubobs_binom=stats::rbinom(n = dplyr::n(),
                          size = 15,
                          prob = 0.5+0.25*(to_shift*(cov_overlap_inp))),
      ubobs_norm=stats::rnorm(n = dplyr::n(),
                         mean = 0+5*(to_shift*(cov_overlap_inp)),
                         sd = 10) ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-to_shift)
  )
  
}

assign_treat_time<-function(treat_tib_inp, treat_start, num_periods_inp){
  #Geometric distribution parameter selected based on time
  #Goal is to have most treatment assigned before end date
  geom_prob=dplyr::case_when(num_periods_inp-treat_start<=20~0.25,
                             num_periods_inp-treat_start<=30~0.2,
                             num_periods_inp-treat_start<=40~0.15,
                             TRUE~0.1 )
  #draw treatment period as first date plus geometric random variable
  return(treat_tib_inp %>% 
           dplyr::mutate(
             treatment_period=treat_start+
               stats::rgeom(dplyr::n(), geom_prob)) %>%
           dplyr::group_by(entry) %>%
           dplyr::mutate(treatment_period=
                           min(treatment_period,num_periods_inp)) %>%
           dplyr::ungroup())
  
}



format_for_est<-function(data_ouput){
  #For estimation, we want a dummy for the treated units in treated times
  reformatted_data=data_ouput %>%
    dplyr::mutate(treatperiod_0=(treated==1)*(post_treat_t >=0) ) %>%
    dplyr::rename(potential_post_t=post_treat_t,
                  potential_treatment_period=treatment_period,
                  potential_impact_t=impact_t,
                  target=y,
                  counter_factual=y0,
                  period=time) %>%
    dplyr::select(period, entry, treated, treatperiod_0 , target,
                  counter_factual)
  
  return(reformatted_data)
  
}

