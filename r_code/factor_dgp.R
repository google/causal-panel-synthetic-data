


factor_synthetic_dgp<-function(num_periods=36,
                               num_entries=2000,
                               num_factors=3,
                               date_start="2017-01-01",
                               first_treat="2018-07-01",
                               date_end="2020-01-01",
                               freq=c("daily", "weekly", "monthly", "yearly"),
                               prop_treated=0.4,
                               rho=0.9,
                               rho_scale=1,
                               overlap_scale=0,
                               treat_impact=0.1,
                               treat_decay=0.7,
                               treat_decay_scale=1,
                               selection=c("random", "observables", 
                                           "unobservables"),
                               seed=19){
  
  
  
  
  
  
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
  #require overlap_scale to be between -1 (shift treat down) and 1
  stopifnot(overlap_scale <= 1 & overlap_scale >= -1)
  
  #given the dates and frequency, identify total number of periods
  N_time=time_interval_calculator(start_t=date_start, end_t=date_end,
                                  freq_t=freq)
  #Store how long after start date the treatment begins
  treat_start_int=time_interval_calculator(start_t=date_start, 
                                           end_t=first_treat, freq_t=freq)
  #Stop if there are too few pre treat periods
  stopifnot(treat_start_int < 0.8*N_time)
  
  #assign covariates and treatment given selection and overlap
  synth_data_unit=assign_treat(N_inp=num_entries,
                               type_inp=selection, 
                               overlap_inp=overlap_scale,
                               prop_treated_inp=prop_treated)
  
  #next, work on time varying characteristics
  return(synth_data_unit)
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




assign_treat<-function( N_inp, type_inp, overlap_inp,prop_treated_inp){
  #Creates a tibble with treatment assignment and time constant covariates
  #Covariate distribution can differ by treatment depending on overlap_inp,
  #and selection into treatment is modelled.
  
  #generate covariates with varying levels of overlap
  covariate_tib=gen_covariates(N_inp=N_inp, overlap_inp=overlap_inp,
                               frac_shifted=prop_treated_inp)
  if (type_inp == "random") {
    #Randomly assign treatment to prop_treated_inp
    return(tibble::tibble(
      entry = seq_len(N_inp),
      treated = as.numeric(dplyr::percent_rank(
        runif(n = N_inp,min = 0,max = 1)
        )>1-prop_treated_inp)) %>% 
      inner_join(covariate_tib, by="entry") )
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
    return( prob_treat %>%
              dplyr::mutate( u = runif(n = N_inp, min = 0, max = 1),
                      p_new = prop_score+u,
                      treated=as.numeric(dplyr::percent_rank(p_new)>
                                           1-prop_treated_inp),
                entry = seq_len(N_inp)) %>% 
              dplyr::select(entry, treated) %>%
              inner_join(covariate_tib, by="entry") )
  }
}




gen_covariates<-function(N_inp, overlap_inp, frac_shifted){
  #Generate several covariates, both observed and unobserved.
  #Shift means according to overlap_inp and frac_shifted 
  
  #first, create a set of correlated variables, mean zero, var/covar Sigma
  Sigma <- matrix(c(16 ,  4,  -4.8,
                    4  , 25,    9,
                    -4.8,  9,    9),3,3)
  corr_xs= MASS::mvrnorm(n = N_inp, rep(0, 3), Sigma)
  colnames(corr_xs)=c("obs_mvnorm1", "unobs_mvnorm", "obs_mvnorm2")
  
  return( 
    x_tib=tibble::tibble(
      entry = seq_len(N_inp),
      to_shift = as.numeric(dplyr::percent_rank(
        runif(n = N_inp,min = 0,max = 1)
                         )>1-frac_shifted)) %>%
    cbind(corr_xs) %>% 
    dplyr::group_by(to_shift) %>%
    dplyr::mutate(
      obs_beta=stats::rbeta(n = dplyr::n(),
                         shape1 = 6,
                         shape2 = 5+3*(to_shift*-(overlap_inp))),
      obs_binom=stats::rbinom(n = dplyr::n(),
                         size = 15,
                         prob = 0.5+0.25*(to_shift*(overlap_inp))),
      obs_norm=stats::rnorm(n = dplyr::n(),
                          mean = 0+5*(to_shift*(overlap_inp)),
                          sd = 10),
      unobs_beta=stats::rbeta(n = dplyr::n(),
                         shape1 = 6,
                         shape2 = 5+3*(to_shift*-(overlap_inp))),
      ubobs_binom=stats::rbinom(n = dplyr::n(),
                          size = 15,
                          prob = 0.5+0.25*(to_shift*(overlap_inp))),
      ubobs_norm=stats::rnorm(n = dplyr::n(),
                         mean = 0+5*(to_shift*(overlap_inp)),
                         sd = 10) ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-to_shift)
  )
  
}






