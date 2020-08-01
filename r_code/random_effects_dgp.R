library(dplyr)
library(furrr)
library(lubridate)
plan(multiprocess)
set.seed(1982)






helper_gen_xs<-function(N_input){
  #these are x variables drawn from a number of different distributions
  #there are a total of 24 variables, so that the Xs are Nx24
  Sigma <- matrix(c(16 ,  4,  -4.8,
                     4  , 25,    9,
                     -4.8,  9,    9),3,3)
  vs_all <- MASS::mvrnorm(n = N_input, rep(0, 3), Sigma)
  vs <- vs_all[,1:2]
  v3 <- stats::rbeta(n = N_input,
              shape1 = 2,
              shape2 = 2)
  v4 <- stats::rbinom(n = N_input, size = 10, prob = 0.2)
  v5 <- matrix(stats::rnorm(N_input * 20), N_input, 20)
  xs <- cbind(vs, v3, v4, v5)
  colnames(xs) <- glue::glue("v{1:24}")
  
  unobs_corr=vs_all[,3]
  unobs_rand1=stats::rbeta(n = N_input,
                   shape1 = 2,
                   shape2 = 2)
  unobs_rand2=stats::rbinom(n = N_input, size = 10, prob = 0.2)
  unobs_xs=cbind(unobs_corr,unobs_rand1,unobs_rand2 )
  colnames(unobs_xs) <- glue::glue("unobs_v{1:3}")
  
  xs <- tibble::as_tibble(cbind(xs, unobs_xs)) %>%
    dplyr::mutate(division_country = 1:N_input)

    return(xs)
}



helper_gen_vert_re<-function(){
  #Q for Ignacio: why this rather than creating a tibble outright?
  return(
    structure(
      list(
        #name of verticals
        vertical = c("A", "B", "C", "D", "E", "F", "G",
                     "H", "I", "J", "K", "L"),
        #probability of belonging in each vertical
        p = c(
          0.0629,
          0.039,
          0.0343,
          0.1943,
          0.0848,
          0.0895,
          0.12,
          0.0343,
          0.2352,
          0.001,
          0.0619,
          0.0428
        ),
        #random effect associated with each verical
        vertical_re = c(
          `vertical[1]` = 0.2212,
          `vertical[2]` = -0.3044,
          `vertical[3]` = 0.3183,
          `vertical[4]` = -0.6042,
          `vertical[5]` = -0.0863,
          `vertical[6]` = 0.5209,
          `vertical[7]` = -0.3082,
          `vertical[8]` = -0.0151,
          `vertical[9]` = 0.0912,
          `vertical[10]` = 0.0544,
          `vertical[11]` = 0.2699,
          `vertical[12]` = -0.1517
        )
      ),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA,-12L)
    )
  )
}


helper_assign_treat<-function(covariate_inp,unobs_inp, type_inp, N_inp, prop_treated_inp){
  #Type determines whether the selection into treatement is random or not
  if (type_inp == "random") {
    #if randome assignment of treatment
    #create a tibble with N rows (one per unit)
    #and sample treatment (0 if never, 1 if ever) using argument prop_treated
     return( tibble::tibble(
      division_country = 1:N_inp,
      treated = sample(
        x = c(0, 1),
        size = N_inp,
        prob = c(1 - prop_treated_inp, prop_treated_inp),
        replace = TRUE  )  )  )
  } 
  else{ #if treatment is based on covariates,
    #select a subset of the variables (v1,3,5,6,8,9,11,12,13,14,16,18,22,23) -- why these (ask ignacio)?
    #and store them in a matrix
    z <- covariate_inp %>%
      dplyr::select(-v2,
             -v4,-v7,
             -v10,
             -v15,
             -v17,-v20,
             -v21,
             -v24) %>%
      as.matrix()
    
    if(type_inp=="unobservables"){
      z <- unobs_inp %>% as.matrix()
    }
    
    
    #draw a vec of standard normal rvs, with length of the num of variables in z
    beta <- stats::rnorm(n = ncol(z),
                  mean = 0,
                  sd = 1)
    #compute the "z_hat" by multiplying the variables for each invidual by their beta
    #this will be N (num units) by 1 dimensional
    zbeta <- z %*% beta
    
    colnames(zbeta) <- 'zbeta'
    #zbeta can range from negative to positive, so we rescale it using the
    #logistic function to map it into the probability space (0,1)
    p <- zbeta %>% tibble::as_tibble() %>%
      dplyr::mutate(p = 1 / (1 + exp(zbeta)))
    
    #Ask ignacio: could we just take the top "prop_treated" from the logistic transformation above?
    #or, to have more randomness, subtract a random uniform from each and take the top frac?
    return( p %>% dplyr::ungroup() %>%
      dplyr::mutate(
        u = runif(n = N_inp, min = 0, max = 1),
        p_new = p+u,
        treated=as.numeric(percent_rank(p_new)>1-prop_treated_inp),
        division_country = 1:N_inp
        ) %>% dplyr::select(division_country, treated)
    
    )
    
    
    #while the average probability is different than our ultimate goal from prop_treated
    #adjust each individual p by the difference.
    #if the average p is below our goal, we will wind up adding to p. if above, we will subtract
    while (abs(mean(p$p) - prop_treated_inp) > 0.0001) {
      add_p <- prop_treated_inp - mean(p$p)
      p <- p %>% rowwise() %>%
        dplyr::mutate(p = min(1, p + add_p))
    }
    #to actually assign treatment, we take our adjusted p from above
    #we compare it to a random uniform (for each indiv)
    #if it's larger, we say they are treated, else not treated
    #store the division (just the row, or unit) and the treated status
    return( p %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        u = runif(n = N_inp, min = 0, max = 1),
        treated = case_when(p > u ~ 1,
                            p <= u ~ 0),
        division_country = 1:N_inp
      ) %>%
      dplyr::select(division_country, treated) )
  }
  
}

helper_first_period<-function(N_inp, N_time_inp,treat_start_inp, treated_tib){
  
  last_treat=ceiling(median(c(N_time_inp, treat_start_inp)))
  return(tibble::tibble(division_country = 1:N_inp,
         first = sample(
           x = treat_start_inp:(N_time_inp - 5),
           size = N_inp,
           replace = TRUE
         )) %>% 
    dplyr::inner_join(treated_tib, by = "division_country") %>%
    dplyr::mutate(first = case_when(treated == 1 ~ first,
                             TRUE ~ NA_integer_)) %>%
    dplyr::select(-treated))
}


helper_gen_countries<-function(){
  return(
    structure(
      list(
        country = 1:9,
        country_re = c(
          -0.9231,
          0.3533,
          0.6499,
          -0.0572,
          -1.1556,
          -0.8078,
          0.4037,
          0.0955,
          1.434
        ),
        p = c(
          0.0457,
          0.1019,
          0.1038,
          0.0067,
          0.0124,
          0.0152,
          0.0305,
          0.0048,
          0.679
        )
      ),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA,-9L)
    ))
}


helper_country_t_re_gen<-function(N_countries_inp, N_time_inp){
  time_country_grid <-
    tidyr::expand_grid(country = 1:N_countries_inp, time = 1:N_time_inp) 
  #tibble for each country, in time period 1, with a country specific normal shock
  r_one <-
    tibble::tibble(
      country  = 1:N_countries_inp,
      time = 1,
      r = stats::rnorm(n = N_countries_inp, mean = 0, sd = 1)
    )
  
  #merge in the period 1 shock (country specific) into the country x time tibble
  r <-
    full_join(time_country_grid, r_one, by = c("country", "time")) %>%
    dplyr::group_by(country)
  
  #now, fill out the rest of the country x time shocks
  #for each time period and country combination,
  #take the country specific shock from the time period before and add normal noise
  #For a given country, we thus have shock2=shock1+norm_shock
  for (t in 2:N_time_inp) {
    for (c in 1:N_countries_inp) {
      r$r[r$time == t &
            r$country == c] = r$r[r$time == (t - 1) &
                                    r$country == c] + stats::rnorm(n = 1, mean = 0, sd = 0.1)
    }
  }
  
  return(r)
}

#TODO(alexdkellogg): Treatment Effect and Decay adjustment by frequency
helper_treatment_effect_seq<-function(fake_data_inp, tau_one_zero_inp, tau_decay_inp){
  #find the longest post treat period of pitch 1 for any unit in data
  tau_one <- rep(NA, max(fake_data_inp$m + 1, na.rm = T))
  #define a treatment decay function: for period one, tau_one is an argument
  #then it decays by tau_decay*prior for each following period until the end of time
  for (m in seq_along(tau_one)) {
    if (m == 1) {
      tau_one[m] <- tau_one_zero_inp
    } else{
      tau_one[m] <- tau_one[m - 1] * tau_decay_inp
    }
  }
  #create a tibble that assings the treatment effect to each post treat period
  return(tau_one)
}


#TODO(alexdkellogg): Find an lapply/dplyr/more efficient way to do this
helper_num_pitches<-function(N_time_inp, N_inp, T_freq_inp, fake_data_inp){
  
  fake_data_inp=fake_data_inp %>%
    dplyr::mutate(new_pitches = 0,
           total_pitches = 0)
  #for each time and unit specific combination, we will create a variable for the number of
  #new pitches in that period (and a cumulative sum tracker)
  #no picthes happen until treatment, but a random amount happen every period after treatment
  for (i in 1:N_inp) {
    for (t in 1:N_time_inp) {
      #if the m is NA for a given entry, it will always be NA so we move to next entry
      if(is.na(fake_data_inp$m[fake_data_inp$division_country == i &
                        fake_data_inp$time == t]) ){
        break
      }
      # if we are in the pre treatment periods (m<0), set pitches to 0
      if (fake_data_inp$m[fake_data_inp$division_country == i &
                      fake_data_inp$time == t] < 0) {
        fake_data_inp$new_pitches[fake_data_inp$division_country == i &
                                fake_data_inp$time == t] <- 0
      }
      else if (fake_data_inp$m[fake_data_inp$division_country == i &
                           fake_data_inp$time == t] == 0) {
        # if we are in the treatment period
        # assign a random number of pitches (number of heads out of 10 coin flips) + minimum of 1
        fake_data_inp$new_pitches[fake_data_inp$division_country == i &
                                fake_data_inp$time == t] <- 1 +
          ifelse(T_freq_inp == "monthly", stats::rbinom(n = 1, size = 10, prob = 0.5),
                 ifelse(T_freq_inp == "weekly", stats::rbinom(n = 1, size = 3, prob = 0.4),
                        ifelse(T_freq_inp=="daily", purrr::rbernoulli(1, 2 / 9),
                               stats::rbinom(n = 1, size = 120, prob = 0.5)  )  ) )
        # if we are in the treatment period, total (cumulative pitches) is just the new pitched
        fake_data_inp$total_pitches[fake_data_inp$division_country == i &
                                  fake_data_inp$time == t] <-
          fake_data_inp$new_pitches[fake_data_inp$division_country == i &
                                  fake_data_inp$time == t]
      }
      else {
        # if in the post treat period, sample from the binomial without a min
        fake_data_inp$new_pitches[fake_data_inp$division_country == i &
                                fake_data_inp$time == t] <- ifelse(T_freq_inp == "monthly", stats::rbinom(n = 1, size = 10, prob = 0.5),
                                                                   ifelse(T_freq_inp == "weekly", stats::rbinom(n = 1, size = 3, prob = 0.4),
                                                                          ifelse(T_freq_inp=="daily", purrr::rbernoulli(1, 2 / 9),
                                                                                 stats::rbinom(n = 1, size = 120, prob = 0.5)  )  ) )
        
        # if we are in the proper post treat periods, we add the new pitches to the past cumulative sum
        fake_data_inp$total_pitches[fake_data_inp$division_country == i &
                                  fake_data_inp$time == t] <-
          fake_data_inp$total_pitches[fake_data_inp$division_country == i &
                                    fake_data_inp$time == t - 1] +
          fake_data_inp$new_pitches[fake_data_inp$division_country == i &
                                  fake_data_inp$time == t]
      }
    }
  }
  
  return(fake_data_inp %>%
    dplyr::mutate(total_pitches = case_when(treated == 1 ~ total_pitches,
                                     TRUE ~ 0)))
}



gen_time_constant_data <-function(N_inp, prop_treated_inp, type_inp) {
  
  #Generates tibble of long form panel data along with a number of x variables, the counterfactual outcome, and the target.
  #Each row is time period x unit unique combnation.
  
  #Args
  #N_inp: number of units to be generated
  #N_time_inp: number of total time periods
  #treat_start: first treatment period. For now, N_time_inp-treat_start>5
  #tau_one_zero: initial period treatment impact
  #tau_decay: fraction multiplied by tau_one_zero for each post treatment period after the impact (for all time)
  #thus, if set to 0.9, the impact is originally 1, then 0.9, then 0.81 in the subsequent periods
  #gamma: multiplier on the impact based on the total number of pitches (treatment intensity) -- impact+gamma*log(total_pitches)
  #seed: random number seed
  #prop_treated_inp: proportion of the N_inp entries that should receive treatment, roughly
  
  #rho_y: scales up the sd of the first period TE (normal sd_0=sd_e/(sqrt(1 - rho_y ^ 2)))
  #also determines how much the center of the normal dist for y_it moves
  #given prior errors --  mu_it <- rho_y * (y_t0[i, (t - 1)] - u[i, (t - 1)]) +u[i, t]
  #type_inp: one of "random" or "observables" (to fix -- unobservables?) determining assignment to treatment
  
  
  #Higher level description of data_full:
  #rows of the df represent period-entry combinations (eg N_inp (total num of entry) rows for period t).
  #each row should have a treatment indicator (treat_indicator), a period number (time_var), 
  #an individual ID (id_var), and an outcome (outcome_var)
  #for associated with that period-ID combination
  
  
  #Output
  #placebo_df_long, a dataframe of the same format as data_full, but now entirely
  #consists of donor units, some of which are placebo-treated (based on matching).

  # Generate time-invariant covariates
  
  #these are x variables drawn from a number of different distributions
  #there are a total of 24 variables, so that the Xs are Nx24
  xs_obs_unobs<-helper_gen_xs(N_inp)
  
  
  # Vertical
  #creating a tibble with 12 rows (one for each vertical), 
  #with columns for the name of the vertical, 
  #it's probability (ie prob that an observations falls into that vertical)
  #and the random effect shock associated with that variable
  vertical_re=helper_gen_vert_re()
  
  
  #crosswalk that creates a tibble with a row for each unit (N_inp)
  #and the vertical they have been randomly assigned to
  #assignment is done by sampling within creation, with replacement of course.
  #and according to the pmf defined in vertical_re$p
  verticalxwalk <- tibble::tibble(
    division_country = 1:N_inp,
    vertical = sample(
      x = LETTERS[1:12],
      size = N_inp,
      replace = TRUE,
      prob = vertical_re$p
    )
  )
  
  #each unit assigned a vertical, we map that back to their row in our data (xs) using join
  #Then, given the vertical they are assigned, we map on their vertical random effect 
  #(ASK Ignacio -- dropped immediately, not used for anything? -- sol: merged later)
  #then, we create a new variable, v25, that is 0 unless they are in vertical E or K
  xs_obs_unobs <- dplyr::inner_join(xs_obs_unobs, verticalxwalk, by = "division_country") %>%
    dplyr::inner_join(vertical_re, by = "vertical") %>%
    dplyr::mutate(v25 = case_when(vertical == "E" ~ -3.1,
                           vertical == "K" ~ -2.02,
                           TRUE ~ 0)) %>%
    dplyr::select(-vertical,-p,-vertical_re)
  
  
  
  # Treatment
  x_obs=xs_obs_unobs %>% dplyr::select(-c(division_country,contains("unobs")) )
  x_unobs=xs_obs_unobs %>% dplyr::select(contains("unobs"))
  treated<-helper_assign_treat(covariate_inp=x_obs,unobs_inp=x_unobs ,type_inp, N_inp, prop_treated_inp)
  
  
  # a_i
  #individual random effect for each unit, drawn from normal
  #output is tibble with ID and the indiv specific effect
  a <- tibble::tibble(division_country = 1:N_inp,
              a_i = stats::rnorm(n = N_inp, mean = 0.1748,
                          sd = 1.4678))
  
  
  # beta_i
  #define a random variable, specific to each unit
  beta <- tibble::tibble(division_country = 1:N_inp,
                 beta_i = stats::rnorm(n = N_inp, mean = 0,
                                sd = 0.07789))
  
  #create a tibble of unit, unit random effect (ai and betai),
  #unit vertical (with prob of assignment and the vertical Random Effect)
  #and unit ever treated status
  fake_data_unit_level <-dplyr::inner_join(a, verticalxwalk, by = "division_country") %>%
    dplyr::inner_join(vertical_re, by = "vertical") %>%
    dplyr::inner_join(treated, by = "division_country") %>%
    dplyr::inner_join(beta, by = c("division_country")) %>% 
    dplyr::inner_join(xs_obs_unobs, by= "division_country")
  
  
  return(fake_data_unit_level)
}




helper_date_shocks<-function(T_freq_inp,date_start_inp,
                             N_time_inp, time_constant_data){
  #generate the sequence of dates depending on the frequency
  period_dates=switch(T_freq_inp,
                      "daily"=format(ymd(date_start_inp)+days(0:(N_time_inp-1))),
                      "weekly"=format(ymd(date_start_inp)+weeks(0:(N_time_inp-1))),
                      "monthly"=format(ymd(date_start_inp)+months(0:(N_time_inp-1))),
                      "yearly"=format(ymd(date_start_inp)+years(0:(N_time_inp-1)))
  )
  
  #Identidy the relevant components of the date (day/week/etc)
  date_info_tib=switch (T_freq_inp,
                        "daily" = tibble::tibble(time=seq_len(N_time_inp),
                                         date_t=period_dates,
                                         day_num=day(period_dates), 
                                         week_num=week(period_dates),
                                         month_num=month(period_dates), 
                                         year_num=year(period_dates)),
                        "weekly"=tibble::tibble(time=seq_len(N_time_inp),
                                        date_t=period_dates,
                                        week_num=week(period_dates),
                                        month_num=month(period_dates), 
                                        year_num=year(period_dates)),
                        "monthly"=tibble::tibble(time=seq_len(N_time_inp),
                                         date_t=period_dates,
                                         month_num=month(period_dates), 
                                         year_num=year(period_dates)),
                        "yearly"=tibble::tibble(time=seq_len(N_time_inp),
                                        date_t=period_dates,
                                        year_num=year(period_dates))
  )
  
  
  date_time_grid=tidyr::expand_grid(time = seq_len(N_time_inp),
                                    division_country = time_constant_data$division_country) %>%
    dplyr::inner_join(date_info_tib, by="time")
  
  
  date_levels=factor(c("daily", "weekly", "monthly", "yearly"), ordered = T, levels =c("daily", "weekly", "monthly", "yearly") )
  agg_ids=c("day_num", "week_num", "month_num", "year_num")
  
  #Identify the larger dates
  freq_ind=match(T_freq_inp,date_levels)
  freq_vars=agg_ids[seq(freq_ind, length(agg_ids))]
  #define noise at each level, starting from daily, weekly, monthly, yearly
  freq_sd=c(0.1,0.15, 0.2,0.25)
  #add noise for the fine-ness of the data, as well as coarser dates (assumed independent)
  joint_sd=sum(freq_sd[seq(freq_ind,length(freq_sd))]**2)

  return(date_time_grid %>% dplyr::group_by_at(c(freq_vars)) %>%
           dplyr::mutate(date_shock=stats::rnorm(1,0, joint_sd)))
  
}

#daily

gen_time_varying_data<-function(time_constant_data,N_inp,  N_time_inp, T_freq_inp,
                                treat_start_inp, date_start_inp){
  
  #Generates tibble of long form panel data along with a number of x variables, the counterfactual outcome, and the target.
  #Each row is time period x unit unique combnation.
  
  #Args
  #N_inp: number of units to be generated
  #N_time_inp: number of total time periods
  #treat_start_inp: first treatment period. For now, N_time_inp-treat_start_inp>5
  #tau_one_zero_inp: initial period treatment impact
  #tau_decay: fraction multiplied by tau_one_zero_inp for each post treatment period after the impact (for all time)
  #thus, if set to 0.9, the impact is originally 1, then 0.9, then 0.81 in the subsequent periods
  #gamma: multiplier on the impact based on the total number of pitches (treatment intensity) -- impact+gamma*log(total_pitches)
  #seed: random number seed
  #prop_treated_inp: proportion of the N_inp entries that should receive treatment, roughly
  
  #rho_y: scales up the sd of the first period TE (normal sd_0=sd_e/(sqrt(1 - rho_y ^ 2)))
  #also determines how much the center of the normal dist for y_it moves
  #given prior errors --  mu_it <- rho_y * (y_t0[i, (t - 1)] - u[i, (t - 1)]) +u[i, t]
  #type: one of "random" or "observables" (to fix -- unobservables?) determining assignment to treatment
  
  
  #Higher level description of data_full:
  #rows of the df represent period-entry combinations (eg N_inp (total num of entry) rows for period t).
  #each row should have a treatment indicator (treat_indicator), a period number (time_var), 
  #an individual ID (id_var), and an outcome (outcome_var)
  #for associated with that period-ID combination
  
  
  #Output
  #placebo_df_long, a dataframe of the same format as data_full, but now entirely
  #consists of donor units, some of which are placebo-treated (based on matching).
  
  
  #create a date_time grid with shocks
  date_time_grid=helper_date_shocks(T_freq_inp=T_freq_inp,date_start_inp=date_start_inp,
                                    N_time_inp=N_time_inp, time_constant_data=time_constant_data)
  
  # d_it
  #creates a NxT by 2 tibble; each row is a unique entry (unit) x time combo
  time_grid <-tidyr::expand_grid(time = 1:N_time_inp,
                                 division_country = time_constant_data$division_country)
  
  
  
  ## countries
  #like the verticals, create a tibble of countries (9 total),
  #with their associated probability of assignment for each unit,
  #and the random effect associated with that country
  countries<-helper_gen_countries()
  
  N_countries <- length(countries$country)
  
  
  #map each row (Id) to a country by sampling from the list of countries
  #following the specified probability of each
  #output is a tibble with the Id and country assigned to said ID
  countryxwalk <- tibble::tibble(
    division_country = 1:N_inp,
    country = sample(
      x = countries$country,
      size = N_inp,
      replace = TRUE,
      prob = countries$p
    )
  )
  
  #recreate fake_data tibble
  #adding columns for country and country random effect (And removing assignment prob)
  #merge with time grid as well, to make it a long panel data (row is unit x time)
  #so far, all the vars are time constant and are thus repeated 
  fake_data_time_varying <-
    dplyr::inner_join(time_constant_data, countryxwalk, by = "division_country") %>%
    dplyr::select(-p) %>%
    dplyr::inner_join(countries, by = "country") %>%
    dplyr::select(-p) %>%
    dplyr::inner_join(time_grid, by = "division_country")
  
  
  #for each time period and country combination,
  #take the country specific shock from the time period before and add normal noise
  #For a given country, we thus have shock2=shock1+norm_shock
  r<-helper_country_t_re_gen(N_countries, N_time_inp)
  
  
  #after filling out the grid for each country/time pair, merge into fake data
  fake_data_time_varying <- dplyr::inner_join(fake_data_time_varying, r, by = c("country", "time"))
  
  ## pitches
  
  #First period is a tibble that assigns, for each unit
  #the first period of treatment. output is an Nx2 tibble.
  #Samples periods between 5 and Total_T-5 for each unit, and assigned 
  #9999L if they are in fact never treated
  #ask ignacio: can we parameterize the first treat (5)?
  #Also thinking of using N_inp-1, but the frac of treated for only 1 period would be too high
  #thoughts on sampling from a dist like Alex P?
  treated=time_constant_data %>% dplyr::select(division_country,treated)
  first_period <-helper_first_period(N_inp, N_time_inp,treat_start_inp,treated)
  
  #add to our fake data tibble two columns: the time of first treatment (unit specific)
  #and the relation of the current row to treat time (-1 would be period before treat)
  fake_data_time_varying <- dplyr::inner_join(fake_data_time_varying, first_period, by = "division_country") %>%
    dplyr::mutate(m = time - first)
  
  fake_data_time_varying=helper_num_pitches(N_time_inp, N_inp, T_freq_inp ,fake_data_time_varying)
  
  return(fake_data_time_varying %>% dplyr::inner_join(date_time_grid, by=c("division_country","time"))%>% arrange(time,division_country ) )
  
}




gen_synthetic_output<-function(synthetic_data_vars, tau_one_zero_inp, tau_decay_inp,
                               N_inp,  N_time_inp,
                               gamma_inp, rho_y_inp){
  
  # Initial y
  sd_e <-  0.5929
  sd_0 <- sd_e / (sqrt(1 - rho_y_inp ^ 2))
  
  #Generate period 1 y from a normal distribution
  #centered around the sum of individual random effects (ai and beta_i), 
  #vertical random effects (vertical_re), and time fixed effects (r)
  #y0 stores the unit and initial y value
  y0 <- synthetic_data_vars %>%
    filter(time == 1) %>%
    dplyr::mutate(
      mu_i_zero = a_i + vertical_re + r + beta_i + date_shock,
      y_0 = stats::rnorm(n = N_inp, mean = mu_i_zero, sd = sd_0)
    ) %>%
    dplyr::select(division_country, y_0)
  
  #create an NxTime matrix with the values indicating time relative to first treat (m)
  m_matrix <- synthetic_data_vars %>%
    dplyr::select(division_country, time, m) %>%
    tidyr::pivot_wider(names_from = time, values_from = m) %>%
    dplyr::select(-division_country) %>%
    as.matrix()
  
  # y matrix without treatment
  
  #create a NxTime matrix with values being the ith unit's mean in time j
  #the mean is the sum of random effects, which is where the normal y will be centered
  u <- synthetic_data_vars %>%
    dplyr::mutate(mu_i_t = a_i + vertical_re + r + beta_i * time + date_shock) %>%
    dplyr::select(division_country, time, mu_i_t) %>%
    tidyr::pivot_wider(names_from = time, values_from = mu_i_t) %>%
    dplyr::select(-division_country) %>%
    as.matrix()
  
  #initialize the non-treated outcome matrix, and assign the 
  #first col (time 0, unit 1:N_inp) from y0
  y_t0 <- matrix(nrow = N_inp, ncol = N_time_inp)
  y_t0[, 1] <- y0$y_0
  
  #generate the treated y from the same starting point
  y_t1 <- y_t0
  
  #generate a NxTime matrix, where value in i j is 
  #unit i's total pitches in period j
  total_pitches_matrix <- synthetic_data_vars %>%
    dplyr::select(division_country, time, total_pitches) %>%
    tidyr::pivot_wider(names_from = time, values_from = total_pitches) %>%
    dplyr::select(-division_country) %>%
    as.matrix()
  
  tau_one <-helper_treatment_effect_seq(synthetic_data_vars, tau_one_zero_inp, tau_decay_inp)
  
  #TODO(alexdkellogg): Treatment Effect and Growth adjustment -- way too big
  #for each unit, by each time period
  for (i in 1:N_inp) {
    for (t in 2:N_time_inp) {
      #define the period t mean as the difference between
      #i;s last period outcome and their last period mean 
      #(ie y_t0 has mean u, but is drawn from normal centered around u)
      #then scale that difference by rho, and add it to the current mean,
      #and draw a new outcome centered around this AR mean
      mu_it <- rho_y_inp * (y_t0[i, (t - 1)] - u[i, (t - 1)]) +
        u[i, t]
      y_t0[i, t] <- stats::rnorm(n = 1, mean = mu_it, sd = sd_e)
      
      #if the unit receives treatment in this period -- ie they are in the post period
      #and received treatment
      if (m_matrix[[i, t]] >= 0 & total_pitches_matrix[[i, t]] > 0) {
        #stop(total_pitches_matrix[i,t]<=0)
        #update only the treated matrix -- so that y_t0 is our counterfactual
        #update by adding the decayed impact along with 
        #gamma_inp (param) * the log of total pitches to date
        #Ask ignacio -- should we care about the marginal increase in pitches?
        y_t1[i, t] = y_t0[i, t] +
          tau_one[m_matrix[i, t] + 1]*(1+
                                         gamma_inp * log(total_pitches_matrix[i, t])) 
        
      } else{
        #if untreated, use counterfactual outcome
        y_t1[i, t] = y_t0[i, t]
      }
    }
  }
  
  
  
  colnames(y_t1)  <- colnames(y_t0) <- 1:N_time_inp
  
  #transform the counterfactual and treated matrices of y into tibbles
  #and pivot long. This yields a tibble with each row having a unique combo of 
  #time, unit and the affiliated y0 or y1
  y_t0_long <- tibble::as_tibble(y_t0) %>%
    dplyr::mutate(division_country = 1:N_inp) %>%
    tidyr::pivot_longer(names_to = 'time', cols = -division_country) %>%
    dplyr::mutate(time = as.integer(time)) %>%
    dplyr::rename(y0 = value)
  
  y_t1_long <- tibble::as_tibble(y_t1) %>%
    dplyr::mutate(division_country = 1:N_inp) %>%
    tidyr::pivot_longer(names_to = 'time', cols = -division_country) %>%
    dplyr::mutate(time = as.integer(time)) %>%
    dplyr::rename(y1 = value)
  
  #merge the treated and counterfactual 
  return(dplyr::inner_join(y_t0_long, y_t1_long, by = c("division_country", "time")))
}





gen_data <-function(N = 500,  N_time = 100, treat_start=45,
                   tau_one_zero = 0.02, tau_decay=0.9, gamma = 0.0009,
           seed = 1982, prop_treated = 0.5,  rho_y = 0.4, 
           type = c("random", "observables")) {
  
  #Generates tibble of long form panel data along with a number of x variables, the counterfactual outcome, and the target.
  #Each row is time period x unit unique combnation.
  
  #Args
  #N: number of units to be generated
  #N_time: number of total time periods
  #treat_start: first treatment period. For now, N_time-treat_start>5
  #tau_one_zero: initial period treatment impact
  #tau_decay: fraction multiplied by tau_one_zero for each post treatment period after the impact (for all time)
  #thus, if set to 0.9, the impact is originally 1, then 0.9, then 0.81 in the subsequent periods
  #gamma: multiplier on the impact based on the total number of pitches (treatment intensity) -- impact+gamma*log(total_pitches)
  #seed: random number seed
  #prop_treated: proportion of the N entries that should receive treatment, roughly
  
  #rho_y: scales up the sd of the first period TE (normal sd_0=sd_e/(sqrt(1 - rho_y ^ 2)))
  #also determines how much the center of the normal dist for y_it moves
  #given prior errors --  mu_it <- rho_y * (y_t0[i, (t - 1)] - u[i, (t - 1)]) +u[i, t]
  #type: one of "random" or "observables" (to fix -- unobservables?) determining assignment to treatment
  
  
  #Higher level description of data_full:
  #rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  #each row should have a treatment indicator (treat_indicator), a period number (time_var), 
  #an individual ID (id_var), and an outcome (outcome_var)
  #for associated with that period-ID combination
  
  
  #Output
  #placebo_df_long, a dataframe of the same format as data_full, but now entirely
  #consists of donor units, some of which are placebo-treated (based on matching).
  
  set.seed(seed)
  if (missing(type)) {
    type <- "random"
  } else{
    type <- match.arg(type)
  }
  
  stopifnot(prop_treated > 0 & prop_treated <= 1)
  stopifnot(N_time-treat_start > 5)
  
  fake_data_time_cons=gen_time_constant_data(N_inp=N, prop_treated_inp=prop_treated,
                                             type_inp=type)
  
  
  fake_data=gen_time_varying_data(time_constant_data=fake_data_time_cons,N_inp=N,
                                  N_time_inp=N_time, 
                                  treat_start_inp=treat_start)
  
  
  
  y_long=gen_synthetic_output(synthetic_data_vars=fake_data, tau_one_zero_inp=tau_one_zero,
                              tau_decay_inp=tau_decay,
                              N_inp=N,  N_time_inp=N_time, 
                              gamma_inp=gamma, rho_y_inp=rho_y)
    

    
    
    #merge the treated/counterfactual outcomes into fake_data
    #exponentiate all outcomes, and multiply by 1000
    #define y as the "observed" outcome in potential outcome framework
    fake_data <-
      dplyr::inner_join(fake_data, y_long, by = c("division_country", "time")) %>%
      dplyr::mutate(y = case_when(treated == 1 ~ y1 ,
                           treated == 0 ~ y0))  %>%
      dplyr::mutate(y0 = exp(y0) * 1000,
             y1 = exp(y1) * 1000,
             y = exp(y) * 1000)
    
    # # True lift
    # #q1 is the first quarter after treatment -- post treat period 0,1,2
    # q1 <- fake_data %>%
    #   filter(m %in% c(0, 1, 2),
    #          treated == 1)
    # 
    # #aggregate Y over that period for all indiv, time (and the counterfactual version)
    # sum_y0 <- sum(q1$y0)
    # sum_y1 <- sum(q1$y1)
    # 
    # #aggregate percent lift in q1
    # true_lift <- (sum_y1 / sum_y0) - 1
    
    ## Output data
    synthetic_data_full=fake_data %>%
      dplyr::rename(period=time,post_treat_t=m, Treatment_Period=first) %>%
      dplyr::mutate(treatperiod_0=case_when(
        period>=Treatment_Period~1,
        period<Treatment_Period~0
      )) %>% 
      dplyr::rename(entry=division_country) %>% arrange(period, entry) %>%
      dplyr::select(period, entry, Treatment_Period, new_pitches, total_pitches, 
             y0, y1, y, post_treat_t, everything())
    
  return(synthetic_data_full)
    
    # my_fake_data <- synthetic_data_full %>% dplyr::mutate(
    #   true_lift = true_lift,
    #   N = N,
    #   tau_one_zero = tau_one_zero,
    #   gamma = gamma,
    #   seed = seed,
    #   prop_treated = mean(synthetic_data_full$treated),
    #   rho_y = rho_y
    # ) %>%
    #   tidyr::nest(data = -c(true_lift, tau_one_zero, gamma, N, seed, prop_treated, rho_y))
    
}




gen_data_lubr <-function(N = 100,  date_start="2019-09-04",
                         date_end="2020-07-04",
                         first_treat="2020-01-04", 
                         T_freq=c("daily","weekly","monthly", "yearly"),
                    tau_one_zero = 0.02, tau_decay=0.9, gamma = 0.0009,
                    seed = 1982, prop_treated = 0.5,  rho_y = 0.4, 
                    type = c("random", "observables", "unobservables")) {
  
  #can allow all the way to the end?
  #unobservables added
  
  set.seed(seed)
  if (missing(type)) {
    type <- "random"
  } else{
    type <- match.arg(type)
  }
  
  if (missing(T_freq)) {
    T_freq <- "monthly"
  } else{
    T_freq <- match.arg(T_freq)
  }
  
  stopifnot(prop_treated > 0 & prop_treated <= 1)
  
  #given the dates and frequency, identify total number of periods
  N_time=switch(T_freq,
                "daily"=ceiling(lubridate::interval(ymd(date_start), ymd(date_end))/days(1)+1),
                "weekly"=ceiling(lubridate::interval(ymd(date_start), ymd(date_end))/weeks(1)+1),
                "monthly"=ceiling(lubridate::interval(ymd(date_start), ymd(date_end))/months(1)+1),
                "yearly"=ceiling(lubridate::interval(ymd(date_start), ymd(date_end))/years(1)+1)
  )
  treat_start_int=switch(T_freq,
                         "daily"=ceiling(lubridate::interval(ymd(date_start), ymd(first_treat))/days(1)+1),
                         "weekly"=ceiling(lubridate::interval(ymd(date_start), ymd(first_treat))/weeks(1)+1),
                         "monthly"=ceiling(lubridate::interval(ymd(date_start), ymd(first_treat))/months(1)+1),
                         "yearly"=ceiling(lubridate::interval(ymd(date_start), ymd(first_treat))/years(1)+1)
  )
  
  #Stop if there are too few pre or post treat periods
  stopifnot(treat_start_int > 10, 
            N_time-treat_start_int>5)

  
  fake_data_time_cons=gen_time_constant_data(N_inp=N, prop_treated_inp=prop_treated,
                                             type_inp=type)

  
  fake_data=gen_time_varying_data(time_constant_data=fake_data_time_cons,N_inp=N,
                                  date_start_inp=date_start,
                                  N_time_inp=N_time, T_freq_inp = T_freq,
                                  treat_start_inp=treat_start_int)
  
  
  y_long=gen_synthetic_output(synthetic_data_vars=fake_data, tau_one_zero_inp=tau_one_zero,
                              tau_decay_inp=tau_decay,
                              N_inp=N,  N_time_inp=N_time, 
                              gamma_inp=gamma, rho_y_inp=rho_y)
  
  
  
  
  #merge the treated/counterfactual outcomes into fake_data
  #exponentiate all outcomes, and multiply by 1000
  #define y as the "observed" outcome in potential outcome framework
  fake_data <-
    dplyr::inner_join(fake_data, y_long, by = c("division_country", "time")) %>%
    dplyr::mutate(y = case_when(treated == 1 ~ y1 ,
                         treated == 0 ~ y0))  %>%
    dplyr::mutate(y0 = exp(y0) * 1000,
           y1 = exp(y1) * 1000,
           y = exp(y) * 1000)
  
  # # True lift
  # #q1 is the first quarter after treatment -- post treat period 0,1,2
  # q1 <- fake_data %>%
  #   filter(m %in% c(0, 1, 2),
  #          treated == 1)
  # 
  # #aggregate Y over that period for all indiv, time (and the counterfactual version)
  # sum_y0 <- sum(q1$y0)
  # sum_y1 <- sum(q1$y1)
  # 
  # #aggregate percent lift in q1
  # true_lift <- (sum_y1 / sum_y0) - 1
  
  ## Output data
  synthetic_data_full=fake_data %>% 
    dplyr::rename(period=time,post_treat_t=m, Treatment_Period=first) %>%
    dplyr::mutate(treatperiod_0=case_when(
      is.na(Treatment_Period)~0,
      period>=Treatment_Period~1,
      period<Treatment_Period~0
    )) %>% 
    dplyr::rename(entry=division_country) %>% arrange(period, entry) %>%
    dplyr::select(date_t, period, entry, Treatment_Period, new_pitches, total_pitches, 
           y0, y1, y, post_treat_t, everything())
  
  return(synthetic_data_full)
  
  # my_fake_data <- synthetic_data_full %>% dplyr::mutate(
  #   true_lift = true_lift,
  #   N = N,
  #   tau_one_zero = tau_one_zero,
  #   gamma = gamma,
  #   seed = seed,
  #   prop_treated = mean(synthetic_data_full$treated),
  #   rho_y = rho_y
  # ) %>%
  #   tidyr::nest(data = -c(true_lift, tau_one_zero, gamma, N, seed, prop_treated, rho_y))
  
}





#Aggregate per week, month, quarter grouping by year/month/week, year/month, year/Quarter (quarter())
aggregate_data_by_date<-function(disagg_data, to="monthly", from="daily"){
  
  date_levels=factor(c("daily", "weekly", "monthly", "yearly"), ordered = T, levels =c("daily", "weekly", "monthly", "yearly") )
  #if no aggregation is happening, break
  stopifnot(match.arg(to,date_levels)>match.arg(from,date_levels))
  
  agg_ids=c("day_num", "week_num", "month_num", "year_num")
  
  #Identify the larger dates
  to_ind=match(to,date_levels)
  to_vars=agg_ids[seq(to_ind, length(agg_ids))]
  
  new_date_name=paste0(to, "_date_end", sep="")
  agg_data=disagg_data %>% dplyr::group_by_at(c("entry",to_vars)) %>%
    summarise(new_pitches=sum(new_pitches), total_pitches=sum(total_pitches),
              y0=sum(y0), y1=sum(y1), y=sum(y),
              !!as.name(new_date_name):=max(date_t),
              treatperiod_0=max(treatperiod_0)
              ) %>% dplyr::ungroup() 
  
  agg_data=agg_data %>% 
    dplyr::mutate(period=agg_data %>% dplyr::group_by_at(to_vars) %>% dplyr::select(tidyselect::all_of(to_vars)) %>% dplyr::group_indices()) 
  
  agg_data_covariates=agg_data %>% dplyr::inner_join(
    disagg_data %>% dplyr::select(-c(period, date_t, Treatment_Period, 
                              treatperiod_0, post_treat_t, 
                              r, date_shock, setdiff(agg_ids,to_vars),
                              new_pitches, total_pitches,
                              y0,y1, y)) %>%
      distinct_at(c("entry",to_vars), .keep_all = T), 
    by=c("entry", to_vars)
  ) %>% dplyr::group_by(entry) %>%
    dplyr::mutate(Treatment_Period=case_when(
      treated==0~NA_real_,
      treated==1~(length(treatperiod_0)-sum(treatperiod_0)+1)
    ),
    post_treat_t=period-Treatment_Period) %>% dplyr::ungroup()
  
  return(agg_data_covariates)

  
}




