library(dplyr)
library(furrr)
plan(multiprocess)
set.seed(1982)

gen_data <-function(N = 800,  N_time = 60, treat_start=15,
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
    # Generate time-invariant covariates
    
    #these are x variables drawn from a number of different distributions
    #there are a total of 24 variables, so that the Xs are Nx24
    Sigma <- matrix(c(10, 3, 3, 2), 2, 2)
    vs <- MASS::mvrnorm(n = N, rep(0, 2), Sigma)
    v3 <- rbeta(n = N,
                shape1 = 2,
                shape2 = 2)
    v4 <- rbinom(n = N, size = 10, prob = 0.2)
    v5 <- matrix(rnorm(N * 20), N, 20)
    #colnames(v5) <- glue::glue("v{5:24}") #redundant
    xs <- cbind(vs, v3, v4, v5)
    colnames(xs) <- glue::glue("v{1:24}")
    xs <- as_tibble(xs) %>%
      mutate(division_country = 1:N)
    
    
    # Vertical
    #creating a tibble with 12 rows (one for each vertical), 
    #with columns for the name of the vertical, 
    #it's probability (ie prob that an observations falls into that vertical)
    #and the random effect shock associated with that variable
    
    #Q for Ignacio: why this rather than creating a tibble outright?
    vertical_re <-
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
    
    #crosswalk that creates a tibble with a row for each unit (N)
    #and the vertical they have been randomly assigned to
    #assignment is done by sampling within creation, with replacement of course.
    #and according to the pmf defined in vertical_re$p
    verticalxwalk <- tibble(
      division_country = 1:N,
      vertical = sample(
        x = LETTERS[1:12],
        size = N,
        replace = TRUE,
        prob = vertical_re$p
      )
    )
    
    #each unit assigned a vertical, we map that back to their row in our data (xs) using join
    #Then, given the vertical they are assigned, we map on their vertical random effect 
    #(ASK Ignacio -- dropped immediately, not used for anything? -- sol: merged later)
    #then, we create a new variable, v25, that is 0 unless they are in vertical E or K
    xs <- inner_join(xs, verticalxwalk, by = "division_country") %>%
      inner_join(vertical_re, by = "vertical") %>%
      mutate(v25 = case_when(vertical == "E" ~ -3.1,
                             vertical == "K" ~ -2.02,
                             TRUE ~ 0)) %>%
      select(-vertical,-p,-vertical_re)
    
    
    
    # Treatment
    #Type determines whether the selection into treatement is random or not
    if (type == "random") {
      #if randome assignment of treatment
      #create a tibble with N rows (one per unit)
      #and sample treatment (0 if never, 1 if ever) using argument prop_treated
      treated <- tibble(
        division_country = 1:N,
        treated = sample(
          x = c(0, 1),
          size = N,
          prob = c(1 - prop_treated, prop_treated),
          replace = TRUE
        )
      )
    } else{ #if treatment is based on observables,
      #select a subset of the variables (v1,3,5,6,8,9,11,12,13,14,16,18,22,23) -- why these (ask ignacio)?
      #and store them in a matrix
      z <- xs %>%
        select(-division_country,
               -v2,
               -v4,-v7,
               -v10,
               -v15,
               -v17,-v20,
               -v21,
               -v24) %>%
        as.matrix()
      #draw a vec of standard normal rvs, with length of the num of variables in z
      beta <- rnorm(n = ncol(z),
                    mean = 0,
                    sd = 1)
      #compute the "z_hat" by multiplying the variables for each invidual by their beta
      #this will be N (num units) by 1 dimensional
      zbeta <- z %*% beta
      
      colnames(zbeta) <- 'zbeta'
      
      #zbeta can range from negative to positive, so we rescale it using the
      #logistic function to map it into the probability space (0,1)
      p <- zbeta %>% as_tibble() %>%
        mutate(p = 1 / (1 + exp(zbeta)))
      
      #Ask ignacio: could we just take the top "prop_treated" from the logistic transformation above?
      #or, to have more randomness, subtract a random uniform from each and take the top frac?
      
      #while the average probability is different than our ultimate goal from prop_treated
      #adjust each individual p by the difference.
      #if the average p is below our goal, we will wind up adding to p. if above, we will subtract
      while (abs(mean(p$p) - prop_treated) > 0.0001) {
        add_p <- prop_treated - mean(p$p)
        p <- p %>% rowwise() %>%
          mutate(p = min(1, p + add_p))
      }
      #to actually assign treatment, we take our adjusted p from above
      #we compare it to a random uniform (for each indiv)
      #if it's larger, we say they are treated, else not treated
      #store the division (just the row, or unit) and the treated status
      treated <- p %>%
        ungroup() %>%
        mutate(
          u = runif(n = N, min = 0, max = 1),
          treated = case_when(p > u ~ 1,
                              p <= u ~ 0),
          division_country = 1:N
        ) %>%
        select(division_country, treated)
    }
    
    #count up the num control and treated
    N_treated <- sum(treated$treated)
    N_control <- N - N_treated
    
    #First period is a tibble that assigns, for each unit
    #the first period of treatment. output is an Nx2 tibble.
    #Samples periods between 5 and Total_T-5 for each unit, and assigned 
    #9999L if they are in fact never treated
    #ask ignacio: can we parameterize the first treat (5)?
    #Also thinking of using N-1, but the frac of treated for only 1 period would be too high
    #thoughts on sampling from a dist like Alex P?
    first_period <-
      tibble(division_country = 1:N,
             first = sample(
               x = treat_start:(N_time - 5),
               size = N,
               replace = TRUE
             )) %>% 
      inner_join(treated, by = "division_country") %>%
      mutate(first = case_when(treated == 1 ~ first,
                               TRUE ~ 9999L)) %>%
      select(-treated)
    
    
    
    
    # a_i
    #individual random effect for each unit, drawn from normal
    #output is tibble with ID and the indiv specific effect
    a <- tibble(division_country = 1:N,
                a_i = rnorm(n = N, mean = 0.1748,
                            sd = 1.4678))
    
    #create a tibble of unit, unit random effect,
    #unit vertical (with prob of assignment and the vertical Random Effect)
    #and unit ever treated status
    fake_data <-
      inner_join(a, verticalxwalk, by = "division_country") %>%
      inner_join(vertical_re, by = "vertical") %>%
      inner_join(treated, by = "division_country")
    
    # d_it
    #creates a NxT by 2 tibble; each row is a unique entry (unit) x time combo
    time_grid <-
      tidyr::expand_grid(time = 1:N_time,
                         division_country = fake_data$division_country)
    
    mu_d = 0.01
    
    ## countries
    #like the verticals, create a tibble of countries (9 total),
    #with their associated probability of assignment for each unit,
    #and the random effect associated with that country
    countries <- structure(
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
    )
    
    N_countries <- length(countries$country)
    
    
    #map each row (Id) to a country by sampling from the list of countries
    #following the specified probability of each
    #output is a tibble with the Id and country assigned to said ID
    countryxwalk <- tibble(
      division_country = 1:N,
      country = sample(
        x = countries$country,
        size = N,
        replace = TRUE,
        prob = countries$p
      )
    )
    
    #recreate fake_data tibble
    #adding columns for country and country random effect (And removing assignment prob)
    #merge with time grid as well, to make it a long panel data (row is unit x time)
    #so far, all the vars are time constant and are thus repeated 
    fake_data <-
      inner_join(fake_data, countryxwalk, by = "division_country") %>%
      select(-p) %>%
      inner_join(countries, by = "country") %>%
      select(-p) %>%
      inner_join(time_grid, by = "division_country")
    
    
    ## Time
    #create a time random effect table -- each time period
    #is assigned a random shock from the normal dist below
    time_re <-
      tibble(time = 1:N_time,
             time_re = rnorm(n = N_time, mean = 0, sd = 0.141))
    #merge the random effect over time.
    #each individual in a given time period will face the same shock.
    fake_data <- inner_join(fake_data, time_re, by = "time")
    
    
    
    # time country interaction
    #redundant
    N_countries <- length(countries$country)
    
    #create a tibble that, in each row, has a unique combo of time and country
    time_country_grid <-
      tidyr::expand_grid(country = 1:N_countries, time = 1:N_time)
    #tibble for each country, in time period 1, with a country specific normal shock
    r_one <-
      tibble(
        country  = 1:N_countries,
        time = 1,
        r = rnorm(n = N_countries, mean = 0, sd = 1)
      )
    
    #merge in the period 1 shock (country specific) into the country x time tibble
    r <-
      full_join(time_country_grid, r_one, by = c("country", "time")) %>%
      group_by(country)
    
    #now, fill out the rest of the country x time shocks
    #for each time period and country combination,
    #take the country specific shock from the time period before and add normal noise
    #For a given country, we thus have shock2=shock1+norm_shock
    for (t in 2:N_time) {
      for (c in 1:N_countries) {
        r$r[r$time == t &
              r$country == c] = r$r[r$time == (t - 1) &
                                      r$country == c] + rnorm(n = 1, mean = 0, sd = 0.1)
      }
    }
    #after filling out the grid for each country/time pair, merge into fake data
    fake_data <- inner_join(fake_data, r, by = c("country", "time"))
    
    ## pitches
    #add to our fake data tibble two columns: the time of first treatment (unit specific)
    #and the relation of the current row to treat time (-1 would be period before treat)
    fake_data <-
      inner_join(fake_data, first_period, by = "division_country") %>%
      mutate(m = time - first)
    #find the longest post treat period of pitch 1 for any unit in data
    tau_one <- rep(NA, max(fake_data$m + 1))
    #define a treatment decay function: for period one, tau_one is an argument
    #then it decays by tau_decay*prior for each following period until the end of time
    for (m in 1:max(fake_data$m + 1)) {
      if (m == 1) {
        tau_one[m] <- tau_one_zero
      } else{
        tau_one[m] <- tau_one[m - 1] * tau_decay
      }
    }
    #create a tibble that assings the treatment effect to each post treat period
    tau <- tibble(m = 0:max(fake_data$m),
                  tau_one = tau_one)
    
    fake_data <- fake_data %>%
      mutate(new_pitches = 0,
             total_pitches = 0)
    
    #for each time and unit specific combination, we will create a variable for the number of
    #new pitches in that period (and a cumulative sum tracker)
    #no picthes happen until treatment, but a random amount happen every period after treatment
    for (t in 1:N_time) {
      for (i in 1:N) {
        #if we are in the pre treatment periods (m<0), set pitches to 0
        if (fake_data$m[fake_data$division_country == i &
                        fake_data$time == t] < 0) {
          fake_data$new_pitches[fake_data$division_country == i &
                                  fake_data$time == t] <- 0
        } else{
          #if we are in the treatment period (or post treat periods)
          #assign a random number of pitches (number of heads out of 10 coin flips)
          #this is for each post treat period
          fake_data$new_pitches[fake_data$division_country == i &
                                  fake_data$time == t] <- rbinom(n = 1,
                                                                 size = 10,
                                                                 prob = 0.5)
        }
        #if we are in the treatment period, total (cumulative pitches) is just the new pitched
        if (fake_data$m[fake_data$division_country == i &
                        fake_data$time == t] == 0) {
          fake_data$total_pitches[fake_data$division_country == i &
                                    fake_data$time == t] <-
            fake_data$new_pitches[fake_data$division_country == i &
                                    fake_data$time == t]
        }
        #if we are in the proper post treat periods, we add the new pitches to the past cumulative sum
        else if (fake_data$m[fake_data$division_country == i &
                               fake_data$time == t] > 0) {
          fake_data$total_pitches[fake_data$division_country == i &
                                    fake_data$time == t] <-
            fake_data$total_pitches[fake_data$division_country == i &
                                      fake_data$time == t - 1] +
            fake_data$new_pitches[fake_data$division_country == i &
                                    fake_data$time == t]
        }
        
        
      }
    }
    
    fake_data <- fake_data %>%
      mutate(total_pitches = case_when(treated == 1 ~ total_pitches,
                                       TRUE ~ 0))
    
    
    
    # beta_i
    #define a random variable, specific to each unit
    beta <- tibble(division_country = 1:N,
                   beta_i = rnorm(n = N, mean = 0,
                                  sd = 0.07789))
    #add this random effect to the data
    fake_data <-
      inner_join(fake_data, beta, by = c("division_country"))
    
    
    
    # Initial y
    sd_e <-  0.5929
    sd_0 <- sd_e / (sqrt(1 - rho_y ^ 2))
    
    #Generate period 1 y from a normal distribution
    #centered around the sum of individual random effects (ai and beta_i), 
    #vertical random effects (vertical_re), and time fixed effects (r)
    #y0 stores the unit and initial y value
    y0 <- fake_data %>%
      filter(time == 1) %>%
      mutate(
        mu_i_zero = a_i + vertical_re + r + beta_i,
        y_0 = rnorm(n = N, mean = mu_i_zero, sd = sd_0)
      ) %>%
      select(division_country, y_0)
    
    #create an NxTime matrix with the values indicating time relative to first treat (m)
    m_matrix <- fake_data %>%
      select(division_country, time, m) %>%
      tidyr::pivot_wider(names_from = time, values_from = m) %>%
      select(-division_country) %>%
      as.matrix()
    
    # y matrix without treatment
    
    #create a NxTime matrix with values being the ith unit's mean in time j
    #the mean is the sum of random effects, which is where the normal y will be centered
    u <- fake_data %>%
      mutate(mu_i_t = a_i + vertical_re + r + beta_i * time) %>%
      select(division_country, time, mu_i_t) %>%
      tidyr::pivot_wider(names_from = time, values_from = mu_i_t) %>%
      select(-division_country) %>%
      as.matrix()
    
    #initialize the non-treated outcome matrix, and assign the 
    #first col (time 0, unit 1:N) from y0
    y_t0 <- matrix(nrow = N, ncol = N_time)
    y_t0[, 1] <- y0$y_0
    
    #generate the treated y from the same starting point
    y_t1 <- y_t0
    
    #generate a NxTime matrix, where value in i j is 
    #unit i's total pitches in period j
    total_pitches_matrix <- fake_data %>%
      select(division_country, time, total_pitches) %>%
      tidyr::pivot_wider(names_from = time, values_from = total_pitches) %>%
      select(-division_country) %>%
      as.matrix()
    
    
    #for each unit, by each time period
    for (i in 1:N) {
      for (t in 2:N_time) {
        #define the period t mean as the difference between
        #i;s last period outcome and their last period mean 
        #(ie y_t0 has mean u, but is drawn from normal centered around u)
        #then scale that difference by rho, and add it to the current mean,
        #and draw a new outcome centered around this AR mean
        mu_it <- rho_y * (y_t0[i, (t - 1)] - u[i, (t - 1)]) +
          u[i, t]
        y_t0[i, t] <- rnorm(n = 1, mean = mu_it, sd = sd_e)
        
        #if the unit receives treatment in this period -- ie they are in the post period
        #and received treatment
        if (m_matrix[[i, t]] >= 0 & total_pitches_matrix[[i, t]] > 0) {
          #stop(total_pitches_matrix[i,t]<=0)
          #update only the treated matrix -- so that y_t0 is our counterfactual
          #update by adding the decayed impact along with 
          #gamma (param) * the log of total pitches to date
          #Ask ignacio -- should we care about the marginal increase in pitches?
          y_t1[i, t] = y_t0[i, t] +
            tau_one[m_matrix[i, t] + 1] +
            gamma * log(total_pitches_matrix[i, t])
        } else{
          #if untreated, use counterfactual outcome
          y_t1[i, t] = y_t0[i, t]
        }
      }
    }
    
    
    
    colnames(y_t1)  <- colnames(y_t0) <- 1:N_time
    
    #transform the counterfactual and treated matrices of y into tibbles
    #and pivot long. This yields a tibble with each row having a unique combo of 
    #time, unit and the affiliated y0 or y1
    y_t0_long <- as_tibble(y_t0) %>%
      mutate(division_country = 1:N) %>%
      tidyr::pivot_longer(names_to = 'time', cols = -division_country) %>%
      mutate(time = as.integer(time)) %>%
      rename(y0 = value)
    
    y_t1_long <- as_tibble(y_t1) %>%
      mutate(division_country = 1:N) %>%
      tidyr::pivot_longer(names_to = 'time', cols = -division_country) %>%
      mutate(time = as.integer(time)) %>%
      rename(y1 = value)
    
    #merge the treated and counterfactual 
    y_long <-
      inner_join(y_t0_long, y_t1_long, by = c("division_country", "time"))
    
    #merge the treated/counterfactual outcomes into fake_data
    #exponentiate all outcomes, and multiply by 1000
    #define y as the "observed" outcome in potential outcome framework
    fake_data <-
      inner_join(fake_data, y_long, by = c("division_country", "time")) %>%
      mutate(y = case_when(treated == 1 ~ y1 ,
                           treated == 0 ~ y0))  %>%
      mutate(y0 = exp(y0) * 1000,
             y1 = exp(y1) * 1000,
             y = exp(y) * 1000)
    
    # True lift
    #q1 is the first quarter after treatment -- post treat period 0,1,2
    q1 <- fake_data %>%
      filter(m %in% c(0, 1, 2),
             treated == 1)
    
    #aggregate Y over that period for all indiv, time (and the counterfactual version)
    sum_y0 <- sum(q1$y0)
    sum_y1 <- sum(q1$y1)
    
    #aggregate percent lift in q1
    true_lift <- (sum_y1 / sum_y0) - 1
    
    ## Output data
    #just the x variables, removing v25 (the randomly effected based on verticals)
    xs <- xs %>%
      select(-v25)
    #Tibble of the relevant output data 
    # synthetic_data <- fake_data %>%
    #   dplyr::select(division_country,
    #                 time,
    #                 treated,
    #                 country,
    #                 vertical,
    #                 total_pitches,
    #                 y)  %>%
    #   inner_join(xs, by = "division_country")
    
    synthetic_data_full=fake_data %>%
      dplyr::select(-y1)  %>% 
      rename(period=time,post_treat_t=m, Treatment_Period=first,
             counterfactual=y0, target=y) %>%
      inner_join(xs, by = "division_country") %>% 
      rename(entry=division_country) %>% arrange(period, entry)
    

    
    # my_fake_data <- synthetic_data_full %>% mutate(
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



reformat_random_effects_dgp<-function(re_fake_data, id_var="entry", time_var="period",
                                      treat_time_var="Treatment_Period", treat_indicator_var="treatperiod_0"){
  #Transforms data from the Random Effects DGP (which has a standard output tibble) to a format usaable with causal_panel_benchmark_functions
  
  #Args
  #re_fake_data: long-form dataframe output from "gen_data" function in random_effects_dgp.R.
  #id_var: column name for the new numeric, unique ID representing the entry (unit) 
  #time_var:column name for the new numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #treat_indicator:column name of the new binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)
  #treat_time_var: column name for the new numeric variable indicating, for each entry, it's first treatment exposure
  
  #output: renaming and reformatting the output from gen_data so that the columns work with future estimation functions in causal_panel_benchmark_functions
  
  reformatted_data=re_fake_data %>% rename(!!as.name(time_var):=time, !!as.name(id_var):=division_country) %>% 
    select(!!as.name(time_var), !!as.name(id_var), treated, total_pitches, y) %>% group_by(!!as.name(id_var)) %>%
    mutate(
      !!as.name(treat_time_var):=case_when(
        treated==1~(length(total_pitches)-sum(total_pitches>0)+1),
        treated==0~NA_real_),
      !!as.name(treat_indicator_var):=case_when(
        is.na(!!as.name(treat_time_var))~0,
        !!as.name(time_var)<!!as.name(treat_time_var)~0,
        !!as.name(time_var)>=!!as.name(treat_time_var)~1
      )
    ) 
  
  return(reformatted_data)
  
}