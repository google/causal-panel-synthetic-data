#Try to write functions for everything:

#ingest data can be a function

#Run the TE can be a function of 
#AA or AB (which type of synthetic analysis is the goal)
#How many types of data do you want to run? -- is there a way to do this in parallel??
#maybe want to create helper functions so that we estimate (and return) Gsynth separately from
#CausalImpact, and SCDID (once I have the proper one)
#That way, I can call the furr r package and map each dataset to the gsynth_estimation etc.

#can then create print plot methods
#to print the Density of ATT, the Gap plots, as well as printing tables


library(CausalImpact)
library(Matrix)

library(gsynth)
library(augsynth)
library(tidyverse)
library(panelView)
library(synthdid)
library(resample)
library(mvtnorm)
library(janitor)
library("qpcR")
library(dplyr)
library(furrr)
library(gridExtra)
library(quadprog)



############################################
#Helper functions for AA testing: selecting the placebo treated

#Potential Improvements/Options: 
#First -- find K possible matches for each treated, without replacement. 
#If there are ties, break them based on how good the second best is, 

#Second -- threshold cutoff? If the closest placebo-treated is far,
#just throw out that treated.
###########################################


nearest_ts_euclidean=function(ts_tomatch, ts_rest){
  #Helper function to matching_without_replacement
  #Measures the euclidean distance between ts_tomatch
  #and each of the rows in the dataframe ts_rest.
  
  #Args
  #ts_tomatch: treated time series dataframe, rows identify ID and columns Time
  #ts_rest: donor time series dataframe, rows identify ID and columns Time.
  
  #Output
  #Vector of indices indicating the row in ts_rest that are best match 
  #(min L2 norm distance) to the row in ts_tomatch. 
  #(Vector length equal to number of rows in ts_tomatch, ie treated entries)
  #If ts_tomatch has more than 1 entry, the matching is done with replacement
  #(multiple treated can have the same match).
  
  apply(ts_tomatch,1,function(ts_tomatch) {
    which.min(
      apply(ts_rest,1,function(ts_rest,ts_tomatch) {
      dist(rbind(ts_rest,ts_tomatch))
    },
    ts_tomatch)
    )
  }
  )
}


matching_without_replacement=function(treated_block, control_block, id_var="entry"){
  #finds the nearest match for each unit in the treated subsample (treated_block)
  #among the donor pool (control_block) without replacement via individual calls
  #to nearest_ts_euclidean.
  
  #Args
  #treated_block: treated time series dataframe, rows identify ID and columns Time
  #control_block: donor time series dataframe, rows identify ID and columns Time.
  
  #Output
  #df_toreturn, a dataframe containing a column for the placebo-treated unit ID numbers,
  #the treated unit it was the nearest match to, and the time that treated unit was actually
  #treated (num_rows of the dataframe equal to num_rows of treated_block)
  
  #store an empty vector for the donor IDs that match
  already_matched=c()
  #Store the time of treatment and treated ID for the true treated units
  placebo_treat_period=treated_block %>% pull(Treatment_Period)
  treatment_unit=treated_block %>% pull(!!as.name(id_var))
  for(i in 1:nrow(treated_block)){
    if(i==1){
      #if we are searching for the match of our first treated unit, we can search across all donors
      temp_match=nearest_ts_euclidean(treated_block %>% slice(i) %>% dplyr::select(-c(!!as.name(id_var), Treatment_Period)), control_block %>% dplyr::select(-!!as.name(id_var)))
      already_matched[i]=control_block %>% slice(temp_match) %>% pull(!!as.name(id_var))
    }
    
    if(i!=1){
      #If we have already found a match, restrict the search for future matches to the subset of currently unmatched donors
      temp_match=nearest_ts_euclidean(treated_block %>% slice(i) %>% 
                                        dplyr::select(-c(!!as.name(id_var), Treatment_Period)), control_block %>% 
                                        filter(!!as.name(id_var) %in% setdiff(!!as.name(id_var), already_matched)) %>% dplyr::select(-!!as.name(id_var)))
      
      already_matched[i]=control_block %>% 
        filter(!!as.name(id_var) %in% setdiff(!!as.name(id_var), already_matched)) %>% 
        slice(temp_match) %>% pull(!!as.name(id_var))
    }
  }
  #Store the resulting vectors in a dataframe for output
  df_toreturn=data.frame(temp_id=already_matched, Treatment_Period=placebo_treat_period, Treatment_Unit=treatment_unit)
  df_toreturn=df_toreturn%>% rename(!!as.name(id_var):=temp_id )
  return(df_toreturn)
}




#To-Do: expand for the case when we have covariates? This could be useful in 1) matching the placebo-treated,
#2) improving performance and range of tests available to use

#To-Do: if covariates are allowed, must adapt functions (pivot wide/long, specifically)
create_placebo_df <- function(data_full, id_var="entry", time_var="period", treat_indicator="treatperiod_0",
                              outcome_var="target", counterfac_var="counter_factual"){
  #Generates a placebo only dataframe, using matching methods to select 
  #placebo-treated entries as those most similar to truly-treated.
  
  #Args
  #data_full: long-form dataframe with both treated and control entries.
  #id_var: column name of numeric, unique ID representing the entry (unit) 
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #treat_indicator:column name of binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)
 
   #Higher level description of data_full:
  #rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  #each row should have a treatment indicator (treat_indicator), a period number (time_var), 
  #an individual ID (id_var), and an outcome (outcome_var)
  #for associated with that period-ID combination
  
  
  #Output
  #placebo_df_long, a dataframe of the same format as data_full, but now entirely
  #consists of donor units, some of which are placebo-treated (based on matching).
  

  
  #Split the dataset based on whether they are ever treated
  tr.entries  <- data_full %>% filter(!!as.name(treat_indicator)>0) %>% distinct( !!as.name(id_var) ) %>% pull()
  ct.entries  <- setdiff(data_full %>% distinct( !!as.name(id_var) ) %>% pull(), tr.entries)
  #Create a dataframe of the subset of control units
  cd  <- data_full %>% filter(!!as.name(id_var) %in% ct.entries)
  cd  <- merge(cd, data.frame(entry = ct.entries, rank = seq_along(ct.entries)))
  
  #Pivot the long cd dataframe to wide..now, each row represents an id_var, with columns for the outcome_var at each time period 
  cd_for_match=pivot_wider(data=cd %>% arrange(!!as.name(time_var), !!as.name(id_var)), names_from = !!as.name(time_var), id_cols = c(!!as.name(id_var)), values_from=c(!!as.name(outcome_var)))
  
  #After identifyin the treated observations (for now, one by one), want to find the placebo match
  treated_to_match  <- data_full %>% filter(!!as.name(id_var) %in% tr.entries)

  #We will want to make the data wide???
  #want a treatment data indicator before doing this as well
  treated_to_match=treated_to_match %>% group_by(!!as.name(id_var)) %>% mutate(Treatment_Period=length(!!as.name(treat_indicator))-sum(!!as.name(treat_indicator))+1) %>% ungroup()

  #for pivoting, potential issues arise if we have several time varying covariates
  #we'd have to take the values_from each of them, and for any constant args we'd presumably have to add them to id_cols
  data_wide_m <- pivot_wider(data=treated_to_match, names_from = !!as.name(time_var), id_cols = c(!!as.name(id_var), Treatment_Period), values_from=c(!!as.name(outcome_var)))  #%>% as.matrix() 
  
  #Supply the treated units for matching (in wide form) and the donor units (in wide form) to the match function
  matched_placebo_df_temp=matching_without_replacement(data_wide_m, cd_for_match, id_var)
  
  #We now have the set of control units that form the placebo set, along with their placebo Treat Period and their correspondning Treat Unit
  #Merge this DF into the cd_for_match using id_var as the key
  
  placebo_df_wide=cd_for_match %>% left_join(matched_placebo_df_temp, by=id_var)

  #Pivot_long it and recreate the treat_indicator indicator based on Treatment_Period, so that the ouput matches the input
  #If there are covariates or other variables which we do not want to be pivoting, append them to non_time_vars vectors
  #may have to add values_to as well if we have multiple time varying variables
  non_time_vars=c(id_var, "Treatment_Period", "Treatment_Unit")
  placebo_df_long=placebo_df_wide %>% pivot_longer(-all_of(non_time_vars), names_to = time_var, values_to=outcome_var) %>%
    mutate(!!as.name(time_var):=as.numeric(!!as.name(time_var))) %>%
    arrange(!!as.name(time_var), !!as.name(id_var)) %>% mutate(!!as.name(treat_indicator):= case_when(
      is.na(Treatment_Period)~0,
      !!as.name(time_var)<Treatment_Period~0,
      !!as.name(time_var)>=Treatment_Period~1)) %>% left_join(
        data_full %>% 
          dplyr::select(id_var, time_var, counterfac_var), by=c(id_var, time_var)
      )
  
  return(placebo_df_long)
}




############################################
#Functions for estimating Treatment Effects

#Potential Improvements/Options: 
#Could add a units_to_estimate parameter to control directly how mnay units to run through
###########################################


estimate_causalimpact_series <- function(data_full, id_var="entry", time_var="period", 
                                         treat_indicator="treatperiod_0", outcome_var="target",
                                         counterfac_var="counter_factual"){
  #Estimates CausalImpact treatment effects given a long form data set, outputting a dataframe
  #consisting of a series of treatments effects for each id_var by time_var in all  periods
  
  #Args
  #data_full: long-form dataframe with both treated and control entries, sorted by period and entry (ascending).
  #id_var: column name of numeric, unique ID representing the entry (unit) 
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #treat_indicator:column name of binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)
  
  #Higher level description of data_full:
  #rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  #each row should have a treatment indicator (treat_indicator), a period number (time_var), 
  #an individual ID (id_var), and an outcome (outcome_var)
  #for associated with that period-ID combination
  
  
  #Output
  #Dataframe containing the full series of the outcome (all T), as well as predicted (counterfactual) outcome
  #and the associated effects by id_Var, time_var
 
  
  #Split the dataset based on whether they are ever treated
  tr.entries  <- data_full %>% filter(!!as.name(treat_indicator)>0) %>% distinct( !!as.name(id_var) ) %>% pull() %>% sort()
  ct.entries  <- setdiff(data_full %>% distinct( !!as.name(id_var) ) %>% pull(), tr.entries)
  #Create a dataframe of the subset of control units
  cd  <- data_full %>% filter(!!as.name(id_var) %in% ct.entries)
  cd  <- merge(cd, data.frame(entry = ct.entries, rank = seq_along(ct.entries)))
  
  
  
  # construct control data matrix
  #For all the control data: Row of matrix (index i) indicates the time, Column (j) indicates observation 
  #matrix entry represents the Target
  donor_outcome_matrix <- as.matrix(sparseMatrix(x = cd[[outcome_var]], i  = cd[[time_var]], j = cd$rank))
  
  
  
  #To output Parallelized process, we gather the inputs as a series of lists
  #first, we want the full treated dataframe to manipulate
  td_all  <- data_full %>% filter(!!as.name(id_var) %in% tr.entries) %>% arrange(!!as.name(time_var), !!as.name(id_var)) #NEW, MIGHT CAUSE ERROR
  
  #Parallelized computation of the data matrix, by id_var, required for causal impact
  #specifically,  append treated data as the first column in the 
  #donor_outcome_matrix matrix (described above, col=obs, row=time, entry=outcome)
  #matrix is now TREATED_DATA in col 1 for all of time, and then all control data 
  list_of_input_data=td_all%>%
    split(.[[id_var]]) %>% future_map(~ .[[outcome_var]]) %>% future_map(~cbind(.,donor_outcome_matrix, deparse.level = 0))
  
  #Parallelized computation of the pre-treatment range, by id_var
  list_of_pretreat_ranges=td_all%>%
    split(.[[id_var]]) %>% future_map(~ .[[treat_indicator]]) %>% future_map(~range(1,sum(.==0)))
  #Parallelized computation of the post-treatment range, by id_var
  list_of_posttreat_ranges=td_all%>%
    split(.[[id_var]]) %>% future_map(~ .[[treat_indicator]]) %>% future_map(~range(sum(.==0)+1,sum(.==0)+sum(.==1) ))
  
  #Parallelized computation of the causal impact
  list_of_causalimpact_series=future_pmap(list(list_of_input_data, list_of_pretreat_ranges, list_of_posttreat_ranges), CausalImpact) %>% future_map(~.$series) %>% 
    future_map(~as.data.frame(.)) 
  
  period_entry_rowlabs_df=data.frame(tempid=rep(sort(tr.entries), each=max(data_full[[time_var]])), temp_t=rep(1:max(data_full[[time_var]]), times=length(tr.entries)))
  names(period_entry_rowlabs_df)=c(id_var, time_var)
  
  causalimpact_series_output= list_of_causalimpact_series%>% do.call(bind_rows, .) %>% dplyr::select(c(response,point.pred, point.effect, point.effect.lower, point.effect.upper)) %>%
    cbind(., period_entry_rowlabs_df, deparse.level = 0 ) %>% 
    inner_join(td_all %>% 
                 group_by(!!as.name(id_var)) %>% mutate(Treatment_Period=length(!!as.name(treat_indicator))-sum(!!as.name(treat_indicator))+1) %>% 
                 ungroup() %>%
                 dplyr::select(id_var, time_var, Treatment_Period), by=c(id_var, time_var)
    ) %>% arrange(!!as.name(time_var), !!as.name(id_var))
  
  
  if(!is.null(counterfac_var)){
    causalimpact_series_output=causalimpact_series_output %>%inner_join(
      td_all %>% 
        dplyr::select(id_var, time_var, counterfac_var), by=c(id_var, time_var)
    )%>% mutate(
      cf_point.effect=(response-!!as.name(counterfac_var)),
      cf_pct.effect=(response/!!as.name(counterfac_var))-1
    )
  }
  
  #add a column with relative (pct) effect
  causalimpact_series_output=causalimpact_series_output %>% mutate(
    pct.effect=(response/point.pred)-1
  )
  
  rownames(causalimpact_series_output)=NULL
  
  
  
  
  return(causalimpact_series_output)

}


# Function to compute aggregate effects by period-- perhaps give the options for MSE, MAE
# compute_causalimpact_aggregates<-function(causalimpact_series_df, time_var){
#   
#   aggregated_causal_impact_placebo_series[[synth_type]]=causalimpact_series_forbootstrap[[synth_type]] %>% group_by(period) %>%
#     summarise(mean_abs_effect=mean(point.effect),
#               mean_abs_lower=mean(point.effect.lower),
#               mean_abs_upper=mean(point.effect.upper),
#               med_abs_effect=median(point.effect),
#               med_abs_lower=median(point.effect.lower),
#               med_abs_upper=median(point.effect.upper),
#               mean_percent_effect=mean(pct.effect),
#               med_percent_effect=median(pct.effect)  ) %>% ungroup() 
#   
#   names(aggregated_causal_impact_placebo_series)[synth_type]=names(data_list)[synth_type]
# }



#Potential Improvement: can I use ... as an argument so that any other gsynth arguments (ex, k=5) can be passed? right now, no way to do that for user
estimate_gsynth_series <- function(data_full, id_var="entry", time_var="period", treat_indicator="treatperiod_0", outcome_var="target", X_in=NULL,
                                   counterfac_var="counter_factual",se_est=TRUE, num_boots=1000, inference_type="parametric",
                                   factor_range=c(0,5),force_FE="unit", cross_val=TRUE, 
                                   EM_flag=FALSE, estimator_type="ife", 
                                   parallel_boot=FALSE){
  #Estimates Gsynth treatment effects given a long form data set, outputting a dataframe
  #consisting of a series of treatments effects for each id_var by time_var in all periods
  
  #Args
  #data_full: long-form dataframe with both treated and control entries, sorted by period and entry (ascending).
  #id_var: column name of numeric, unique ID representing the entry (unit) 
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #outcome_var: the y var for the time series
  #treat_indicator:column name of binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)
  #X_in: time-varying covariates
  #note: a panel reg formula is implicitly defined as outcome_var~treat_indicator+X_in
  #for the remaining parameters, these are options into Gsynth (use ?gsynth for info)
  #se_est: boolean as to whether uncertainty estimates are provided/estimated
  #num_boots: number of bootstraps run to obtain SE estimates (only applies if se_est=TRUE)
  #inference_type: string -- can be parametric, non-parametric, or jackknife. Parametric is recommended if treatment units are few (approx 40).
  #factor_range: the sequence of unobservable factors to estimate, selected by cross-validation if cross_val=TRUE
  #force_FE: string (unit, time, two-way, none) indicating the type of fixed effects to estimate
  #cross_val: indicate whether to use cross validation to select the optimal number of factors (or the hyperparameter if matrix completion)
  #EM_flag: boolean indicating whether EM aglorithm from (Gobillon and Magnac 2016) is to be used for estimating factors
  #estimator_type: string controlling the estimation method -- either Interactive Fixed Effects "ife" or Matrix Completion ("mc")
  #parallel_boot: boolean for whether parallel computing is to be used for bootstrap/se. (FAILS FOR ME....)
  
  #Higher level description of data_full:
  #rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  #each row should have a treatment indicator (treat_indicator), a period number (time_var), 
  #an individual ID (id_var), and an outcome (outcome_var)
  #for associated with that period-ID combination
  
  
  #Output
  
  #estimate the panel SCM 
  gsynth_agg_te_all_t <- gsynth(Y=outcome_var, D=treat_indicator, data = data_full, index = c(id_var,time_var), X=X_in,
                                se=se_est, nboots=num_boots,inference=inference_type,r=factor_range, force=force_FE, CV=cross_val,
                                EM=EM_flag, estimator = estimator_type, parallel = parallel_boot)
 
  
  
  if(se_est & inference_type=="parametric"){
    #gets all indiv effects for all T
    #renames the multi-d array output so that we can more easily pivot it
    gsynth_indiv_te_series=gsynth_agg_te_all_t$est.ind[, ,1:(gsynth_agg_te_all_t$Ntr)] %>% 
      as.data.frame() %>% mutate(!!as.name(time_var):=1:nrow(.))  %>% clean_names() %>%
      setNames(ifelse(str_count(names(.), "_" )>1, str_replace(names(.),"[_]","" ), names(.) )) 
    
    
    #Pivot the wide dataframe of effects (and CI) per entry with one row per period into a long data set
    #where each period-entry combination has a column for effect, CI.
    #Then, for each period, summarize over all the entries 
    gsynth_series_output=gsynth_indiv_te_series %>%
      pivot_longer(
        -!!as.name(time_var),
        names_to=c(".value", id_var),
        names_sep="_" )  %>% 
      mutate( !!as.name(id_var):=as.numeric(!!as.name(id_var)) ) %>% #here, we manually add on the true target so that we can compute percent change
      inner_join(data_full %>% 
                   group_by(!!as.name(id_var)) %>% mutate(Treatment_Period=length(!!as.name(treat_indicator))-sum(!!as.name(treat_indicator))+1) %>% 
                   ungroup() %>%
                   filter(!is.na(Treatment_Period)) %>%
                   dplyr::select(id_var, time_var, outcome_var, Treatment_Period), by=c(id_var, time_var)
      ) %>% rename(response=outcome_var, point.effect=eff) %>% mutate(point.pred=response-point.effect)
    
    
  }
  else{
    gsynth_series_output= (gsynth_agg_te_all_t$Y.tr-gsynth_agg_te_all_t$Y.ct) %>% as.data.frame() %>% rownames_to_column(var=time_var) %>%
      pivot_longer(
        -!!as.name(time_var),
        names_to=id_var,
        values_to="point.effect"
      ) %>% mutate(!!as.name(id_var):=as.numeric(!!as.name(id_var)),
                   !!as.name(time_var):=as.numeric(!!as.name(time_var))) %>% 
    inner_join(data_full %>% 
                 group_by(!!as.name(id_var)) %>% mutate(Treatment_Period=length(!!as.name(treat_indicator))-sum(!!as.name(treat_indicator))+1) %>% 
                 ungroup() %>%
                 filter(!is.na(Treatment_Period)) %>%
                 dplyr::select(id_var, time_var, outcome_var, Treatment_Period), by=c(id_var, time_var)
    ) %>% rename(response=outcome_var) %>% mutate(point.pred=response-point.effect)
  }
  
  if(!is.null(counterfac_var)){
    gsynth_series_output=gsynth_series_output %>%inner_join(
      data_full %>% 
        dplyr::select(id_var, time_var, counterfac_var), by=c(id_var, time_var)
    )%>% mutate(
      cf_point.effect=(response-!!as.name(counterfac_var)),
      cf_pct.effect=(response/!!as.name(counterfac_var))-1
    )
  }
  
  #add a column with relative (pct) effect
  gsynth_series_output=gsynth_series_output %>% mutate(
    pct.effect=(response/point.pred)-1
      )
  
  return(gsynth_series_output)
}
  
 
  
  
  
  # gsynth_aggregated_te_series[[synth_type]]=gsynth_series_forbootstrap[[synth_type]] %>% group_by(period) %>% summarise(
  #   mean_abs_effect=mean(eff),
  #   mean_abs_lower=mean(cilower),
  #   mean_abs_upper=mean(ciupper),
  #   med_abs_effect=median(eff),
  #   med_abs_lower=median(cilower),
  #   med_abs_upper=median(ciupper),
  #   mean_percent_effect=mean(percent_eff),
  #   med_percent_effect=median(percent_eff)) %>% ungroup() 
  # 
  # names(gsynth_aggregated_te_series)[synth_type]=names(data_list)[synth_type]
  # #to extract each treated units series, go through this, where the treated unit ID is from dimnames(ind_effects)[[3]][i] (will be string)
  # #ind_effects[,,i]
  # #dimnames(ind_effects)[[3]][i]
  # #ppool_placebo <- multisynth(target ~ treatperiod_0, entry, period, nu = 0.5, placebo_df_long)
  # 
  # 
  # #The above Gsynth aggregation takes uncertainty from the Individual Series, which is much higher
  # #Here, we take the estimate of the ATT directly, along with its CI
  # gsynth_att_series[[synth_type]]=gsynth_pred_placebo$est.att %>% as.data.frame() %>% rownames_to_column(var="period") %>% mutate(period=as.numeric(period))
  # names(gsynth_att_series)[synth_type]=names(data_list)[synth_type]






ScWeight <- function(M, target, zeta = 1) {
  
  if (nrow(M) != length(target)) {
    stop("invalid dimensions")
  }
  
  # solve.QP cannot have 0 penalty for quadratic term
  if (zeta == 0) { zeta = 1e-06 }
  # we solve a QP with parameters [weights, imbalance]
  # where we use an equality constraint to impose that
  # imbalance = M * weights - target. Our objective is encoded as
  # zeta*||weights||^2 + || imbalance ||^2 / length(target)
  # = [ weights, imbalance]' * [zeta*I, 0; 0, (1/length(target)) I]
  # * [weights, imbalance] in our call to solve.QP, the parameter
  # Dmat is this block-diagonal matrix, and we pass dvec=0 because we
  # have no linear term
  
  Dmat <- diag(c(rep(zeta, ncol(M)), rep(1 / length(target), nrow(M))))
  dvec <- rep(0, ncol(M) + nrow(M))
  
  # our first nrow(M)+1 constraints are equality constraints
  # the first nrow(M) impose that M*weights - imbalance = target
  # the next imposes that sum(weights)=1
  # and the remaining constraints impose the positivity of our weights
  
  meq <- nrow(M) + 1
  AT <- rbind(cbind(M, diag(1, nrow(M))),
              c(rep(1, ncol(M)), rep(0, nrow(M))),
              cbind(diag(1, ncol(M)), matrix(0, ncol(M), nrow(M))))
  bvec <- c(target, 1, rep(0, ncol(M)))
  soln <- solve.QP(Dmat, dvec, t(AT), bvec, meq = meq)
  gamma <- soln$solution[1 : ncol(M)]
  
  return(gamma)
  
}



#' Gets SDID prediction for one treated advertiser and one post treatment
#' period.
#'
#' @param Y the matrix of control entries and one treated entry over all
#'  time period prior to the treatment period and one post treatment period
#'  for which the counterfactual prediction is made.
#' @param T_0 the treatment period.
#' @param pre.periods Number of periods right before treatment period exlcuded
#'  from sunthetic control.
#' @param post.periods Number of periods past adoption used in counterfactual
#'  prediction. Default = NULL. means all periods after are predicted.
#' @return a scalar estimate of the counterfactual prediction for one treated
#'  entry and the specified post treatment periods.
#' @export

SdidPredict <- function(Y, T_0,  pre.periods,  post.periods, zeta = var(as.numeric(Y))) {
  
  # The unit weights are only estimated once, but time weights are estimated
  # for each period.
  NN <- nrow(Y)
  pre <- Y[NN, ]
  
  if (is.null(post.periods)) {
    post.periods <- ncol(Y) - T_0
  }
  
  end.t <- min(ncol(Y), T_0 + post.periods)
  start.t <- max(T_0 - pre.periods, 1)
  

  omega.weight <- ScWeight(t(Y[-NN, seq_len(start.t - 1)]),  Y[NN, seq_len(start.t - 1)], zeta = zeta)

  for (t in start.t : end.t)  {
    Yt <- Y[, c(seq_len(start.t - 1), t)]
    TT <- ncol(Yt)
    lambda.weight <- ScWeight(Yt[-NN, -TT], Yt[-NN, TT], zeta = zeta)
    SC.transpose.est <- sum(lambda.weight * Yt[NN, -TT])
    SC.est <- sum(omega.weight * Yt[-NN, TT])
    interact.est <- omega.weight %*% Yt[-NN, -TT] %*% lambda.weight
    pre[t] <- SC.est + SC.transpose.est - interact.est
    
  }
  
  return(pre)
  
}



#' Main function for predictions by NSDID.
#'
#' @param y a vector of history values before treatment.
#' @param ct.m matrix with time series entries
#' @param treatperiod treatment period.
#' @param pre.periods number of periods excluded from training before the
#'  treatment period.
#' @param post.periods number of periods for which prediction is done after the
#'  treatment period. When it is NULL, prediction is done for all periods.
#' @param nnsize Nearest neighbour size to be selected for each treated entry.
#'  default value is NULL.
#' @param scale scaling the entries of the matrix by a constant value to help
#'  the optimization problem as it often fails to encompass large values (so
#'  first scale down to smaller values and then after computation scale
#'  up to large values).
#' @param period training period as history before treatperiod - period is
#'  ignored.
#' @return prediction for the given entry.
#' @export

NSDIDPrediction <-function(y, ct.m, treatperiod, pre.periods = 0, post.periods = 20,nnsize = NULL,   scale = 100, period = 52) {
  
  Y_con <- ct.m / scale
  nc <- ncol(Y_con)
  np <- nrow(Y_con)
  Y_pre <- c(y / scale, rep(0, np - treatperiod + 1))
  
  # tv is a linear weight vector. Periods before treatperiod - period use
  # zero weights.
  # nnsize neighbours are identified using weighted distances before the
  # treatment period.
  
  start.period <- max(1, (treatperiod - 1) - period - 1)
  tv <- c(rep(0, start.period), seq_len(np - start.period))
  tmc <- t(matrix(tv, nrow = nc, ncol = length(tv), byrow = TRUE))
  wtY_con <- sqrt(tmc) * Y_con
  wtY_pre <- sqrt(tv) * Y_pre
  
  
  
  # If nnsize is NULL, the number of the neighbours are chosen to be close
  # to the number of periods.
  
  if (is.null(nnsize)) {
    nnsize <- treatperiod
  }
  
  Yr <- t(matrix(wtY_pre, nrow = nc, ncol = np, byrow = TRUE))
  E <- (Yr[seq_len(treatperiod - 1), ] - wtY_con[seq_len(treatperiod - 1), ]) ** 2
  
  ES <- colSums(E)
  ES_order <- order(ES)
  
  Y <- cbind(Y_con[, ES_order[seq_len(min(nnsize, nc))]], Y_pre)
  
  pred <- SdidPredict(t(Y), treatperiod, pre.periods = pre.periods, post.periods = post.periods)
  
  return(pred * scale)
  
}




estimate_scdid_series<-function(data_full, id_var="entry", time_var="period", treat_indicator="treatperiod_0", 
                                outcome_var="target", counterfac_var="counter_factual",
                                pre_SDID = 0, post_SDID = NULL,nn_SDID = NULL,   scale_SDID = 100, period_SDID = 30){
  #Estimates SCDID treatment effects given a long form data set, outputting a dataframe
  #consisting of a series of treatments effects for each id_var by time_var in all  periods
  
  #Args
  #data_full: long-form dataframe with both treated and control entries, sorted by period and entry (ascending).
  #id_var: column name of numeric, unique ID representing the entry (unit) 
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #treat_indicator:column name of binary (0, 1) indicator for whether id_var i in time_var t was treated. Once treated, must always be treated (for now)
  #outcome_var: the y var for the time series
  #counterfac_var: the counterfactual value for the time series, if available (otherwise, NULL)
  #remaining values are from NSDIDPrediction function
  # pre_SDID: number of periods excluded from training before the treatment period.
  #post_SDID: number of periods for which prediction is done after the treatment period. When it is NULL, prediction is done for all periods.
  #nn_SDID: Nearest neighbour size to be selected for each treated entry. default value is NULL.
  #scaleSDID: scaling the entries of the matrix by a constant value to help the optimization problem as it often fails to 
  #encompass large values (so first scale down to smaller values and then after computation scale up to large values).
  #period_SDID: number determining how many pre-treat periods to weight in the estimation (for synthetic control)
  
  #Higher level description of data_full:
  #rows of the df represent period-entry combinations (eg N (total num of entry) rows for period t).
  #each row should have a treatment indicator (treat_indicator), a period number (time_var), 
  #an individual ID (id_var), and an outcome (outcome_var)
  #for associated with that period-ID combination
  
  
  #Output
  #Dataframe containing the full series of the outcome (all T), as well as predicted (counterfactual) outcome
  #and the associated effects by id_Var, time_var
  
  
  #Split the dataset based on whether they are ever treated
  tr.entries  <- data_full %>% filter(!!as.name(treat_indicator)>0) %>% distinct( !!as.name(id_var) ) %>% pull() %>% sort()
  ct.entries  <- setdiff(data_full %>% distinct( !!as.name(id_var) ) %>% pull(), tr.entries)

  
  #create control data frame, with a new id for the sake of ordering observations later
  control_data  <- data_full %>% filter(!!as.name(id_var) %in% ct.entries) 
  #n0, number of control entries, is just the number of unique entries in cd
  n0=control_data %>% distinct(!!as.name(id_var)) %>% nrow()
  
  
  #In the loop, we also want to compute the SCDID Estimates as they require a single observation at a time
  #Because SCDID cannot handle staggered adoption, introduce one treated unit at a time
  #THIS GOES WITHIN Causal Impact For Loop
  treat_data   <- data_full %>% filter(!!as.name(id_var) %in% tr.entries)%>% mutate(new_id=group_indices(., c(!!as.name(id_var)))+n0 ) %>% 
     arrange(!!as.name(time_var), new_id) %>% group_by(new_id) %>%
    mutate(Treatment_Period=length(!!as.name(treat_indicator))-sum(!!as.name(treat_indicator))+1) %>% ungroup()
  
  
  #create the control matrix once, which is an input to SDID estimator
  control_matrix <- spread(control_data %>% dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(outcome_var)),
                           !!as.name(time_var), !!as.name(outcome_var))%>% dplyr::select(-!!as.name(id_var))  %>% as.matrix() %>% t()
  
  
  #for each treated unit, find when it was treated
  list_of_treat_times=treat_data%>%
    split(.[[id_var]]) %>% future_map(~ max(.[["Treatment_Period"]]))
  
  
  #split the treated units into individual vectors
  list_of_treat_data=treat_data%>%
    split(.[[id_var]]) %>% future_map(~ .[[outcome_var]]) %>% future_map(~as.matrix(.) %>% t())
  

   list_of_scdid_series=future_pmap(list(list_of_treat_data, rep(list(control_matrix), length(tr.entries)),list_of_treat_times,
                                        rep(list(pre_SDID),length(tr.entries) ), rep(list(post_SDID),length(tr.entries)), 
                                        rep(list(nn_SDID),length(tr.entries)), rep(list(scale_SDID),length(tr.entries)),
                                        rep(list(period_SDID),length(tr.entries))),
                                   NSDIDPrediction) 
  
  
  # list_of_scdid_series=future_pmap(list(list_of_treat_data, control_matrix,list_of_treat_times,
  #                                       pre_SDID,tr.entries, post_SDID,tr.entries, 
  #                                       nn_SDID,scale_SDID,
  #                                      period_SDID),
  #                                  NSDIDPrediction) 
  #compute the TE by subtracting the matrix of predictions (list_of_scdid_series) from the outcome_var
  #Reformat so that the output is the same as the other functions, namely,
  #each row represents a period-unit combination, with the outcome_var, the prediction, the effect (and the treatment time/indicator)
  
  df_scdid_series=list_of_scdid_series %>% as.data.frame() %>% rownames_to_column(var= time_var) %>%
    mutate(!!as.name(time_var):=as.numeric(!!as.name(time_var))) %>%
    pivot_longer(
     cols= starts_with("X"),
      names_to= "temp_id",
     names_transform=list(temp_id=readr::parse_number),
     values_to="point.pred")  %>% rename(!!as.name(id_var):=temp_id) %>% 
    inner_join(
      treat_data %>% dplyr::select(!!as.name(id_var), !!as.name(time_var), !!as.name(outcome_var), Treatment_Period) %>%
        rename(response=outcome_var), by=c(id_var, time_var)
    ) %>% mutate(point.effect=response-point.pred)
  
  
  if(!is.null(counterfac_var)){
    df_scdid_series=df_scdid_series %>%left_join(
      treat_data %>% 
        dplyr::select(id_var, time_var, counterfac_var), by=c(id_var, time_var)
    ) %>% mutate(
      cf_point.effect=(response-!!as.name(counterfac_var)),
      cf_pct.effect=(response/!!as.name(counterfac_var))-1
    )
  }
  
  #add a column with relative (pct) effect
  df_scdid_series=df_scdid_series %>% mutate(
    pct.effect=(response/point.pred)-1
    
  )

  
  return(df_scdid_series)
  
}
  

#Could add a pure Interactive Fixed Effects Model (interFE) fairly easily






#lots of edits:
#weight constraints to be positive? 
#multiple methods?
ensemble_placebo_weights<-function(method1_estimated_df, method2_estimated_df, method3_estimated_df,
                                   time_var="period", id_var="entry", treat_time_var="Treatment_Period",
                                   pred_var="point.pred", counterfac_var="counter_factual", pos_coef_constr=F){
  
  #Estimates the weights for an ensemble using linear regression on a placebo set
  #This method is, as of now, highly specific -- intended for a placebo data set
  #with fake treatment (but zero treatment effect). Given this structure, we have fit several SCM methods
  #on this placebo set to estimate the TE (hopefully close to zero)
  #This method ignores any pre-treat information and outputs weights on the post-treatment point predictions
  #from each method that yields the closest value to the truth. Thus, we must know the truth to compute these.
  
  
  #Args
  #method1_estimated_df: long-form dataframe of the estimated and counterfactual point predictions for a given method
  #method2_estimated_df, method3_estimated_df distinct methods of the same form.
  #id_var: column name of numeric, unique ID representing the entry (unit) 
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #treat_time_var: string variable for the column indicating when that unit was treated
  #pred_var: string name of the point estimate column in the methodX dataframe
  #counterfac_var: name of var for the counterfactual value for the time series, if available (otherwise, NULL)
 #pos_coef_constr: Boolean flag for whether the weights should be non-negative
  
  #Output
  #weights, as a 3x1 (num methods X 1) vector, from an unconstrained linear reg
  
  combined_methods_df=method1_estimated_df %>%
    dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(treat_time_var), 
                  !!as.name(pred_var), !!as.name(counterfac_var)  ) %>% 
    rename(m1_pred=!!as.name(pred_var)) %>%
    left_join(
      method2_estimated_df %>%
        dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(pred_var) ) %>%
        rename(m2_pred=!!as.name(pred_var)), by=c(id_var, time_var)
    ) %>% left_join(
      method3_estimated_df %>%
        dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(pred_var) ) %>%
        rename(m3_pred=!!as.name(pred_var)), by=c(id_var, time_var)
    )
  
  #Find the weights of the three methods that best fit the data in the post period
  #For now, we do this with simple linear regression and no constraints
  post_treat_combined_df=combined_methods_df %>% filter(!!as.name(time_var)>=!!as.name(treat_time_var))
  
  if(pos_coef_constr==T){
    X=as.matrix(cbind(combined_methods_df %>% pull(m1_pred),
                      combined_methods_df %>% pull(m2_pred),
                      combined_methods_df %>% pull(m3_pred)))
    
    
    Rinv <- solve(chol(t(X) %*% X));
    C <- cbind(rep(1,3), diag(3))
    b <- c(1,rep(0,3))
    
    d <- t(combined_methods_df %>% pull(!!as.name(counterfac_var))) %*% X 
    nn2 = sqrt(norm(d,"2"))
    
    constr_weights_sol=solve.QP(Dmat = Rinv*nn2, factorized = TRUE, dvec = d/(nn2^2), Amat = C, bvec = b, meq = 1)
    weight_vec=constr_weights_sol$solution
  }
  else{
    lm_formula=paste(counterfac_var,"~+offset(m1_pred)+I(m2_pred-m1_pred)+I(m3_pred-m1_pred)+0")
    constrained_linreg_ensemble_weights=lm(as.formula(lm_formula),data=post_treat_combined_df)
    
    weight_vec=c(1-constrained_linreg_ensemble_weights$coefficients[[1]]-constrained_linreg_ensemble_weights$coefficients[[2]],
                 constrained_linreg_ensemble_weights$coefficients[[1]], constrained_linreg_ensemble_weights$coefficients[[2]])
    
  }
  

  

  return(weight_vec)
}




ensembled_predictor<-function(method1_estimated_df, method2_estimated_df, method3_estimated_df, est_weights=c(1/3,1/3,1/3),
                              time_var="period", id_var="entry", treat_time="Treatment_Period",
                              pred_var="point.pred", counterfac_var="counter_factual",
                              outcome_var="response"){
  #Estimates the ensemble using weights and predictions from 3 methods
  
  
  #Args
  #method1_estimated_df: long-form dataframe of the estimated and counterfactual point predictions for a given method
  #method2_estimated_df, method3_estimated_df distinct methods of the same form.
  #est_weights: 3 (num methods) by 1 vector of numeric weights to be placed on the predictions of each method
  #id_var: column name of numeric, unique ID representing the entry (unit) 
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #treat_time_var: string variable for the column indicating when that unit was treated
  #pred_var: string name of the point estimate column in the methodX dataframe
  #counterfac_var: name of var for the counterfactual value for the time series, if available (otherwise, NULL)
  #outcome_var: the original y variable the methods were estimating
  
  #Output
  #Dataframe with the ensemble predictions by id, period for all time, as well as point effects and counterfactual effects
  
  #should work even if counter_fac is null
  ensemble_output=method1_estimated_df %>%
    dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(treat_time),
                  !!as.name(outcome_var)) %>% 
    mutate(
      point.pred=as.vector(as.matrix(
        cbind(method1_estimated_df %>% pull(!!as.name(pred_var)),
              method2_estimated_df %>% pull(!!as.name(pred_var)),
              method3_estimated_df %>% pull(!!as.name(pred_var)))
      ) %*% est_weights
      ) )%>% 
    mutate(point.effect=response-point.pred)
  
  if(!is.null(counterfac_var)){
    ensemble_output=ensemble_output %>%inner_join(
      method1_estimated_df %>% 
        dplyr::select(id_var, time_var, counterfac_var), by=c(id_var, time_var)
    )%>% mutate(
      cf_point.effect=(response-!!as.name(counterfac_var)),
      cf_pct.effect=(response/!!as.name(counterfac_var))-1
    )
  }
  
  return(ensemble_output)
}









# #Uncertainty Bootstrap
# compute_ci_bounds_jackknife<-function(estimated_series_df, time_var="period",treat_period_var="Treatment_Period", boot_var="point.effect", alpha_ci=0.95){
#   
#   effect_series_postonly=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
#     filter(post_period_t>=0)
#   #http://www.math.montana.edu/jobo/thainp/jack.pdf
#   #according to this, we want to compute the jackknife() of our statistic,
#   #and then remove the bias ($jack.bias) from the statistic on the dataset to get the
#   #jackknife value
#   theta_quant_low<-function(temp_data_input, a=alpha_ci){
#     quantile(temp_data_input, probs=c( 1-a))
#   }
#   
#   theta_quant_high<-function(temp_data_input, a=alpha_ci){
#     quantile(temp_data_input, probs=c( a))
#   }
#   
#   #Parallelized computation of the pre-treatment range, by id_var
#    low_quantile_bias_by_postperiod=effect_series_postonly%>%
#      split(.[["post_period_t"]]) %>% future_map(~ .[[boot_var]]) %>% future_map(~jackknife(.,theta_quant_low)) %>% future_map(~ .$"jack.bias") %>%
#      do.call( bind_rows, . ) %>% pull()
#    
#    high_quantile_bias_by_postperiod=effect_series_postonly%>%
#      split(.[["post_period_t"]]) %>% future_map(~ .[[boot_var]]) %>% future_map(~jackknife(.,theta_quant_high)) %>% future_map(~ .$"jack.bias") %>%
#      do.call( bind_rows, . ) %>% pull()
#    
#    orig_lowquantiles_by_postperiod=effect_series_postonly%>%
#      split(.[["post_period_t"]]) %>% future_map(~ .[[boot_var]]) %>% future_map(~theta_quant_low(.)) %>% do.call( bind_rows, . ) %>% pull()
#    
#    orig_highquantiles_by_postperiod=effect_series_postonly%>%
#      split(.[["post_period_t"]]) %>% future_map(~ .[[boot_var]]) %>% future_map(~theta_quant_high(.)) %>% do.call( bind_rows, . ) %>% pull()
#   
#   debiased_quantiles_df=tibble("post_treat_t"=effect_series_postonly %>% distinct(post_period_t) %>% pull(),
#                                    "jackknife_upper_ci"=orig_highquantiles_by_postperiod -high_quantile_bias_by_postperiod,
#                                    "jackknife_lower_ci"=orig_lowquantiles_by_postperiod -low_quantile_bias_by_postperiod)
#   
#   return(debiased_quantiles_df)
# }



compute_tot_se_jackknife<-function(estimated_series_df, time_var="period",treat_period_var="Treatment_Period", 
                                   boot_var="point.effect", stat_in="mean", alpha_ci=0.95,
                                   compute_cf_eff=T, counterfac_eff=NULL){
  #Computes jackknife estimates (from resample package) of the treatment effect on the treated, by period
  
  #Args
  #estimated_series_df: long-form dataframe with the estimated effects, aka output from one of the methods above.
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #treat_period_var:column name of the period for which the unit's treatment starts (assumes it ends at T max)
  #boot_var: string indicating the effect we wish to bootstrap, either "point.effect" or "pct.effect" works
  #stat_in: string indicating the statistic which we hope to get a bootstrap mean and se for, either "mean" or "median" typical
  # alpha_ci: number between 0 and 1, indicating the confidence interval desires
  #compute_cf_eff: boolean flag for whether a column of counterfactual effects should be appended to output
  #counterfac_eff: string variable name of the counterfactual effect to estimate, typically "cf_point.effect or "cf_pct.effect" from above
  
  #Output
  #Tibble containing, by post_treat_period, the sample stat and bootstrapped mean of boot_var, as well as upper and lower bounds
  #that make up the alpha_ci*100% confidence interval
  
  effect_series_postonly=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
    filter(post_period_t>=0)

  if(compute_cf_eff){
    
    if(is.null(counterfac_eff)){
      counterfac_eff=paste("cf_", boot_var, sep="")
    }
    
    cf_tot_df=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
      filter(post_period_t>=0)
    if(tolower(stat_in)=="mean"){
      cf_tot_df=cf_tot_df %>% group_by(post_period_t) %>%
      summarise(mean_cf_tot=mean(!!as.name(counterfac_eff))) %>% ungroup()
    }
    if(tolower(stat_in)=="median"){
      cf_tot_df=cf_tot_df %>% group_by(post_period_t) %>%
        summarise(median_cf_tot=median(!!as.name(counterfac_eff))) %>% ungroup()
    }
  }
  
  #store the bootstrap metric by post treatment period
  if(tolower(stat_in)=="mean"){
    #Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp=effect_series_postonly%>%
      split(.[["post_period_t"]]) %>% future_map(~ .[[boot_var]]) %>% future_map(~jackknife(.,mean))
  }
  if(tolower(stat_in)=="median"){
    #Parallelized computation of the pre-treatment range, by id_var
    jackknife_comp=effect_series_postonly%>%
      split(.[["post_period_t"]]) %>% future_map(~ .[[boot_var]]) %>% future_map(~jackknife(.,median))
  }
  

  statname=paste("jackknife_", stat_in,"_tot", sep="")
  jackknife_comp_ci= future_pmap(list(jackknife_comp, list(probs=c(0.5-alpha_ci/2,0.5+alpha_ci/2 ))), CI.t) %>% do.call(rbind, .) %>%
    as_tibble() %>% mutate(post_period_t=effect_series_postonly %>% distinct(post_period_t) %>% pull())  %>% mutate(
      !!as.name(statname):= jackknife_comp%>% future_map(~.[["stats"]]) %>% future_map(~.[["Mean"]]) %>% 
      unlist()
  )  

  names(jackknife_comp_ci)[1:2]=c("jackknife_lb_tot", "jackknife_ub_tot")

  #Manual version of the above code
  # jackknife_mean_bias= jackknife_comp%>% future_map(~ .$"stats") %>% future_map(~ .$"Mean") %>%
  #   unlist() %>% tibble("de_biased_mean"=.,  "post_period_t"=effect_series_postonly %>% distinct(post_period_t) %>% pull())
  # 
  # jackknife_se=jackknife_comp%>% future_map(~ .$"stats") %>% future_map(~ .$"SE") %>%
  #   unlist() %>% tibble("se"=.,"post_period_t"=effect_series_postonly %>% distinct(post_period_t) %>% pull())
  
  if(tolower(stat_in)=="mean"){
    stat_tot_sample=effect_series_postonly%>% group_by(post_period_t) %>%
      summarise(sample_mean_tot=mean(!!as.name(boot_var)),
                treated_n=n()) %>% ungroup()
  }
  
  if(tolower(stat_in)=="median"){
    stat_tot_sample=effect_series_postonly%>% group_by(post_period_t) %>%
      summarise(sample_median_tot=median(!!as.name(boot_var)),
                treated_n=n()) %>% ungroup()
  }
  
  if(compute_cf_eff){
    jackknife_ci_by_postperiod=stat_tot_sample %>% left_join(
      jackknife_comp_ci, by="post_period_t"
    )  %>% left_join(
      cf_tot_df, by="post_period_t"
    )
  }
  if(!compute_cf_eff){
    jackknife_ci_by_postperiod=stat_tot_sample %>% left_join(
      jackknife_comp_ci, by="post_period_t"
    ) 
  }

  
  return(jackknife_ci_by_postperiod)
}
  


#Uncertainty Bootstrap
compute_ci_bounds_bootstrap<-function(estimated_series_df, time_var="period",treat_period_var="Treatment_Period", 
                                      boot_var="point.effect", stat_in="mean", alpha_ci=0.95,
                                      compute_cf_eff=T,counterfac_eff="cf_point.effect",
                                      nboots=10000){
  #Computes bootstrap estimates (from resample package) of the treatment effect on the treated, by period
  
  #Args
  #estimated_series_df: long-form dataframe with the estimated effects, aka output from one of the methods above.
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #treat_period_var:column name of the period for which the unit's treatment starts (assumes it ends at T max)
  #boot_var: string indicating the effect we wish to bootstrap, either "point.effect" or "pct.effect" works
  #stat_in: string indicating the statistic which we hope to get a bootstrap mean and se for, either "mean" or "median" typical
  # alpha_ci: number between 0 and 1, indicating the confidence interval desires
  #compute_cf_eff: boolean flag for whether a column of counterfactual effects should be appended to output
  #counterfac_eff: string variable name of the counterfactual effect to estimate, typically "cf_point.effect or "cf_pct.effect" from above
  #nboots, number of bootstrap samples
  
  #Output
  #Tibble containing, by post_treat_period, the sample stat and bootstrapped mean of boot_var, as well as upper and lower bounds
  #that make up the alpha_ci*100% confidence interval
  
  effect_series_postonly=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
    filter(post_period_t>=0)
  
  if(compute_cf_eff){
    cf_tot_df=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
      filter(post_period_t>=0)
    if(tolower(stat_in)=="mean"){
      cf_tot_df=cf_tot_df %>% group_by(post_period_t) %>%
        summarise(mean_cf_tot=mean(!!as.name(counterfac_eff))) %>% ungroup()
    }
    if(tolower(stat_in)=="median"){
      cf_tot_df=cf_tot_df %>% group_by(post_period_t) %>%
        summarise(median_cf_tot=median(!!as.name(counterfac_eff))) %>% ungroup()
    }
  }
  
  
  if(tolower(stat_in)=="mean"){
    #Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp=effect_series_postonly%>%
      split(.[["post_period_t"]]) %>% future_map(~ .[[boot_var]]) %>% future_map(~bootstrap(.,mean, R=nboots))
  }
  if(tolower(stat_in)=="median"){
    #Parallelized computation of the pre-treatment range, by id_var
    bootstrap_comp=effect_series_postonly%>%
      split(.[["post_period_t"]]) %>% future_map(~ .[[boot_var]]) %>% future_map(~bootstrap(.,median, R=nboots))
  }
  
  statname=paste("bootstrap_", stat_in,"_tot" ,sep="")
  bootstrap_comp_ci= future_pmap(list(bootstrap_comp, list(probs=c(0.5-alpha_ci/2,0.5+alpha_ci/2 ))), CI.t) %>% do.call(rbind, .) %>%
    as_tibble() %>% mutate(post_period_t=effect_series_postonly %>% distinct(post_period_t) %>% pull()) %>%
    mutate(
      !!as.name(statname):= bootstrap_comp%>% future_map(~.[["stats"]]) %>% future_map(~.[["Mean"]]) %>% 
        unlist()
    )  
  names(bootstrap_comp_ci)[1:2]=c("bootstrap_lb_tot", "bootstrap_ub_tot")
  
  if(tolower(stat_in)=="mean"){
    stat_tot_sample=effect_series_postonly%>% group_by(post_period_t) %>%
    summarise(sample_mean_tot=mean(!!as.name(boot_var)),
              treated_n=n()) %>% ungroup()
  }
  if(tolower(stat_in)=="median"){
    stat_tot_sample=effect_series_postonly%>% group_by(post_period_t) %>%
      summarise(sample_median_tot=median(!!as.name(boot_var)),
                treated_n=n()) %>% ungroup()
  }
  
  
  if(compute_cf_eff){
    bootstrap_ci_by_postperiod=stat_tot_sample %>% left_join(
      bootstrap_comp_ci, by="post_period_t"
    )   %>% left_join(
      cf_tot_df, by="post_period_t"
    )
  }
  if(!compute_cf_eff){
    bootstrap_ci_by_postperiod=stat_tot_sample %>% left_join(
      bootstrap_comp_ci, by="post_period_t"
    ) 
  }
  
  
  
  return(bootstrap_ci_by_postperiod)
}
  



#################################################################
#Plots
#################################################################



#print PDF Likely does not work with future_map???
create_gap_ci_plot<-function(bootstrapped_effects_df, time_var="post_period_t", effect_var="jackknife_mean_tot",
                             upper_ci="jackknife_ub_tot",  lower_ci="jackknife_lb_tot", 
                             cf_plot=T, cf_var="mean_cf_tot", 
                             print_to_pdf=NULL, plot_title=NULL, plot_x_lab=NULL, plot_y_lab=NULL){
  
  #Computes bootstrap estimates (from resample package) of the treatment effect on the treated, by period
  
  #Args
  #bootstrapped_effects_df: data on the bootstrapped (and potentially counterfactual) effects for each post treat period, 
  # from compute_ci_bounds_bootstrap or jackknife version
  #time_var:column name of numeric period number indicating the post treatment time period, in increasing order (eg 0 is the first time)
  #effect_var: string indicating the effect we wish to plot, either "jackknife_mean_tot" or "jackknife_median_tot" works, depending on input
  #upper_ci: string indicating thename of the upper CI variable in bootstrapped_effects_df
  #lower_ci: string indicating thename of the lower CI variable in bootstrapped_effects_df
  #cf_plot: boolean flag for whether a counterfactual effect should be plotted
  #cf_var: string variable name of the counterfactual effect to plot, typically "mean_cf_tot or "median_cf_tot" from bootstrapped_effects_df
  #print_to_pdf: string, file path where the pdfs should be printed. 
  #plot_title: title of the plot to be printed
  #plot_x_lab: title of the x axis label
  #plot_y_lab: title of the y axis label
  
  #Output
  #gap plot, potentially with counterfactual predictions
  
  if(!is.null(print_to_pdf)){
    pdf(print_to_pdf)
  }
  
  if(cf_plot){
    plot_out=bootstrapped_effects_df%>% 
      ggplot(aes(x=!!as.name(time_var) , y=!!as.name(effect_var), color="Estimate")) + 
      geom_line()+ geom_ribbon(aes(ymin=!!as.name(lower_ci),ymax=!!as.name(upper_ci)),alpha=0.3, color=NA)+ 
      geom_line(aes(y=!!as.name(cf_var), color="True" ))+
      ggtitle(ifelse(is.null(plot_title), paste("Gap Plot of", effect_var), plot_title ))+
      labs(x=ifelse(is.null(plot_x_lab), time_var, plot_x_lab ), 
           y=ifelse(is.null(plot_y_lab), effect_var, plot_y_lab ))+
      scale_colour_manual(name = "", 
                          values =c('Estimate'='black','True'='red'), labels = c('Estimate','True'))
    
      
  }
  
  if(!cf_plot){
    plot_out=bootstrapped_effects_df%>% 
      ggplot(aes(x=!!as.name(time_var) , y=!!as.name(effect_var))) + 
      geom_line()+ geom_ribbon(aes(ymin=!!as.name(lower_ci),ymax=!!as.name(upper_ci)),alpha=0.3)+ 
    ggtitle(ifelse(is.null(plot_title), paste("Gap Plot of", effect_var), plot_title ))+
      labs(x=ifelse(is.null(plot_x_lab), time_var, plot_x_lab ), 
           y=ifelse(is.null(plot_y_lab), effect_var, plot_y_lab ))
    
  }
  
  
  if(!is.null(print_to_pdf)){
    print(plot_out)
    dev.off()
  }
  return(plot_out)
}



#################################################################
#Metrics
#################################################################

compute_avg_metric_per_t<-function(estimated_series_df, time_var="period", outcome_var="response", prediction_var="point.pred",
                                   counterfac_var="counter_factual", treat_period_var="Treatment_Period",
                                   metric_str="mae", pct_eff_flag=F){
  
  #Create the gap plot for each of the post treatment periods
  
  #Args
  #estimated_series_df: long-form dataframe with predicted counterfactual, actual counterfactual, period, and Treatment Period by ID.
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #outcome_var: the true y var for the time series
  #prediction_var: predicted counterfactual by time and id
  #counterfac_var: true counterfactual by time and id
  #treat_period_var: str name of the column which indicates the treatment period for each id_var (at each time_var)
  #metric_str: string indicating desired metric. can be mae, mse, or both
  #pct_eff_flag: binary flag for whether the percent error should be computed
  
  #Output
  #Dataframe containing number of rows equal to the longest post-treat period
  #and for each post treat period (from 0 -- time of treat, to max), an average of the metric for all observations
  #that experienced that particular post treat period (eg treated in t=5, total T=15 means 10 post treat periods)
  
  #absolute error computation
  if(pct_eff_flag==FALSE){
    if(tolower(metric_str)=="mse"){
      avg_metric_series=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_error_sq=(!!as.name(prediction_var)-!!as.name(counterfac_var))^2) %>% 
        group_by(post_period_t) %>% summarise(MSE_TE_byT=mean(indiv_error_sq))
    }
    
    if(tolower(metric_str)=="mae"){
      avg_metric_series=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error=abs(!!as.name(prediction_var)-!!as.name(counterfac_var))) %>% 
        group_by(post_period_t) %>% summarise(MAE_TE_byT=mean(indiv_t_error))
    }
    
    
    if(tolower(metric_str)=="both"){
      avg_metric_series=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error=abs(!!as.name(prediction_var)-!!as.name(counterfac_var)),
               indiv_error_sq=(!!as.name(prediction_var)-!!as.name(counterfac_var))^2) %>% 
        group_by(post_period_t) %>% summarise(MAE_TE_byT=mean(indiv_t_error),
                                              MSE_TE_byT=mean(indiv_error_sq))
    }
  }
  
  #percent error computation
  if(pct_eff_flag){
    if(tolower(metric_str)=="mse"){
      avg_metric_series=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueT- EstTE=(Target/counterFac)-1-(Target/Pred-1)
        mutate(indiv_error_sq=(!!as.name(outcome_var)/!!as.name(counterfac_var) -  (!!as.name(outcome_var)/!!as.name(prediction_var)) )^2) %>% 
        group_by(post_period_t) %>% summarise(MSE_TE_byT=mean(indiv_error_sq))
    }
    
    if(tolower(metric_str)=="mae"){
      avg_metric_series=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error=abs(!!as.name(outcome_var)/!!as.name(counterfac_var) -  (!!as.name(outcome_var)/!!as.name(prediction_var)) )) %>% 
        group_by(post_period_t) %>% summarise(MAE_TE_byT=mean(indiv_t_error))
    }
    
    
    if(tolower(metric_str)=="both"){
      avg_metric_series=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error=abs(!!as.name(outcome_var)/!!as.name(counterfac_var) -  (!!as.name(outcome_var)/!!as.name(prediction_var)) ),
               indiv_error_sq=(!!as.name(outcome_var)/!!as.name(counterfac_var) -  (!!as.name(outcome_var)/!!as.name(prediction_var)))^2) %>% 
        group_by(post_period_t) %>% summarise(MAE_TE_byT=mean(indiv_t_error),
                                              MSE_TE_byT=mean(indiv_error_sq))
    }
  }
  
  
  

  return(avg_metric_series)
}



compute_avg_metric<-function(estimated_series_df, time_var="period", outcome_var="response", prediction_var="point.pred",
                                   counterfac_var="counter_factual", treat_period_var="Treatment_Period",
                                   metric_str="mae", pct_eff_flag=F){
  #Estimates the relevant metric for the predicted Treatment Effect by comparing the predicted and true counterfactuals
  
  #Args
  #estimated_series_df: long-form dataframe with predicted counterfactual, actual counterfactual, period, and Treatment Period by ID.
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #outcome_var: the true y var for the time series
  #prediction_var: predicted counterfactual by time and id
  #counterfac_var: true counterfactual by time and id
  #treat_period_var: str name of the column which indicates the treatment period for each id_var (at each time_var)
  #metric_str: string indicating desired metric. can be mae, mse, or both
  #pct_eff_flag: binary flag for whether the percent error should be computed
  
  #Output
  #tibble containing an average of the metric for all observations, over all T
  
  #absolute error computation
  if(pct_eff_flag==FALSE){
    if(tolower(metric_str)=="mse"){
      avg_metric=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_error_sq=(!!as.name(prediction_var)-!!as.name(counterfac_var))^2) %>% 
        summarise(MSE_TE=mean(indiv_error_sq))
    }
    
    if(tolower(metric_str)=="mae"){
      avg_metric=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error=abs(!!as.name(prediction_var)-!!as.name(counterfac_var))) %>% 
        summarise(MAE_TE=mean(indiv_t_error))
    }
    
    
    if(tolower(metric_str)=="both"){
      avg_metric=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error=abs(!!as.name(prediction_var)-!!as.name(counterfac_var)),
               indiv_error_sq=(!!as.name(prediction_var)-!!as.name(counterfac_var))^2) %>% 
         summarise(MAE_TE=mean(indiv_t_error),
                   MSE_TE=mean(indiv_error_sq))
    }
  }
  
  #percent error computation
  if(pct_eff_flag){
    if(tolower(metric_str)=="mse"){
      avg_metric=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueT- EstTE=(Target/counterFac)-1-(Target/Pred-1)
        mutate(indiv_error_sq=(!!as.name(outcome_var)/!!as.name(counterfac_var) -  (!!as.name(outcome_var)/!!as.name(prediction_var)) )^2) %>% 
        summarise(MSE_TE=mean(indiv_error_sq))
    }
    
    if(tolower(metric_str)=="mae"){
      avg_metric=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error=abs(!!as.name(outcome_var)/!!as.name(counterfac_var) -  (!!as.name(outcome_var)/!!as.name(prediction_var)) )) %>% 
        summarise(MAE_TE=mean(indiv_t_error))
    }
    
    
    if(tolower(metric_str)=="both"){
      avg_metric=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueTE-EstTE=Target-counterFac-(Target-Pred)
        mutate(indiv_t_error=abs(!!as.name(outcome_var)/!!as.name(counterfac_var) -  (!!as.name(outcome_var)/!!as.name(prediction_var)) ),
               indiv_error_sq=(!!as.name(outcome_var)/!!as.name(counterfac_var) -  (!!as.name(outcome_var)/!!as.name(prediction_var)))^2) %>% 
        summarise(MAE_TE=mean(indiv_t_error),
                  MSE_TE=mean(indiv_error_sq))
    }
  }
  
  
  
  
  return(avg_metric)
}


compute_tot_coverage<-function(bootstrapped_series_df, time_var="post_period_t", ub_var="jackknife_ub_tot",
                               lb_var="jackknife_lb_tot", counterfac_var="mean_cf_tot"){
  
  coverage_out=bootstrapped_series_df %>% mutate(
    contains_truth=as.numeric(!!as.name(counterfac_var)<=!!as.name(ub_var) &
                                !!as.name(counterfac_var)>=!!as.name(lb_var))
    ) %>% pull(contains_truth) %>% mean()
    
  return(coverage_out)
}




plot_tot_bias_per_t<-function(estimated_series_df, time_var="period", id_var="entry", outcome_var="response", prediction_var="point.pred",
                                   counterfac_var="counter_factual", treat_period_var="Treatment_Period",
                                   max_post_t=12, pct_eff_flag=F, plot_title=NULL){
  
  #Create the gap plot for each of the post treatment periods
  
  #Args
  #estimated_series_df: long-form dataframe with predicted counterfactual, actual counterfactual, period, and Treatment Period by ID.
  #time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  #id_var
  #outcome_var: the true y var for the time series
  #prediction_var: predicted counterfactual by time and id
  #counterfac_var: true counterfactual by time and id
  #treat_period_var: str name of the column which indicates the treatment period for each id_var (at each time_var)
  #max_post_t: numeric variable, the number of time periods over which we are averaging per unit
  #pct_eff_flag: binary flag for whether the percent error should be computed
  
  #Output
  #Dataframe containing number of rows equal to the longest post-treat period
  #and for each post treat period (from 0 -- time of treat, to max), an average of the metric for all observations
  #that experienced that particular post treat period (eg treated in t=5, total T=15 means 10 post treat periods)
  
  #absolute error computation
  if(pct_eff_flag==FALSE){

    avg_bias_series=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0, post_period_t<=max_post_t) %>%
      group_by(!!as.name(id_var)) %>%
      summarise(bias_per_id=mean(!!as.name(prediction_var) - !!as.name(counterfac_var))) %>% 
      ggplot(aes(x=bias_per_id)) + geom_density(fill="blue", alpha=0.4)+theme_bw()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      geom_vline(xintercept = 0, linetype="dashed", color = "red")+scale_x_continuous(name="ATT by ID")+
      ggtitle(ifelse(is.null(plot_title), paste("ToT Bias by Unit, averaged over", max_post_t, "post-treat Periods"), plot_title ))
  }
  
  #####FIXXXXXXX
  #percent error computation
  if(pct_eff_flag){
   
      avg_bias_series=estimated_series_df %>% mutate(post_period_t=!!as.name(time_var)-!!as.name(treat_period_var)) %>%
        filter(post_period_t>=0) %>% #note that MSE of TE: TrueT- EstTE=(Target/counterFac)-1-(Target/Pred-1)
        mutate(indiv_bias=(!!as.name(outcome_var)/!!as.name(counterfac_var) -  (!!as.name(outcome_var)/!!as.name(prediction_var)) )) %>% 
        group_by(post_period_t) %>% summarise(Bias_byT=mean(indiv_bias))
    
    
  }
  
  
  
  
  return(avg_bias_series)
}


