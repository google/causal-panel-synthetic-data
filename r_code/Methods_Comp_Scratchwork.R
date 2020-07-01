#Methods comparison scratch

#To compute things like MPSE etc, we know the true series and we know the indiv ATT over each time period (except in SCDID)
#Thus, we can just calculate the "True Y" - "Predicted_TE ="Synthetic Y", and compare the errors
#But, since we know true TE=0, we can also just take the function of the TE to be the MSE, ie sqrt( sum( (hat(TE)-0)**2))

#Will want to loop through this -- maybe by period create a df and then a df with Type (SCDID, CI, Gs) and Error (MSE, MAE)

#Post-Placebo-Treatment MSE

#Gsynth
#First, we do an aggregate over all entries and all periods
gsynth_series_forbootstrap[[1]] %>% filter(period>=0) %>%
  mutate(indiv_t_error_sq=(eff)^2) %>% summarise(MSE_post=mean(indiv_t_error_sq))
#MAE
gsynth_series_forbootstrap[[1]] %>% filter(period>=0) %>%
  mutate(indiv_t_error=abs(eff)) %>% summarise(MAE_post=mean(indiv_t_error))

#then, we break it down by period to see if things get worse farther out
gsynth_series_forbootstrap[[1]] %>% filter(period>=0) %>%
  mutate(indiv_t_error_sq=(eff)^2) %>% 
  group_by(period) %>% summarise(MSE_post_byt=mean(indiv_t_error_sq))

gsynth_series_forbootstrap[[1]] %>% filter(period>=0) %>%
  mutate(indiv_t_error=abs(eff)) %>% 
  group_by(period) %>% summarise(MAE_post_byt=mean(indiv_t_error))


#CausalImpact
causalimpact_series_forbootstrap[[1]]%>%mutate(indiv_t_error_sq=(point.effect)^2) %>%
  summarise(MSE_post=mean(indiv_t_error_sq))

#MAE
causalimpact_series_forbootstrap[[1]] %>%  mutate(indiv_t_error=abs(point.effect)) %>% 
  summarise(MAE_post=mean(indiv_t_error))

#then, we break it down by period to see if things get worse farther out
causalimpact_series_forbootstrap[[1]] %>% mutate(indiv_t_error_sq=(point.effect)^2) %>% 
  group_by(period) %>% summarise(MSE_post_byt=mean(indiv_t_error_sq))

causalimpact_series_forbootstrap[[1]] %>% mutate(indiv_t_error=abs(point.effect)) %>% 
  group_by(period) %>% summarise(MAE_post_byt=mean(indiv_t_error))


#SCDID (there is no By Period for this)
#ISSUE: not apples to apples comparison: summing over total errors of ATT.
#Would be (more closely) equivalent to taking the per-entry mean ATT in CI and Gsyn and 
#then computing MSE/MAE
#Not sure if this is still a fair comparison, but its the correct order of magnitude

scdid_att_placebo_forbootstrap[[1]]%>%mutate(indiv_t_error_sq=abs(att)^2) %>%
  summarise(MSE_post=mean(indiv_t_error_sq))

gsynth_series_forbootstrap[[1]] %>% filter(period>=0) %>% group_by(entry) %>% 
  summarise(ATT_post_idagg=mean(eff)) %>%
  pull(ATT_post_idagg) %>% sapply(., function(x){ x^2}) %>% mean()

causalimpact_series_forbootstrap[[1]] %>% filter(period>=0) %>% group_by(entry) %>% 
  summarise(ATT_post_idagg=mean(point.effect)) %>%
  pull(ATT_post_idagg) %>% sapply(., function(x){ x^2}) %>% mean()



scdid_att_placebo_forbootstrap[[1]] %>%  mutate(indiv_t_error=abs(att)) %>% 
  summarise(MAE_post=mean(indiv_t_error))

gsynth_series_forbootstrap[[1]] %>% filter(period>=0) %>% group_by(entry) %>% 
  summarise(ATT_post_idagg=mean(eff)) %>%
  pull(ATT_post_idagg) %>% abs() %>% mean()

causalimpact_series_forbootstrap[[1]] %>% filter(period>=0) %>% group_by(entry) %>% 
  summarise(ATT_post_idagg=mean(point.effect)) %>%
  pull(ATT_post_idagg) %>% abs() %>% mean()








#Why is the NSDID counterfactual always exactly right? Just not stored?
#WE IMPOSE EQUALITY CONSTRAINTS IS WHY

ScWeight_test <- function(M, target, zeta = 1) {
  
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

SdidPredict_test <- function(Y, T_0,  pre.periods,  post.periods, zeta = var(as.numeric(Y))) {
  
  # The unit weights are only estimated once, but time weights are estimated
  # for each period.
  NN <- nrow(Y)
  pre <- Y[NN, ]
  
  print(Y)
  if (is.null(post.periods)) {
    post.periods <- ncol(Y) - T_0
  }
  
  end.t <- min(ncol(Y), T_0 + post.periods)
  start.t <- max(T_0 - pre.periods, 1)
  
  print(start.t)
  print(end.t)
  omega.weight <- ScWeight_test(t(Y[-NN, seq_len(start.t - 1)]),  Y[NN, seq_len(start.t - 1)], zeta = zeta)

  for (t in start.t : end.t)  {
    Yt <- Y[, c(seq_len(start.t - 1), t)]
    TT <- ncol(Yt)
    lambda.weight <- ScWeight_test(Yt[-NN, -TT], Yt[-NN, TT], zeta = zeta)
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

NSDIDPrediction_test <-function(y, ct.m, treatperiod, pre.periods = 0, post.periods = 20,nnsize = NULL,   scale = 100, period = 52) {
  
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
  
  pred <- SdidPredict_test(t(Y), treatperiod, pre.periods = pre.periods, post.periods = post.periods)
  
  return(pred * scale)
  
}


test_scdid_series<-function(data_full, id_var="entry", time_var="period", treat_indicator="treatperiod_0", 
                                outcome_var="target", counterfac_var="counter_factual",
                                pre_SDID = 0, post_SDID = NULL,nn_SDID = NULL,   scale_SDID = 100, period_SDID = 30){

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
                                        rep(list(60),length(tr.entries))),
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
  


inp_to_nsdid=test_scdid_series(data_list[[1]], id_var="entry",
                 time_var="period",
                 (treat_indicator="treatperiod_0"), 
                 (outcome_var="target"),
                 (counterfac_var="counter_factual"),
                 (pre_SDID = 0),
                 (post_SDID = NULL),
                 (nn_SDID = NULL),
                 (scale_SDID = 100),
                 (period_SDID = 30)) 




NSDIDPrediction_test(inp_to_nsdid[[1]]$`1`, inp_to_nsdid[[2]], inp_to_nsdid[[3]],
                inp_to_nsdid[[3]]+1, inp_to_nsdid[[5]], inp_to_nsdid[[6]],
                inp_to_nsdid[[7]], inp_to_nsdid[[8]])











#Combining the model predictions via ensemble

#lots of edits:
#weight constraints to be positive? 
#multiple methods?
ensemble_placebo_weights<-function(method1_est, method2_est, method3_est,
                                      time_var="period", id_var="entry", treat_time="Treatment_Period",
                                      pred_var="point.pred", counterfac_var="counter_factual"){
  
  combined_methods_df=method1_est %>%
    dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(treat_time), 
                  !!as.name(pred_var), !!as.name(counterfac_var)  ) %>% 
    rename(m1_pred=!!as.name(pred_var)) %>%
    left_join(
      method2_est %>%
        dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(pred_var) ) %>%
        rename(m2_pred=!!as.name(pred_var)), by=c(id_var, time_var)
    ) %>% left_join(
      method3_est %>%
        dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(pred_var) ) %>%
        rename(m3_pred=!!as.name(pred_var)), by=c(id_var, time_var)
    )
  
  #Find the weights of the three methods that best fit the data in the post period
  #For now, we do this with simple linear regression and no constraints
  post_treat_combined_df=combined_methods_df %>% filter(!!as.name(time_var)>=!!as.name(treat_time))
  
  lm_formula=paste(counterfac_var,"~+m1_pred+m2_pred+m3_pred+0")
  unconstrained_linreg_ensemble_weights=lm(as.formula(lm_formula),data=post_treat_combined_df)
  
  
  
  
  return(unconstrained_linreg_ensemble_weights$coefficients)
  
  
}


temp_reg_in=ensemble_placebo_weights(gsynth_ife_AA[[1]], scdid_results_AA[[1]], causalimpact_results_AA[[1]])



ensembled_predictor<-function(method1_est, method2_est, method3_est, est_weights=c(1/3,1/3,1/3),
                              time_var="period", id_var="entry", treat_time="Treatment_Period",
                              pred_var="point.pred", counterfac_var="counter_factual",
                              outcome_var="response"){
  
  #should work even if counter_fac is null
  ensemble_output=method1_est %>%
    dplyr::select(!!as.name(time_var), !!as.name(id_var), !!as.name(treat_time),
                 !!as.name(outcome_var)) %>% 
    mutate(
      point.pred=as.vector(as.matrix(
        cbind(method1_est %>% pull(!!as.name(pred_var)),
              method2_est %>% pull(!!as.name(pred_var)),
              method3_est %>% pull(!!as.name(pred_var)))
      ) %*% est_weights
    ) )%>% 
    mutate(point.effect=response-point.pred)
  
  if(!is.null(counterfac_var)){
    ensemble_output=ensemble_output %>%inner_join(
      method1_est %>% 
        dplyr::select(id_var, time_var, counterfac_var), by=c(id_var, time_var)
    )%>% mutate(
      cf_point.effect=(response-!!as.name(counterfac_var)),
      cf_pct.effect=(response/!!as.name(counterfac_var))-1
    )
  }
  
  return(ensemble_output)
}


temp_ens=ensembled_predictor(gsynth_ife_AB[[1]], scdid_results_AB[[1]], causalimpact_results_AB[[1]],temp_reg_in)

error_ensem=compute_tot_se_jackknife(temp_ens)
create_gap_ci_plot(error_ensem)
compute_avg_metric(temp_ens, metric_str = "both")






