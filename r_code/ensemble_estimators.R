pacman::p_load(dplyr, ggplot2, quadprog, purrr, furrr, tidyr, glue, tibble)

reduce_pred_list <- function(pred_list, time_var = "period", 
                                 id_var = "entry", treat_time_var = "Treatment_Period",
                                 pred_var = "point.pred",
                             outcome_var="response",
                             counterfac_var = "counter_factual"){
  tib_out=pred_list %>% 
    purrr::reduce(dplyr::left_join, 
                  by=c(id_var, time_var, treat_time_var,outcome_var)) %>%
    dplyr::rename(
      !!purrr::set_names(tidyselect::contains(pred_var), 
                         c(glue::glue("m_pred{i}", 
                                      i=seq_along(pred_list))))) %>%
    dplyr::select(
      tidyselect::all_of(c(time_var, id_var,treat_time_var,outcome_var,
                           c(glue::glue("m_pred{i}", 
                                        i=seq_along(pred_list))))))
  if(!is.null(counterfac_var)){
    tib_out=tib_out %>% 
      dplyr::inner_join(pred_list[[1]] %>% 
                          dplyr::select(tidyselect::all_of(c(counterfac_var,
                                                           time_var,id_var,
                                                           "cf_point.effect",
                                                           "cf_pct.effect"))),
                        by=c(time_var, id_var))
  }
  
  return(tib_out)
  
}


weight_solver<-function(post_treat_combined_df,intercept_allowed,constrained,
                        method_df_names,counterfac_var ){
  x=post_treat_combined_df %>% 
    dplyr::select(tidyselect::contains("m_pred")) %>%
    as.matrix()
  
  if(intercept_allowed) x=cbind(1,x)
  
  r_inv <- tryCatch({
    solve(chol(t(x) %*% x))
  }, error=function(e){
    "Failure"
  })
  
  c=cbind(rep(1, length(method_df_names)), diag(length(method_df_names)))
  
  if(intercept_allowed) c=t(cbind(0, rbind(1, diag(length(method_df_names)))))
  
  b <- c(1, rep(0, length(method_df_names)))
  d <- t(post_treat_combined_df %>% pull(!!as.name(counterfac_var))) %*% x
  nn2 <- sqrt(norm(d, "2"))
  weight_vec=tryCatch( 
    {
      constr_weights_sol=solve.QP(Dmat = r_inv * nn2, factorized = TRUE, 
                                  dvec = d / (nn2^2), Amat = c, bvec = b, meq = 1)
      
      
      if(!constrained) constr_weights_sol$unconstrained.solution
      else constr_weights_sol$solution
      
      },
    error=function(e){
      if(intercept_allowed){
        error_out=c(0,rep(1/length(method_df_names), length(method_df_names)))
      } else  rep(1/length(method_df_names), length(method_df_names))
    })
  
 
  
  
  return(weight_vec)
}



ensemble_placebo_weights <- function(method_names, combined_methods_df, indiv_weights,
                                     constrained, intercept_allowed, time_var, 
                                     id_var, treat_time_var, counterfac_var) {
  
  # Estimates the weights for an ensemble using linear regression on a placebo set
  # This method is, as of now, highly specific -- intended for a placebo data set
  # with fake treatment (but zero treatment effect). Given this structure, we have fit several SCM methods
  # on this placebo set to estimate the TE (hopefully close to zero)
  # This method ignores any pre-treat information and outputs weights on the post-treatment point predictions
  # from each method that yields the closest value to the truth. Thus, we must know the truth to compute these.
  
  
  # Args
  # method1_estimated_df: long-form dataframe of the estimated and counterfactual point predictions for a given method
  # method2_estimated_df, method3_estimated_df distinct methods of the same form.
  # id_var: column name of numeric, unique ID representing the entry (unit)
  # time_var:column name of numeric period number indicating the time period, in increasing order (eg 0 is the first time, 120 is the last)
  # treat_time_var: string variable for the column indicating when that unit was treated
  # pred_var: string name of the point estimate column in the methodX dataframe
  # counterfac_var: name of var for the counterfactual value for the time series, if available (otherwise, NULL)
  # pos_coef_constr: Boolean flag for whether the weights should be non-negative
  # sum_to_one: Boolean flag for whether the coefficients (without intercept) sum to 1
  
  # Output
  # weights, as a 4x1 (num methods X 1) vector, from an unconstrained linear reg
  #first, some cleaning to keep only the relevant columns
    
  # Find the weights of the three methods that best fit the data in the post period
  # For now, we do this with simple linear regression and no constraints
  post_treat_combined_df <- combined_methods_df %>% filter(!!as.name(time_var) >= !!as.name(treat_time_var))
  
  #estimate one set of weights for each method, not unit specific
  if(indiv_weights){
    list_to_ensemble=post_treat_combined_df %>% split(.[[id_var]])
    tib_result=furrr::future_map(.x=list_to_ensemble,
                      .f=~weight_solver(.x,intercept_allowed,
                                        constrained,
                                        method_names,counterfac_var)) %>%
      dplyr::bind_rows() %>% 
      t() %>% 
      tibble::as_tibble(rownames = id_var, .name_repair="unique")  %>% 
      setNames(c("entry",paste0("V", seq_len(ncol(.)))))
    
    if(intercept_allowed){
      tib_result=tib_result %>% dplyr::rename(Intercept="V1")
    }
    
    return(tib_result %>%
      dplyr::rename(
        !!purrr::set_names(tidyselect::contains("V"), 
                           c(glue::glue("m_weight{i}", 
                                        i=seq_along(method_names))))))
    
    
   
  }else{
    df_result=weight_solver(post_treat_combined_df,intercept_allowed,constrained,
                  method_names,counterfac_var) %>% 
      as.list() %>%
      data.frame()
    
    tib_result=df_result %>% 
      tibble::tibble(.name_repair = "universal") %>% 
      setNames(c(paste0("V", seq_len(ncol(.)))))
      

      if(intercept_allowed){
        tib_result=tib_result %>% dplyr::rename(Intercept="V1")
      }
    return(tib_result %>%
             dplyr::rename(
               !!purrr::set_names(tidyselect::contains("V"), 
                                  c(glue::glue("m_weight{i}", 
                                               i=seq_along(method_names))))))
  }
}


placebo_estimation <- function(method_name, placebo_data,...){
  #define the estimator call, and apply. MC is the only exception to the rule
  if(method_name=="mc"){
    return(estimate_gsynth_series(placebo_data, estimator="mc",...))
  }
  if(method_name=="scdid_uncon"){
    return(estimate_scdid_series(placebo_data, constrained=F,...))
  }
  estimator_call=paste("estimate",method_name,"series", sep="_")
  return(do.call(estimator_call, list(placebo_data,...)))
}

#TODO(alexdkellogg): Potential improvement - focus ensemble weights on k periods
#    after treatment, rather than the full post-treat time period
estimate_ensemble <- function(method_names, true_data, pred_list,
                              constrained = T, intercept_allowed = T,  
                              indiv_weights=F, time_var = "period", 
                              id_var = "entry", outcome_var="response",
                              treat_time_var = "Treatment_Period",
                              pred_var = "point.pred", 
                              counterfac_var = "counter_factual"){
  
  #First step, create a placebo dataset from true_data via matching
  placebo_data=create_placebo_df(true_data)
  
  unit_mapping=placebo_data %>% 
    dplyr::filter(!is.na(Treatment_Period)) %>%
    dplyr::distinct(entry, Treatment_Unit)
  #Next, estimate the counterfactual series for each method
  #due to a bug in future library, must call these by name to tell future
  #this is a global object
  estimate_bart_series
  estimate_gfoo_series
  estimate_gsynth_series
  estimate_scdid_series
  estimate_scm_series
  
  estimates_list=furrr::future_map(.x=method_names, .f=placebo_estimation,
                        placebo_data = placebo_data)
  # estimates_list=lapply(method_names, placebo_estimation,
  #                       placebo_data = placebo_data)
  names(estimates_list)=method_names
  
  #Combine the list of estimates (one per method) into a single tibble
  combined_pred<- reduce_pred_list(estimates_list)
  
  #Estimate the weights on each of the methods
  method_weights<- 
    ensemble_placebo_weights(method_names, combined_methods_df=combined_pred, 
                             indiv_weights=indiv_weights,
                           constrained+constrained, 
                           intercept_allowed=intercept_allowed,
                           time_var = time_var, 
                           id_var = id_var, treat_time_var = treat_time_var,
                           counterfac_var = counterfac_var)
  

  #combine the existing estimates on true data into 1 tibble to extract preds
  
  combined_true_methods=reduce_pred_list(pred_list)
  if(indiv_weights){
    matched_weights=unit_mapping %>% 
      dplyr::inner_join(method_weights %>% dplyr::mutate_all((as.numeric)),
                        by=id_var) %>%
      dplyr::select(-tidyselect::all_of(id_var)) %>%
      dplyr::rename(entry="Treatment_Unit")
    
    combined_true_methods=combined_true_methods %>%
      dplyr::left_join(matched_weights, by=id_var)
    
    return(weighted_predictions(combined_true_methods, intercept_allowed))
    
  }else{
    combined_true_methods=combined_true_methods %>% 
      dplyr::bind_cols(method_weights)
    
    
    return(weighted_predictions(combined_true_methods, intercept_allowed))
  }
  
  
  
}


weighted_predictions <- function(data, intercept){
  method_preds=data %>% 
    dplyr::select(tidyselect::contains("m_pred")) %>%
    as.matrix()
  
  if(intercept) method_preds=cbind(1,method_preds)
  
  method_weights=data %>% 
    dplyr::select(tidyselect::matches("Intercept|m_weight")) %>%
    as.matrix()
  
  data=data %>% 
    dplyr::mutate("point.pred"=c(rowSums(method_preds * method_weights)),
                  "point.effect"=response-point.pred)
  
  return(data)
  
}

preprocess_ensemble_input <- function(estimated_list){
  if(inherits(estimated_list[[1]], "list")){
    estimated_list <- lapply(seq_len(min(lengths(estimated_list))), 
                         function(x) return(lapply(estimated_list, `[[`, x)))
  }

  return(estimated_list)
}



