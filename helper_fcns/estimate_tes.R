### This script implements various estimators for the treatment effect among
### compliers in the target population. They include:
###  - EIF-based estimators
###  - IPW-based estimators
###  - OM-based estimators

### Each of the above three categories have their own functions

estimate_eif_te <- function(data){
  
  ### We'll break this down into two pieces: 
  ### 1. Estimation of E[Y(1)|R=0, C(1)=1, C(0)=0], and
  ### 2. Estimation of E[Y(0)|R=0, C(1)=1, C(0)=0]
  
  ### Our final estimates will be formed by taking a contrast of these
  ### individual estimates. Note that we consider both normalized and 
  ### unnormalized estimators.
  
  ### In our data, each id already has an appropriate cross fit estimator for
  ### every nuisance function. To construct a cross-fit estimator of the whole
  ### functional, we take a weighted average of the estimators in each of 2
  ### groups. Because the group-specific estimators are themselves means, 
  ### this is equivalent to taking the mean of id-specific cross-fit estimates
  ### over our entire sample at once, which is how we proceed below.
  
  ### Setting A, C, and Y equal to zero for those in the target so that
  ### when we take the mean over terms A, C, Y, data from those with 
  ### R=0 contribute nothing.
  eif_est_data <- copy(data)[R == 0, `:=` (A=0, C=0, Y=0)]
  
  ### Start by estimating pieces of the denominator, which are the same across
  ### both parts of the contrast
  
  eif_denom_est_data <- copy(eif_est_data)
  eif_denom_est_data[, rho_resid_denom := (ps_1_cf_pred-ps_0_cf_pred)*((1-R-(1-rho_cf_pred)))]
  eif_denom_est_data[, p_1_resid_denom := (1-rho_cf_pred)*((A*R*(C-ps_1_cf_pred))/(pi_cf_pred*rho_cf_pred))]
  eif_denom_est_data[, p_0_resid_denom := (1-rho_cf_pred)*(((1-A)*R*(C-ps_0_cf_pred)))/((1-pi_cf_pred)*rho_cf_pred)]
  plug_in_outcome_denom <- data[, mean((ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred))]
  
  eif_denom_est_data[, total_resid := (1-rho_cf_pred)*(ps_1_cf_pred-ps_0_cf_pred)-plug_in_outcome_denom]
  
  denom_uncentered_eif <- eif_denom_est_data[, mean(rho_resid_denom)]+
    eif_denom_est_data[, mean(p_1_resid_denom-p_0_resid_denom)]+
    plug_in_outcome_denom
  
  normalized_denom_uncentered_eif <- eif_denom_est_data[, sum(p_1_resid_denom)/sum((A*R)/(pi_cf_pred*rho_cf_pred))]-
    eif_denom_est_data[, sum(p_0_resid_denom)/sum(((1-A)*R)/((1-pi_cf_pred)*rho_cf_pred))]+
    eif_denom_est_data[, mean(rho_resid_denom)]+
    plug_in_outcome_denom
  
  # Estimation of E[Y(1)|R=0, C(1)=1, C(0)=0] --------------------------------
  
  ### Constructing components of the EIF
  
  eif_1_est_data <- copy(eif_est_data)
  eif_1_est_data[, rho_resid := (ps_1_cf_pred-ps_0_cf_pred)*(mu_11_cf_pred)*((1-R)-(1-rho_cf_pred))]
  eif_1_est_data[, mu_resid := (1-rho_cf_pred)*(ps_1_cf_pred-ps_0_cf_pred)*((A*C*R*(Y-mu_11_cf_pred))/(pi_cf_pred*rho_cf_pred*ps_1_cf_pred))]
  eif_1_est_data[, p_1_resid := (1-rho_cf_pred)*(mu_11_cf_pred)*((A*R*(C-ps_1_cf_pred))/(pi_cf_pred*rho_cf_pred))]
  eif_1_est_data[, p_0_resid := (1-rho_cf_pred)*(mu_11_cf_pred)*(((1-A)*R*(C-ps_0_cf_pred)))/((1-pi_cf_pred)*rho_cf_pred)]
  
  plug_in_outcome_num_1 <- data[, mean((ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred)*(mu_11_cf_pred))]
  
  eif_1_est_data[, total_resid := (ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred)*(mu_11_cf_pred)-plug_in_outcome_num_1]
  
  num_1_uncentered_eif <- eif_1_est_data[, mean(mu_resid)]+
    eif_1_est_data[, mean(p_1_resid-p_0_resid)]+
    eif_1_est_data[, mean(rho_resid)]+
    plug_in_outcome_num_1
  
  normalized_num_1_uncentered_eif <- eif_1_est_data[, sum(mu_resid)/sum(A*C*R/(pi_cf_pred*rho_cf_pred*ps_1_cf_pred))]+
    eif_1_est_data[, sum(p_1_resid)/sum((A*R)/(pi_cf_pred*rho_cf_pred))]-
    eif_1_est_data[, sum(p_0_resid)/sum(((1-A)*R)/((1-pi_cf_pred)*rho_cf_pred))]+
    eif_1_est_data[, mean(rho_resid)]+
    plug_in_outcome_num_1
  
  y_1_eif_estimator <- num_1_uncentered_eif/denom_uncentered_eif
  
  normalized_y_1_eif_estimator <- normalized_num_1_uncentered_eif/normalized_denom_uncentered_eif
  
  # Estimation of E[Y(0)|R=0, C(1)=1, C(0)=0] --------------------------------
  
  eif_0_est_data <- copy(eif_est_data)
  eif_0_est_data[, rho_resid := (ps_1_cf_pred-ps_0_cf_pred)*(mu_00_cf_pred)*((1-R)-(1-rho_cf_pred))]
  eif_0_est_data[, mu_resid := (1-rho_cf_pred)*(ps_1_cf_pred-ps_0_cf_pred)*(((1-A)*(1-C)*R*(Y-mu_00_cf_pred))/((1-pi_cf_pred)*(1-ps_0_cf_pred)*rho_cf_pred))]
  eif_0_est_data[, p_1_resid := (1-rho_cf_pred)*(mu_00_cf_pred)*((A*R*(C-ps_1_cf_pred))/(pi_cf_pred*rho_cf_pred))]
  eif_0_est_data[, p_0_resid := (1-rho_cf_pred)*(mu_00_cf_pred)*(((1-A)*R*(C-ps_0_cf_pred))/((1-pi_cf_pred)*(rho_cf_pred)))]
  
  plug_in_outcome_num_0 <- data[, mean((ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred)*(mu_00_cf_pred))]
  
  eif_0_est_data[, total_resid := (ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred)*(mu_00_cf_pred)-plug_in_outcome_num_0]
  
  num_0_uncentered_eif <- eif_0_est_data[, mean(mu_resid)]+
    eif_0_est_data[, mean(p_1_resid-p_0_resid)]+
    eif_0_est_data[, mean(rho_resid)]+
    plug_in_outcome_num_0
  
  normalized_num_0_uncentered_eif <- eif_0_est_data[, sum(mu_resid)/sum(((1-A)*(1-C)*R)/((1-pi_cf_pred)*(1-ps_0_cf_pred)*rho_cf_pred))]+
    eif_0_est_data[, sum(p_1_resid)/sum((A*R)/(pi_cf_pred*rho_cf_pred))]-
    eif_0_est_data[, sum(p_0_resid)/sum(((1-A)*R)/((1-pi_cf_pred)*rho_cf_pred))]+
    eif_0_est_data[, mean(rho_resid)]+
    plug_in_outcome_num_0
  
  y_0_eif_estimator <- num_0_uncentered_eif/denom_uncentered_eif
  
  normalized_y_0_eif_estimator <- normalized_num_0_uncentered_eif/normalized_denom_uncentered_eif
  
  # Estimation of E[Y(1)-Y(0)|R=0, C(1)=1, C(0)=0] --------------------------------
  
  te_eif_estimator <- y_1_eif_estimator-y_0_eif_estimator
  
  normalized_te_eif_estimator <- normalized_y_1_eif_estimator-normalized_y_0_eif_estimator
  
  ### Last piece: variance estimation. We're gonna see how the asymptotic approximation
  ### works for obtaining confidence intervals. This requires estimating the variance of
  ### the EIF-based plug-in estimator, which we accomplish by computing the empirical
  ### variance of the uncentered influence function, derived elsewhere. Note that this
  ### empirical variance requires an estimate of the treatment effect itself, so
  ### we'll consider both normalized and un-normalized versions of the variance estimator 
  ### according to whether we use normalized or un-normalized versions of the treatment
  ### effect estimator, i.e., `te_eif_estimator` or `normalized_te_eif_estimator`.
  
  var_eif_est_data <- list(eif_0_est_data[, .(id, centered_eif_0 = mu_resid+(p_1_resid-p_0_resid)+rho_resid+total_resid)], 
                           eif_1_est_data[, .(id, centered_eif_1 = mu_resid+(p_1_resid-p_0_resid)+rho_resid+total_resid)], 
                           eif_denom_est_data[, .(id, centered_eif_denom = rho_resid_denom+(p_1_resid_denom-p_0_resid_denom)+total_resid)]) %>%
    Reduce(function(x, y){merge(x, y, by = "id")}, x = .) %>%
    .[, centered_eif_num := centered_eif_1-centered_eif_0] %>%
    .[, centered_eif := (centered_eif_num-te_eif_estimator*centered_eif_denom)/(denom_uncentered_eif)] %>%
    .[, normalized_centered_eif := (centered_eif_num-normalized_te_eif_estimator*centered_eif_denom)/(normalized_denom_uncentered_eif)]
  
  var_eif_estimator <- var_eif_est_data[, var(centered_eif)]/nrow(data)
  
  normalized_var_eif_estimator <- var_eif_est_data[, var(normalized_centered_eif)]/nrow(data)
  
  return(data.table(te_eif_estimator = te_eif_estimator, 
                    var_eif_estimator = var_eif_estimator,
                    normalized_te_eif_estimator = normalized_te_eif_estimator, 
                    normalized_var_eif_estimator = normalized_var_eif_estimator))
  
}

estimate_om_te <- function(data){
  
  ### Cross-fit version
  
  num_1_plug_in_cross_fit <- data[, mean((ps_1_cf_pred-ps_0_cf_pred)*(1-R)*(mu_11_cf_pred))]
  num_0_plug_in_cross_fit <- data[, mean((ps_1_cf_pred-ps_0_cf_pred)*(1-R)*(mu_00_cf_pred))]
  denom_plug_in_cross_fit <- data[, mean((ps_1_cf_pred-ps_0_cf_pred)*(1-R))]
  
  te_plug_in_cross_fit <- (num_1_plug_in_cross_fit-num_0_plug_in_cross_fit)/(denom_plug_in_cross_fit)
  
  ### Non cross-fit version
  
  num_1_plug_in <- data[, mean((ps_1_full_data_pred-ps_0_full_data_pred)*(1-R)*(mu_11_full_data_pred))]
  num_0_plug_in <- data[, mean((ps_1_full_data_pred-ps_0_full_data_pred)*(1-R)*(mu_00_full_data_pred))]
  denom_plug_in <- data[, mean((ps_1_full_data_pred-ps_0_full_data_pred)*(1-R))]
  
  te_plug_in <- (num_1_plug_in-num_0_plug_in)/(denom_plug_in)
  
  return(data.table(te_om_cf_estimator = te_plug_in_cross_fit, 
                    te_om_full_data_estimator = te_plug_in))
  
}

estimate_plug_in_te <- function(data){
  
  ### Cross-fit version
  
  num_1_plug_in_cross_fit <- data[, mean((ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred)*(mu_11_cf_pred))]
  num_0_plug_in_cross_fit <- data[, mean((ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred)*(mu_00_cf_pred))]
  denom_plug_in_cross_fit <- data[, mean((ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred))]
  
  te_plug_in_cross_fit <- (num_1_plug_in_cross_fit-num_0_plug_in_cross_fit)/(denom_plug_in_cross_fit)
  
  ### Non cross-fit version
  
  num_1_plug_in <- data[, mean((ps_1_full_data_pred-ps_0_full_data_pred)*(1-rho_full_data_pred)*(mu_11_full_data_pred))]
  num_0_plug_in <- data[, mean((ps_1_full_data_pred-ps_0_full_data_pred)*(1-rho_full_data_pred)*(mu_00_full_data_pred))]
  denom_plug_in <- data[, mean((ps_1_full_data_pred-ps_0_full_data_pred)*(1-rho_full_data_pred))]
  
  te_plug_in <- (num_1_plug_in-num_0_plug_in)/(denom_plug_in)
  
  return(data.table(te_plug_in_cf_estimator = te_plug_in_cross_fit, 
                    te_plug_in_full_data_estimator = te_plug_in))
  
}

estimate_ipw_te <- function(data){
  
  ipw_est_data <- copy(data)[R==0, `:=` (A=0, C=0, Y=0)]
  
  ### Cross-fit version
  
  ipw_est_data[, wtd_y_1_cf := ((ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred)*(A*C*R*Y)/(pi_cf_pred*ps_1_cf_pred*rho_cf_pred))]
  ipw_est_data[, wtd_y_0_cf := ((ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred)*((1-A)*(1-C)*R*Y)/((1-pi_cf_pred)*(1-ps_0_cf_pred)*(rho_cf_pred)))]
  denom_ipw_cf <- data[, mean((ps_1_cf_pred-ps_0_cf_pred)*(1-rho_cf_pred))]
  
  ## Normalized
  
  norm_ipw_est_cf_y_1 <- ipw_est_data[, sum(wtd_y_1_cf)/sum((A*C*R)/(pi_cf_pred*ps_1_cf_pred*rho_cf_pred))]
  norm_ipw_est_cf_y_0 <- ipw_est_data[, sum(wtd_y_0_cf)/sum(((1-A)*(1-C)*R)/((1-pi_cf_pred)*(1-ps_0_cf_pred)*(rho_cf_pred)))]
  
  norm_ipw_est_cf <- (norm_ipw_est_cf_y_1-norm_ipw_est_cf_y_0)/denom_ipw_cf
  
  ## Non-Normalized
  
  ipw_est_cf <- (ipw_est_data[, mean(wtd_y_1_cf)]-ipw_est_data[, mean(wtd_y_0_cf)])/(denom_ipw_cf)
  
  ### Non cross-fit version
  
  ipw_est_data[, wtd_y_1 := ((ps_1_full_data_pred-ps_0_full_data_pred)*(1-rho_full_data_pred)*(A*C*R*Y)/(pi_full_data_pred*ps_1_full_data_pred*rho_full_data_pred))]
  ipw_est_data[, wtd_y_0 := ((ps_1_full_data_pred-ps_0_full_data_pred)*(1-rho_full_data_pred)*((1-A)*(1-C)*R*Y)/((1-pi_full_data_pred)*(1-ps_0_full_data_pred)*(rho_full_data_pred)))]
  denom_ipw <- data[, mean((ps_1_full_data_pred-ps_0_full_data_pred)*(1-rho_full_data_pred))]
  
  ## Normalized
  
  norm_ipw_est_y_1 <- ipw_est_data[, sum(wtd_y_1)/sum((A*C*R)/(pi_full_data_pred*ps_1_full_data_pred*rho_full_data_pred))]
  norm_ipw_est_y_0 <- ipw_est_data[, sum(wtd_y_0)/sum(((1-A)*(1-C)*R)/((1-pi_full_data_pred)*(1-ps_0_full_data_pred)*(rho_full_data_pred)))]
  
  norm_ipw_est <- (norm_ipw_est_y_1-norm_ipw_est_y_0)/denom_ipw
  
  ## Non-Normalized
  
  ipw_est <- (ipw_est_data[, mean(wtd_y_1)]-ipw_est_data[, mean(wtd_y_0)])/(denom_ipw)
  
  ### Bringing them all together:
  
  ipw_est_table <- data.table(te_ipw_cf_estimator = ipw_est_cf, 
                              te_ipw_full_data_estimator = ipw_est, 
                              normalized_te_ipw_cf_estimator = norm_ipw_est_cf, 
                              normalized_te_ipw_full_data_estimator = norm_ipw_est)
  
  return(ipw_est_table)
  
}
