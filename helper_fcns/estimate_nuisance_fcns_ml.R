### This script constructs estimates of various nuisance functions, including:

### 1. The principal score
### 2. The outcome model
### 3. The participation probability
### 4. The treatment probability

### It does so via four separate functions for each such model. All of those
### functions take as arguments some version of the following: whether the
### model should be estimated via cross-fitting, whether the model
### is correctly specified and the simulated data used to fit the model. 

### Each function returns a `data.table` with sample ids and corresponding
### nuisance function estimates which, in the main simulation function, are
### ultimately merged back with the larger dataset.

### Again, each model returns two sets of estimates: one fit separately for
### each fold, another fit on the entire sample

### Update: 2025-08-03: We're adding ML-type estimators to perform the 
### nuisance function estimation. This was in response to a reviewer's request.

estimate_ps <- function(data, ps, n_folds){
  
  if (ps == TRUE){
    model_string <- C ~ A*(c_1_tilde + c_2_tilde + c_3_tilde + c_4_tilde)
    preds <- c("A", "c_1_tilde", "c_2_tilde", "c_3_tilde", "c_4_tilde")
  } else if (ps == FALSE){
    model_string <- C ~ A*(c_1_wrong + c_2_wrong + c_3_wrong + c_4_wrong)
    preds <- c("A", "c_1_wrong", "c_2_wrong", "c_3_wrong", "c_4_wrong")
  }
  
  ### Fitting models to each fold and then predicting:
  
  pred_data_cols <- c("id", "cf_fold", setdiff(preds, "A"))
  
  fold_preds_list <- vector(mode = "list", length = n_folds)
  
  names(fold_preds_list) <- paste0("fold_", 1:n_folds, "_preds")
  
  for (fold in 1:n_folds){
    
    fold_preds_name <- paste0("fold_", fold, "_preds")
    
    ### Doing machine learning here instead of logistic regression.
    
    fold_model_data <- data[cf_fold != fold] %>%
      ### For SuperLearner, removing missing data
      .[!is.na(C)]
    fold_model_outcomes <- as.numeric(fold_model_data$C)
    fold_model_preds <- fold_model_data[, ..preds]
    
    fold_model <- SuperLearner(Y = fold_model_outcomes, 
                               X = fold_model_preds, 
                               family = binomial(), 
                               SL.library = c("SL.glmnet", 
                                              "SL.ranger"))
    
    ### Note: For SuperLearner, I think we can only pass in `newdata` 
    ### a data frame with the same columns as appeared when fitting the model, i.e., 
    ### we can't include `cf_fold` or `id` in `newdata`. I remove those columns
    ### before calling predict.SuperLearner, but ultimately we'll want those
    ### predicted values re-associated with the IDs. Hence saving the IDs here
    ### and sorting by ID below.
    
    fold_ids <- data[cf_fold==fold, sort(id)]
    
    p_0_preds <- copy(data)[cf_fold==fold, ..pred_data_cols] %>%
      .[order(id)] %>%
      .[, `:=` (id = NULL, cf_fold = NULL)] %>%
      .[, A := 0] %>%
      .[, ps_0_cf_pred := predict.SuperLearner(fold_model, newdata = .SD)$pred]
    
    p_1_preds <- copy(data)[cf_fold==fold, ..pred_data_cols] %>%
      .[order(id)] %>%
      .[, `:=` (id = NULL, cf_fold = NULL)] %>%
      .[, A := 1] %>%
      .[, ps_1_cf_pred := predict.SuperLearner(fold_model, newdata = .SD)$pred]
    
    fold_preds <- merge(p_0_preds[, .(id = fold_ids, ps_0_cf_pred)], 
                        p_1_preds[, .(id = fold_ids, ps_1_cf_pred)], 
                        by = "id")
    
    fold_preds_list[[fold_preds_name]] <- fold_preds
    
  }
  
  fold_preds_data <- rbindlist(fold_preds_list) ### Putting all the cross-fitting predictions
                                                ### in one dataset
  
  ### Making predictions using a model fit on the full data:
  
  full_model_data <- data[!is.na(C)]
  
  full_data_outcomes <- as.numeric(full_model_data$C)
  full_data_preds <- full_model_data[, ..preds]
  
  all_ids <- data[, sort(id)]
  
  full_data_model <- SuperLearner(Y=full_data_outcomes, 
                                  X=full_data_preds, 
                                  family = binomial(), 
                                  SL.library=c("SL.glmnet", 
                                               "SL.ranger"))
  
  p_0_preds_full_data <- copy(data)[, ..pred_data_cols] %>%
    .[order(id)] %>%
    .[, `:=` (id = NULL, cf_fold = NULL)] %>%
    .[, A:= 0] %>%
    .[, ps_0_full_data_pred := predict.SuperLearner(full_data_model, newdata = .SD)$pred]
  
  p_1_preds_full_data <- copy(data)[, ..pred_data_cols] %>%
    .[order(id)] %>%
    .[, `:=` (id = NULL, cf_fold = NULL)] %>%
    .[, A:= 1] %>%
    .[, ps_1_full_data_pred := predict.SuperLearner(full_data_model, newdata = .SD)$pred]
  
  full_data_preds <- merge(p_0_preds_full_data[, .(id = all_ids, ps_0_full_data_pred)], 
                           p_1_preds_full_data[, .(id = all_ids, ps_1_full_data_pred)], 
                           by = "id")
  
  all_preds <- merge(fold_preds_data, full_data_preds, by = "id")
  
  return(all_preds)
  

}

estimate_om <- function(data, om, n_folds){
  
  if (om == TRUE){
    c_preds <- c("c_1_tilde", "c_2_tilde", "c_3_tilde", "c_4_tilde", "x_5")
  } else if (om == FALSE){
    c_preds <- c("c_1_wrong", "c_2_wrong", "c_3_wrong", "c_4_wrong", "x_5")
  }
  
  ### Making the model string
  
  model_string <- c(c_preds, paste0(c_preds, ":A"), paste0(c_preds, ":C")) %>%
    paste0(collapse = " + ") %>%
    paste("Y", ., sep = " ~ ") %>%
    ### Getting rid of the intercept term (already implicit in the covars)
    paste(., "-1") %>%
    as.formula()
  
  pred_data_cols <- pred_data_cols <- c("id", "cf_fold", c_preds)
  
  fold_preds_list <- vector(mode = "list", length = n_folds)
  
  names(fold_preds_list) <- paste0("fold_", 1:n_folds, "_preds")
  
  for (fold in 1:n_folds){
    
    fold_preds_name <- paste0("fold_", fold, "_preds")
    
    ### Doing machine learning here instead of logistic regression.
    
    preds <- c("A", "C", c_preds)
    
    fold_model_data <- data[cf_fold != fold] %>%
      ### For SuperLearner, removing missing data
      .[!is.na(Y)]
    fold_model_outcomes <- as.numeric(fold_model_data$Y)
    fold_model_preds <- fold_model_data[, ..preds]
    
    fold_model <- SuperLearner(Y=fold_model_outcomes, 
                               X=fold_model_preds, 
                               family=gaussian, 
                               SL.library=c("SL.glmnet", 
                                            "SL.ranger"))
    
    ### Note: For SuperLearner, I think we can only pass in `newdata` 
    ### a data frame with the same columns as appeared when fitting the model, i.e., 
    ### we can't include `cf_fold` or `id` in `newdata`. I remove those columns
    ### before calling predict.SuperLearner, but ultimately we'll want those
    ### predicted values re-associated with the IDs. Hence saving the IDs here
    ### and sorting by ID below.
    
    ### For the purposes of this simulation, we only need predictions when C=1
    ### and A=1 and when C=0 and A=0
    
    fold_ids <- data[cf_fold==fold, sort(id)]
    
    mu_11_preds <- copy(data)[cf_fold==fold, ..pred_data_cols] %>%
      .[order(id)] %>%
      .[, `:=` (id = NULL, cf_fold = NULL)] %>%
      .[, `:=` (A = 1, C = 1)] %>%
      .[, mu_11_cf_pred := predict.SuperLearner(fold_model, newdata = .SD)$pred]
    
    mu_00_preds <- copy(data)[cf_fold==fold, ..pred_data_cols] %>%
      .[order(id)] %>%
      .[, `:=` (id = NULL, cf_fold = NULL)] %>%
      .[, `:=` (A = 0, C = 0)] %>%
      .[, mu_00_cf_pred := predict.SuperLearner(fold_model, newdata = .SD)$pred]
    
    fold_preds <- merge(mu_11_preds[, .(id = fold_ids, mu_11_cf_pred)], 
                        mu_00_preds[, .(id = fold_ids, mu_00_cf_pred)], 
                        by = "id")
    
    fold_preds_list[[fold_preds_name]] <- fold_preds
    
  }
  
  fold_preds_data <- rbindlist(fold_preds_list)
  
  ### Making predictions using a model fit on the full data:
  
  preds <- c("A", "C", c_preds)
  
  full_model_data <- data[!is.na(Y)]
  
  full_data_outcomes <- as.numeric(full_model_data$Y)
  full_data_preds <- full_model_data[, ..preds]
  
  all_ids <- data[, sort(id)]
  
  full_data_model <- SuperLearner(Y=full_data_outcomes,
                                  X=full_data_preds,
                                  family=gaussian,
                                  SL.library=c("SL.glmnet", 
                                               "SL.ranger"))
  
  mu_11_preds_full_data <- copy(data)[, ..pred_data_cols] %>%
    .[order(id)] %>%
    .[, `:=` (id = NULL, cf_fold = NULL)] %>%
    .[, `:=` (A = 1, C = 1)] %>%
    .[, mu_11_full_data_pred := predict.SuperLearner(full_data_model, newdata = .SD)$pred]
  
  mu_00_preds_full_data <- copy(data)[, ..pred_data_cols] %>%
    .[order(id)] %>%
    .[, `:=` (id = NULL, cf_fold = NULL)] %>%
    .[, `:=` (A = 0, C = 0)] %>%
    .[, mu_00_full_data_pred := predict.SuperLearner(full_data_model, newdata = .SD)$pred]
  
  full_data_preds <- merge(mu_11_preds_full_data[, .(id = all_ids, 
                                                     mu_11_full_data_pred)], 
                           mu_00_preds_full_data[, .(id = all_ids, 
                                                     mu_00_full_data_pred)], 
                           by = "id")
  
  all_preds <- merge(fold_preds_data, full_data_preds, by = "id")
  
  return(all_preds)
}

estimate_pp <- function(data, pp, n_folds){
  
  if (pp == TRUE){
    preds <- c("c_1_tilde", "c_2_tilde", "c_3_tilde")
  } else if (pp == FALSE){
    preds <- c("c_1_wrong", "c_2_wrong", "c_3_wrong")
  }
  
  model_string <- paste("R", paste(preds, collapse = " + "), sep = " ~ ") %>%
    ### Getting rid of the intercept term (already implicit in the covars)
    # paste(., "-1") %>%
    ### coercing it to a formula
    as.formula()
  
  pred_data_cols <- c("id", "cf_fold", preds)
  
  fold_preds_list <- vector(mode = "list", length = n_folds)
  
  names(fold_preds_list) <- paste0("fold_", 1:n_folds, "_preds")
  
  for (fold in 1:n_folds){
    
    fold_preds_name <- paste0("fold_", fold, "_preds")
    
    ### Doing machine learning here instead of logistic regression.
    
    fold_model_data <- data[cf_fold != fold] %>%
      .[!is.na(R)]
    fold_model_outcomes <- as.numeric(fold_model_data$R)
    fold_model_preds <- fold_model_data[, ..preds]
    
    fold_model <- SuperLearner(Y = fold_model_outcomes,
                               X = fold_model_preds,
                               family = binomial(),
                               SL.library = c("SL.glmnet", 
                                              "SL.ranger"))
    
    fold_ids <- data[cf_fold == fold, sort(id)]
    
    rho_preds <- copy(data)[cf_fold==fold, ..pred_data_cols] %>%
      .[order(id)] %>%
      .[, `:=` (id = NULL, cf_fold = NULL)] %>%
      .[, rho_cf_pred := predict.SuperLearner(fold_model, newdata = .SD)$pred]
    
    fold_preds <- rho_preds[, .(id = fold_ids, rho_cf_pred)]
    
    fold_preds_list[[fold_preds_name]] <- fold_preds
    
  }
  
  fold_preds_data <- rbindlist(fold_preds_list)
  
  ### Making predictions using a model fit on the full data:
  
  full_model_data <- data[!is.na(R)]
  
  full_data_outcomes <- as.numeric(full_model_data$R)
  full_data_preds <- full_model_data[, ..preds]
  
  all_ids <- data[, sort(id)]
  
  full_data_model <- SuperLearner(Y=full_data_outcomes, 
                                  X=full_data_preds,
                                  family=binomial(),
                                  SL.library=c("SL.glmnet", 
                                               "SL.ranger"))
  
  full_data_preds <- copy(data)[, ..pred_data_cols] %>%
    .[order(id)] %>%
    .[, `:=` (id = NULL, cf_fold = NULL)] %>%
    .[, rho_full_data_pred := predict.SuperLearner(full_data_model, newdata = .SD)$pred]
  
  all_preds <- merge(fold_preds_data, 
                     full_data_preds[, .(id = all_ids, 
                                         rho_full_data_pred)], by = "id")
  
  return(all_preds)
}

estimate_tp <- function(data, tp, n_folds){
  
  if (tp == TRUE){
    preds <- c("c_1_tilde", "c_2_tilde", "c_3_tilde", "c_4_tilde")
  } else if (tp == FALSE){
    preds <- c("c_1_wrong", "c_2_wrong", "c_3_wrong", "c_4_wrong")
  }
  
  model_string <- paste("A", paste(preds, collapse = " + "), sep = " ~ ") %>%
    ### Getting rid of the intercept term (already implicit in the covars)
    paste(., "-1") %>%
    ### coercing it to a formula
    as.formula()
  
  pred_data_cols <- c("id", "cf_fold", preds)
  
  fold_preds_list <- vector(mode = "list", length = n_folds)
  
  for (fold in 1:n_folds){
    
    fold_preds_name <- paste0("fold_", fold, "_preds")
    
    ### Doing machine learning here instead of logistic regression.
    
    fold_model_data <- data[cf_fold != fold] %>%
      .[!is.na(A)]
    fold_model_outcomes <- as.numeric(fold_model_data$A)
    fold_model_preds <- fold_model_data[, ..preds]
    
    fold_model <- SuperLearner(Y = fold_model_outcomes, 
                               X = fold_model_preds, 
                               family = binomial(), 
                               SL.library = c("SL.glmnet", 
                                              "SL.ranger"))
    
    fold_ids <- data[cf_fold == fold, sort(id)]
    
    pi_preds <- copy(data)[cf_fold==fold, ..pred_data_cols] %>%
      .[order(id)] %>%
      .[, `:=` (id = NULL, cf_fold = NULL)] %>%
      .[, pi_cf_pred := predict.SuperLearner(fold_model, newdata = .SD)$pred]
    
    fold_preds <- pi_preds[, .(id = fold_ids, pi_cf_pred)]
    
    fold_preds_list[[fold_preds_name]] <- fold_preds
    
  }
  
  fold_preds_data <- rbindlist(fold_preds_list)
  
  ### Making predictions using a model fit on the full data:
  
  full_model_data <- data[!is.na(A)]
  
  full_data_outcomes <- as.numeric(full_model_data$A)
  full_data_preds <- full_model_data[, ..preds]
  
  all_ids <- data[, sort(id)]
  
  full_data_model <- SuperLearner(Y=full_data_outcomes, 
                                  X=full_data_preds, 
                                  family = binomial(), 
                                  SL.library=c("SL.glmnet", 
                                               "SL.ranger"))
  
  full_data_preds <- copy(data)[, ..pred_data_cols] %>%
    .[order(id)] %>%
    .[, `:=` (id = NULL, cf_fold = NULL)] %>%
    .[, pi_full_data_pred := predict.SuperLearner(full_data_model, newdata = .SD)$pred]
  
  all_preds <- merge(fold_preds_data, 
                     full_data_preds[, .(id = all_ids, 
                                         pi_full_data_pred)], by = "id")
  
  return(all_preds)
  
}





