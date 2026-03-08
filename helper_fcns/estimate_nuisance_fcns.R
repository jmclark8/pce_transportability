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
    
    fold_model <- glm(model_string, family = binomial(link = "logit"), 
                      data = data[cf_fold != fold])
    
    fold_data <- data[cf_fold == fold, ..pred_data_cols]
    
    p_0_preds <- copy(fold_data)[, A := 0][, ps_0_cf_pred := predict.glm(fold_model, newdata = .SD, 
                                                                         type = "response")]
    p_1_preds <- copy(fold_data)[, A := 1][, ps_1_cf_pred := predict.glm(fold_model,
                                                                         newdata = .SD, 
                                                                         type = "response")]
    fold_preds <- merge(p_0_preds[, .(id, ps_0_cf_pred)], 
                        p_1_preds[, .(id, ps_1_cf_pred)], 
                        by = "id")
    
    fold_preds_list[[fold_preds_name]] <- fold_preds
    
  }
  
  fold_preds_data <- rbindlist(fold_preds_list) ### Putting all the cross-fitting predictions
  ### in one dataset
  
  ### Making predictions using a model fit on the full data:
  
  full_data_model <- glm(model_string, family = binomial(link = "logit"), 
                         data = data)
  p_0_preds_full_data <- copy(data)[, ..pred_data_cols] %>%
    .[, A:= 0] %>%
    .[, ps_0_full_data_pred := predict.glm(full_data_model, 
                                           newdata = .SD, type = "response")] %>%
    .[, .(id, ps_0_full_data_pred)]
  
  p_1_preds_full_data <- copy(data)[, ..pred_data_cols] %>%
    .[, A:= 1] %>%
    .[, ps_1_full_data_pred := predict.glm(full_data_model, 
                                           newdata = .SD, type = "response")] %>%
    .[, .(id, ps_1_full_data_pred)]
  
  full_data_preds <- merge(p_0_preds_full_data, p_1_preds_full_data, 
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
    
    fold_model <- lm(model_string, data = data[cf_fold != fold])
    
    fold_data <- data[cf_fold == fold, ..pred_data_cols]
    
    ### For the purposes of this simulation, we only need predictions when C=1
    ### and A=1 and when C=0 and A=0
    
    mu_11_preds <- copy(fold_data)[, `:=` (A = 1, C = 1)] %>%
      .[, mu_11_cf_pred := predict(fold_model, newdata = .SD)] %>%
      .[, .(id, mu_11_cf_pred)]
    
    mu_00_preds <- copy(fold_data)[, `:=` (A = 0, C = 0)] %>%
      .[, mu_00_cf_pred := predict(fold_model, newdata = .SD)]
    
    fold_preds <- merge(mu_11_preds[, .(id, mu_11_cf_pred)], 
                        mu_00_preds[, .(id, mu_00_cf_pred)], 
                        by = "id")
    
    fold_preds_list[[fold_preds_name]] <- fold_preds
    
  }
  
  fold_preds_data <- rbindlist(fold_preds_list)
  
  ### Making predictions using a model fit on the full data:
  
  full_data_model <- lm(model_string, data = data)
  
  mu_11_preds_full_data <- copy(data)[, ..pred_data_cols] %>%
    .[, `:=` (A = 1, C = 1)] %>%
    .[, mu_11_full_data_pred := predict(full_data_model, newdata = .SD)] %>%
    .[, .(id, mu_11_full_data_pred)]
  
  mu_00_preds_full_data <- copy(data)[, ..pred_data_cols] %>%
    .[, `:=` (A = 0, C = 0)] %>%
    .[, mu_00_full_data_pred := predict(full_data_model, newdata = .SD)] %>%
    .[, .(id, mu_00_full_data_pred)]
  
  full_data_preds <- merge(mu_11_preds_full_data, mu_00_preds_full_data, 
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
    
    fold_model <- glm(model_string, family = binomial(link = "logit"), 
                      data = data[cf_fold != fold])
    
    fold_data <- data[cf_fold == fold, ..pred_data_cols]
    
    rho_preds <- copy(fold_data) %>%
      .[, rho_cf_pred := predict.glm(fold_model, newdata = .SD, 
                                     type = "response")]
    
    fold_preds <- rho_preds[, .(id, rho_cf_pred)]
    
    fold_preds_list[[fold_preds_name]] <- fold_preds
    
  }
  
  fold_preds_data <- rbindlist(fold_preds_list)
  
  ### Making predictions using a model fit on the full data:
  
  full_data_model <- glm(model_string, family = binomial(link = "logit"), 
                         data = data)
  
  full_data_preds <- copy(data)[, ..pred_data_cols] %>%
    .[, rho_full_data_pred := predict.glm(full_data_model, newdata = .SD, 
                                          type = "response")] %>%
    .[, .(id, rho_full_data_pred)]
  
  all_preds <- merge(fold_preds_data, full_data_preds, by = "id")
  
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
    
    fold_model <- glm(model_string, family = binomial(link = "logit"), 
                      data = data[cf_fold != fold])
    
    fold_data <- data[cf_fold == fold, ..pred_data_cols]
    
    pi_preds <- copy(fold_data) %>%
      .[, pi_cf_pred := predict.glm(fold_model, newdata = .SD, 
                                    type = "response")]
    
    fold_preds <- pi_preds[, .(id, pi_cf_pred)]
    
    fold_preds_list[[fold_preds_name]] <- fold_preds
    
  }
  
  fold_preds_data <- rbindlist(fold_preds_list)
  
  ### Making predictions using a model fit on the full data:
  
  full_data_model <- glm(model_string, family = binomial(link = "logit"), 
                         data = data)
  
  full_data_preds <- copy(data)[, ..pred_data_cols] %>%
    .[, pi_full_data_pred := predict.glm(full_data_model, newdata = .SD, 
                                         type = "response")] %>%
    .[, .(id, pi_full_data_pred)]
  
  all_preds <- merge(fold_preds_data, full_data_preds, by = "id")
  
  return(all_preds)
  
}





