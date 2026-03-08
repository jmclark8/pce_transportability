### This function simulates the covariate distributions and assigns rows to
### the trial or target populations based on the assumed data generating process.


sim_pop_samples <- function(total_n, ratio){
  
  ### Generating a sample from the marginal distribution of 
  ### covariates.
  
  marginal_covar_sample_cts <- data.table(x_1 = runif(n = total_n, min = -2, max = 2),
                                          x_2 = runif(n = total_n, min = -2, max = 2),
                                          x_3 = runif(n = total_n, min = -2, max = 2),
                                          x_4 = runif(n = total_n, min = -2, max = 2))
  
  marginal_covar_sample_binary <- data.table(x_5 = rbinom(n = total_n, 
                                                          size = 1, 
                                                          prob = 0.55))
  
  marginal_covar_sample <- cbind(marginal_covar_sample_cts,
                                 marginal_covar_sample_binary) %>%
    ### Adding the transformed covariates we'll use in the model
    .[, c("c_1_tilde", "c_2_tilde", "c_3_tilde", "c_4_tilde") := lapply(.SD, 
                                                                        function(x){(x^2-1)/sqrt(2)}), 
      .SDcols = c("x_1", "x_2", "x_3", "x_4")] %>%
    ### Adding the transformed covariates we'll use in the wrong model
    .[, c("c_1_wrong", "c_2_wrong", "c_3_wrong", "c_4_wrong") := lapply(.SD, 
                                                                        function(x){x-0.25}), 
      .SDcols = c("x_1", "x_2", "x_3", "x_4")] %>%
    .[, id := 1:total_n]
  
  ### Great! Each row represents a draw from the marginal distribution of
  ### covariates. Now we can assign each such draw a probability of being in the
  ### trial population, using the appropriate intercept so our sample
  ### sizes are good to go. Note that the `intercept` dataset needs to be
  ### loaded into the environment for this to work.
  
  int_ratio <- ratio
  
  balanced_intercept <- intercepts[ratio == int_ratio, unique(int_value)]
  
  marginal_covar_sample[, rho := expit(balanced_intercept + c_1_tilde + c_2_tilde + c_3_tilde), 
                        by = "id"]
  
  ### Now assigning rows to trial participation or not:
  
  marginal_covar_sample[, R := rbinom(n = 1, size = 1, prob = rho), 
                        by = "id"]
  
  ### Awesome! So, to summarize, we can treat this data frame as a sample
  ### from the marginal distribution of x of size n. The relative sample sizes
  ### of target and trial participants are driven by the intercept in the logit
  ### model that associates covariate patterns with each source. In essence, 
  ### we've guaranteed that E[rho(X)] is approximately equal to the relative
  ### sample sizes of trial vs. target in our sample.
  
  return(marginal_covar_sample)
  
  
}








