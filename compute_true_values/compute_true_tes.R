library(here)
library(data.table)
library(magrittr)
library(MASS)
library(truncnorm)

expit <- function(x){
  exp(x)/(1+exp(x))
}

### This script computes the true value of our parameter of interest
### using monte carlo simulation.

get_truth <- function(n_reps, pop_ratio){
  
  n <- n_reps
  
  marginal_covar_sample_cts <- data.table(x_1 = runif(n = n, min = -2, max = 2),
                                          x_2 = runif(n = n, min = -2, max = 2),
                                          x_3 = runif(n = n, min = -2, max = 2),
                                          x_4 = runif(n = n, min = -2, max = 2))
  
  marginal_covar_sample_binary <- data.table(x_5 = rbinom(n = n, 
                                                          size = 1, 
                                                          prob = 0.55))
  
  monte_carlo_sample <- cbind(marginal_covar_sample_cts,
                              marginal_covar_sample_binary) %>%
    ### Adding the transformed covariates we'll use in the model
    .[, c("c_1_tilde", "c_2_tilde", "c_3_tilde", "c_4_tilde") := lapply(.SD, 
                                                                        function(x){(x^2-1)/sqrt(2)}), 
      .SDcols = c("x_1", "x_2", "x_3", "x_4")] %>%
    .[, id := 1:n]
  
  balanced_intercept <- intercepts[ratio == pop_ratio, unique(int_value)]
  
  ### Generating all the true values of the nuisance fcns
  
  monte_carlo_sample[, rho := expit(balanced_intercept + c_1_tilde + c_2_tilde + c_3_tilde), 
                     by = "id"]
  monte_carlo_sample[, mu_11 := (1+1+1)/4*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde+x_5)]
  monte_carlo_sample[, mu_00 := (1+0+0)/4*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde+x_5)]
  monte_carlo_sample[, p_1 := expit((2/5)*(2*1-1+(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)+1*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)))]
  monte_carlo_sample[, p_0 := expit((2/5)*(2*0-1+(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)+0*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)))]
  
  ### Computing the "true" E[Y(1)|U=10, R=0] by applying our identification formulae
  
  monte_carlo_sample[, mean((p_1-p_0)*(1-rho)*mu_11)]/monte_carlo_sample[, mean((p_1-p_0)*(1-rho))]
  
  true_y_1 <- monte_carlo_sample[, mean((p_1-p_0)*(1-rho)*mu_11)/mean((p_1-p_0)*(1-rho))]
  
  true_y_0 <- monte_carlo_sample[, mean((p_1-p_0)*(1-rho)*mu_00)/mean((p_1-p_0)*(1-rho))]
  
  true_te <- true_y_1-true_y_0
  
  return(true_te)
}


### Now, applying the above fcn to compute the true value across the 
### different population ratios. 

### Note that these intercepts are generated in "computing_balancing_intercepts.R"

intercepts <- readRDS(file = "balancing_intercepts.RDS")

set.seed(1)

true_tes <- copy(intercepts)[, true_te := get_truth(n_reps = 10000000, 
                                                    pop_ratio = ratio), 
                             by = "ratio"] %>%
  .[, .(ratio, true_te)]

saveRDS(object = true_tes,
        file = "true_te_values.RDS")



### Computing the true variance of the uncentered EIF:

get_true_var <- function(n_reps, pop_ratio){
  
  
  n <- n_reps
  
  marginal_covar_sample_cts <- data.table(x_1 = runif(n = n, min = -2, max = 2),
                                          x_2 = runif(n = n, min = -2, max = 2),
                                          x_3 = runif(n = n, min = -2, max = 2),
                                          x_4 = runif(n = n, min = -2, max = 2))
  
  marginal_covar_sample_binary <- data.table(x_5 = rbinom(n = n, 
                                                          size = 1, 
                                                          prob = 0.55))
  
  monte_carlo_sample <- cbind(marginal_covar_sample_cts,
                              marginal_covar_sample_binary) %>%
    ### Adding the transformed covariates we'll use in the model
    .[, c("c_1_tilde", "c_2_tilde", "c_3_tilde", "c_4_tilde") := lapply(.SD, 
                                                                        function(x){(x^2-1)/sqrt(2)}), 
      .SDcols = c("x_1", "x_2", "x_3", "x_4")] %>%
    .[, id := 1:n]
  
  balanced_intercept <- intercepts[ratio == pop_ratio, unique(int_value)]
  
  ### Generating all the true values of the nuisance fcns
  
  monte_carlo_sample[, rho := expit(balanced_intercept + c_1_tilde + c_2_tilde + c_3_tilde), 
                     by = "id"]
  monte_carlo_sample[, mu_11 := (1+1+1)/4*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde+x_5)]
  monte_carlo_sample[, mu_00 := (1+0+0)/4*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde+x_5)]
  monte_carlo_sample[, p_1 := expit((2/5)*(2*1-1+(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)+1*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)))]
  monte_carlo_sample[, p_0 := expit((2/5)*(2*0-1+(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)+0*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)))]
  
  monte_carlo_sample[, R := rbinom(n = 1, size = 1, prob = rho), by = "id"]
  
  ### Assigning treatment in the trial sample
  
  monte_carlo_sample[, pi := expit(2/5*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde))]
  monte_carlo_sample[R==1, A := rbinom(n = 1, size = 1, prob = pi), by = "id"]
  
  ### Generating observed compliance behavior
  
  monte_carlo_sample[, p := expit((2/5)*(2*A-1+(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)+A*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)))]
  monte_carlo_sample[R==1, C := rbinom(n = 1, size = 1, prob = p), by = "id"]
  
  ### Generating outcomes
  
  monte_carlo_sample[, mu := (1+A+C)/4*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde+x_5)]
  monte_carlo_sample[R==1, Y := mu+rnorm(1), by = "id"]
  
  ### Computing the uncentered EIF for each person (this is to get a monte
  ### carlo estimate of the true variance of the EIF):
  
  monte_carlo_sample[R == 0, `:=`(A=0, C=0, Y=0)]
  
  monte_carlo_sample[, rho_resid_denom := (p_1-p_0)*((1-R-(1-rho)))]
  monte_carlo_sample[, p_1_resid_denom := (1-rho)*((A*R*(C-p_1))/(pi*rho))]
  monte_carlo_sample[, p_0_resid_denom := (1-rho)*(((1-A)*R*(C-p_0)))/((1-pi)*rho)]
  monte_carlo_sample[, total_resid_denom := (p_1-p_0)*(1-rho)-mean((p_1-p_0)*(1-rho))]
  monte_carlo_sample[, eif_denom := rho_resid_denom+p_1_resid_denom+p_0_resid_denom+total_resid_denom]
  
  monte_carlo_sample[, rho_11_resid := (p_1-p_0)*(mu_11)*((1-R)-(1-rho))]
  monte_carlo_sample[, mu_11_resid := (1-rho)*(p_1-p_0)*((A*C*R*(Y-mu_11))/(pi*rho*p_1))]
  monte_carlo_sample[, p_1_11_resid := (1-rho)*(mu_11)*((A*R*(C-p_1))/(pi*rho))]
  monte_carlo_sample[, p_0_11_resid := (1-rho)*(mu_11)*(((1-A)*R*(C-p_0)))/((1-pi)*rho)]
  monte_carlo_sample[, total_11_resid := (1-rho)*mu_11*(p_1-p_0)-mean((1-rho)*mu_11*(p_1-p_0))]
  monte_carlo_sample[, eif_num_11 := rho_11_resid+mu_11_resid+p_1_11_resid-p_0_11_resid+total_11_resid]
  
  monte_carlo_sample[, rho_00_resid := (p_1-p_0)*(mu_00)*((1-R)-(1-rho))]
  monte_carlo_sample[, mu_00_resid := (1-rho)*(p_1-p_0)*(((1-A)*(1-C)*R*(Y-mu_00))/((1-pi)*(1-p_0)*rho))]
  monte_carlo_sample[, p_1_00_resid := (1-rho)*(mu_00)*((A*R*(C-p_1))/(pi*rho))]
  monte_carlo_sample[, p_0_00_resid := (1-rho)*(mu_00)*(((1-A)*R*(C-p_0)))/((1-pi)*rho)]
  monte_carlo_sample[, total_00_resid := (1-rho)*mu_00*(p_1-p_0)-mean((1-rho)*mu_00*(p_1-p_0))]
  monte_carlo_sample[, eif_num_00 := rho_00_resid+mu_00_resid+p_1_00_resid-p_0_00_resid+total_00_resid]
  
  ### Great! Now we can construct the uncentered influence function for each person:
  
  true_ratio <- monte_carlo_sample[, mean((p_1-p_0)*(1-rho)*(mu_11-mu_00))/mean((p_1-p_0)*(1-rho))]
  true_denom <- monte_carlo_sample[, mean((p_1-p_0)*(1-rho))]
  
  monte_carlo_sample[, eif_num := eif_num_11-eif_num_00]
  
  monte_carlo_sample[, uncentered_eif := (1/true_denom)*(eif_num)-(true_ratio/true_denom)*eif_denom]
  
  true_var <- monte_carlo_sample[, mean(uncentered_eif^2)]
  
  return(true_var)
}

set.seed(1)

true_vars <- copy(intercepts)[, true_var := get_true_var(n_reps = 10000000, 
                                                         pop_ratio = ratio), 
                              by = "ratio"] %>%
  .[, .(ratio, true_var)]

saveRDS(object = true_vars,
        file = "true_var_values.RDS")







