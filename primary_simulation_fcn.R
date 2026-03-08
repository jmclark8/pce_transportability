### This is the primary function to run the expanded simulation for our paper.
### It takes in the parameters relevant to a simulation setting
### and spits out the result for that setting.

### Key params:

### 1. Model Mispecification (5 settings)
### -----------------------------------------------
#### This is specified via four arguments, `tp`, `om`, `pp`, and `ps`.
#### The abbreviations used here are "om" for outcome model, 
#### "tp" for treatment probability, "ps" for principal score, 
#### and "pp" for participation probability. Each of them takes value
#### TRUE or FALSE if we'd like the corresponding nuisance function to be
#### specified correctly.

### 2. Sample Sizes (9 settings)
### -----------------------------------------------
#### This is specified via two arguments: trial_size and target_size. They each
#### take numeric arguments that specify the size of the trial and target sample
#### sizes, respectively.

library(data.table)
library(magrittr)
library(MASS)
library(truncnorm)
# 
# 
intercepts <- readRDS("balancing_intercepts.RDS")

source("sim_covariates.R")

source("estimate_nuisance_fcns.R")

source("estimate_tes.R")

true_te_values <- readRDS("true_te_values.RDS")

true_var_values <- readRDS("rue_var_values.RDS")

expit <- function(x){
  exp(x)/(1+exp(x))
}

upper_lower <- function(x, alpha  = 0.05){
  data.frame(lwr = quantile(x, probs = alpha/2, names = FALSE),
             upper = quantile(x, probs = 1-alpha/2, names = FALSE))
}

### This is the primary simulation function to generate the results for a specific
### set of correctly specified nuisance functions and sample sizes. Please see
### below this function for a full specification of our parameter grid.

get_sim_res <- function(tp, om, pp, ps, 
                        trial_size, target_size){
  
  # Data Generation -----------------------------------------------------
  
  ### The basic structure here is that we're gonna iteratively build the simulation
  ### run's dataset `sample_data` by sequentially adding certain data elements.
  
  ### These elements will ultimately include estimated nuisance fcns.
  
  ### Generating covariates and participation indicator
  
  total_n <- trial_size + target_size
  ratio <- trial_size/(trial_size+target_size)
  
  sample_data <- sim_pop_samples(total_n = total_n, 
                                 ratio = ratio)
  
  ### Assigning treatment in the trial sample
  
  sample_data[R==1, pi := expit(2/5*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde))]
  sample_data[R==1, A := rbinom(n = 1, size = 1, prob = pi), by = "id"]
  
  ### Generating observed compliance behavior
  
  sample_data[R==1, p := expit((2/5)*(2*A-1+(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)+A*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde)))]
  sample_data[R==1, C := rbinom(n = 1, size = 1, prob = p), by = "id"]
  
  ### Generating outcomes
  
  sample_data[R==1, mu := (1+A+C)/4*(c_1_tilde+c_2_tilde+c_3_tilde+c_4_tilde+x_5)]
  sample_data[R==1, Y := mu+rnorm(1), by = "id"]
  
  ### Classifying each row in the sample data for cross-fitting purposes
  
  n_folds <- 10
  
  sample_data[, cf_fold := which(rmultinom(n = 1, size = 1, 
                                           prob = rep(1/n_folds, n_folds))==1), 
              by = "id"] %>%
    .[, cf_fold := factor(cf_fold)]
  
  # Estimating Nuisance Functions ----------------------------------------
  
  ### We construct two sets of estimates for each nuisance function: one
  ### set that's estimated using data only in the row's fold and another set
  ### estimated from the entire dataset. Note that we construct both such
  ### estimates for every simulation run. 
  
  ### Principal Scores
  ps_estimates <- estimate_ps(data = sample_data, ps = ps, 
                              n_folds = n_folds)
  
  ### Outcome Model
  om_estimates <- estimate_om(data = sample_data, om = om, 
                              n_folds = n_folds)
  
  ### Participation Probability
  pp_estimates <- estimate_pp(data = sample_data, pp = pp, 
                              n_folds = n_folds)
  
  ### Treatment Probability
  tp_estimates <- estimate_tp(data = sample_data, tp = tp, 
                              n_folds = n_folds)
  
  ### Merging these estimates back into the larger dataset:
  
  sample_data <- list(sample_data, ps_estimates, om_estimates, pp_estimates, 
                      tp_estimates) %>%
    Reduce(function(x, y){merge(x, y, by = "id")}, x = .)
  
  # Estimating treatment effects --------------------------------------------
  
  ### Using the Efficient Influence Function (Note that this also returns EIF-based
  ### variance estimators)
  eif_te_estimates <- estimate_eif_te(data = sample_data)
  
  ### Using the IPW approach
  ipw_te_estimates <- estimate_ipw_te(data = sample_data)
  
  ### Using the OM approach
  om_te_estimates <- estimate_om_te(data = sample_data)
  
  ### Using a pure plug-in approach
  
  plug_in_te_estimates <- estimate_plug_in_te(data = sample_data)
  
  ### Putting them all together:
  
  te_estimates <- lapply(list(eif_te_estimates, ipw_te_estimates, 
                              om_te_estimates, plug_in_te_estimates), transpose, keep.names = "rn") %>%
    rbindlist() %>%
    setnames(old = c("rn", "V1"), 
             new = c("te_est_type", "te_est"))
  
  return(te_estimates)
  
}

### Parameter grid:

all_true <- data.table(tp = TRUE, om = TRUE, pp = TRUE, ps = TRUE)
tp_pp_ps <- data.table(tp = TRUE, om = FALSE, pp = TRUE, ps = TRUE)
tp_pp_om <- data.table(tp = TRUE, om = TRUE, pp = TRUE, ps = FALSE)
om_ps <- data.table(tp = FALSE, om = TRUE, pp = FALSE, ps = TRUE)
all_wrong <- data.table(tp = FALSE, om = FALSE, pp = FALSE, ps = FALSE)

sample_size_settings <- data.table(ratio = c(500/(500+10000), 
                                             500/(500+1000), 
                                             500/(500+500))) %>%
  .[, small_trial_size := 500] %>%
  .[, small_target_size := (small_trial_size-small_trial_size*ratio)/ratio] %>%
  ### Specifying other sample sizes
  .[, `:=` (med_trial_size = small_trial_size*(2)^1, 
            med_target_size = small_target_size*(2)^1)] %>%
  .[, `:=` (large_trial_size = small_trial_size*(2)^2, 
            large_target_size = small_target_size*(2)^2)] %>%
  .[, `:=` (asymp_trial_size = small_trial_size*100, 
            asymp_target_size = small_target_size*100)] %>%
  pivot_longer(cols = ends_with("_size"), 
               names_to = c("size_cat", "source_type"), 
               names_pattern = "(^[a-z]+)_(trial_size|target_size)$") %>%
  pivot_wider(names_from = "source_type", values_from = "value") %>%
  dplyr::select(c(trial_size, target_size)) %>%
  data.table()

param_grid <- rbindlist(list(cbind(sample_size_settings, all_true), 
                             cbind(sample_size_settings, tp_pp_ps), 
                             cbind(sample_size_settings, tp_pp_om), 
                             cbind(sample_size_settings, om_ps), 
                             cbind(sample_size_settings, all_wrong))) %>%
  data.frame()





























