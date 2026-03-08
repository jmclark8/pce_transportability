### This script is not run as part of any given simulation run. It accomplishes
### a stand-alone task of numerically solving for the intercept of the logistic
### regression model for trial participation that reflects our desired marginal
### proportion of target vs. trial participants in the combined population. 

### It relies heavily on the clarifying work of Sarah E. Robertson, Jon A. 
### Steingrimsson, and Issa J. Dahabreh in "Using Numerical Methods to Design 
### Simulations: Revisiting the Balancing Intercept" published in AJE in 
### 2022. 

### Setting up the functions to apply the Robertson et al. method; these
### are taken from the GitHub repository associated with that work 
### (https://github.com/serobertson/BalancingInterceptSolver)

expit <- function(x){
  exp(x)/(1+exp(x))
}


determine_beta <- function(upper = upper_bound, lower = lower_bound)
{
  t_beta <- ((upper + lower) / 2)
  return(t_beta)
}


####################
determine_intercept <- function(intercept = NULL, beta_vec = beta_vec, marg = 0.5, 
                                tolerance = .001, 
                                lower_bound=-10, upper_bound=10, external_dataset=NA)
{
  if(is.null(intercept))
  {
    
    temp_dataset<-external_dataset
    beta_vec <- c((upper_bound+lower_bound)/2, beta_vec)
    count <- 1
    
    #calculate initial p-tilde
    linearformula=paste(beta_vec, colnames(temp_dataset), sep="*")
    linearformula2=paste(linearformula, collapse=" + ")
    lp <- with(temp_dataset, eval(parse(text = linearformula2)))
    mean.oc<-mean(expit(lp))
    
    while(abs(mean.oc - marg) > tolerance && upper_bound>lower_bound && count<=99)
    {
      print(paste("a:",lower_bound, "b:", upper_bound))
      
      if(mean.oc > marg)
      {
        upper_bound <- beta_vec[1]
      } else if(mean.oc < marg)
      {
        lower_bound <- beta_vec[1]
      }
      temp.coefficients <- beta_vec
      beta_vec[1] <- determine_beta(upper_bound, lower_bound)
      print(paste("iteration:",count))
      print(paste("intercept value:", beta_vec[1]))
      linearformula=paste(beta_vec, colnames(temp_dataset), sep="*")
      linearformula2=paste(linearformula, collapse=" + ")
      lp <- with(temp_dataset, eval(parse(text = linearformula2)))
      mean.oc<-mean(expit(lp))
      cat("Outcome marginal expectation:", mean.oc, "\n")
      count <- count + 1      
    }
  } else {
    
    beta_vec <- c(intercept, beta_vec)
    temp.coefficients <- beta_vec
    
    temp_dataset<-external_dataset
    sample_size=nrow(temp_dataset)
    
    linearformula=paste(beta_vec, colnames(temp_dataset), sep="*")
    linearformula2=paste(linearformula, collapse=" + ")
    lp <- with(temp_dataset, eval(parse(text = linearformula2)))
    #use p-hat to evaluate:
    Ysim <- rbinom(sample_size, 1, expit(lp))
    mean.oc<-mean(Ysim)
    cat("Outcome marginal expectation:", mean.oc, "\n")
    count=0 #no solving steps
  }
  #only print this warning when solving for an intercept
  if(abs(mean.oc - marg) > tolerance & is.null(intercept)==T){warning('Algorithm did not converge within tolerance.')}
  return(data.frame(intercept=beta_vec[1], marginal_expectation=mean.oc, iterations=count, 
                    a=lower_bound, b=upper_bound))
  
}

### This is a tiny extension on the above that we'll use when actually computing
### the intercept.

get_int <- function(x){
  
  determine_intercept(external_data = ext_data, 
                      tolerance = .0001,
                      marg = x, 
                      ### Setting the beta coefficients according to 
                      ### our data generating model.
                      beta_vec = c(1, 1, 1),
                      intercept=NULL,
                      lower_bound=-10, upper_bound=10) %>%
    .$`intercept`
  
}


# Generating data for the monte carlo approximation -----------------------

library(MASS)
library(data.table)
library(magrittr)
library(here)
library(truncnorm)

set.seed(123)

n <- 10^6

int <- rep(1, n)

marginal_covar_sample_cts <- data.table(x_1 = runif(n = n, min = -2, max = 2),
                                        x_2 = runif(n = n, min = -2, max = 2),
                                        x_3 = runif(n = n, min = -2, max = 2),
                                        x_4 = runif(n = n, min = -2, max = 2))

### Adding the transformed covariates we'll use in the linear predictor

marginal_covar_sample_cts[, c("c_1_tilde", "c_2_tilde", 
                              "c_3_tilde", "c_4_tilde") := lapply(.SD, 
                                                                  function(x){(x^2-1)/sqrt(2)}), 
                          .SDcols = c("x_1", "x_2", "x_3", "x_4")]

### Adding an intercept to the design matrix

marginal_covar_sample_cts[, int := 1]

ext_data <- marginal_covar_sample_cts[, .(int, c_1_tilde, 
                                          c_2_tilde, c_3_tilde)] %>%
  data.frame()



# Solving for the intercept in each simulation setting --------------------

### We're gonna do this by making a table with each of the various ratios and
### then solving for the intercept by ratio

int_res <- data.table(ratio = c(500/(500+10000), 
                                500/(500+1000),
                                500/(500+500))) %>%
  ### Computing the intercept for each scenario. 
  .[, int_value := get_int(ratio), by = "ratio"] %>%
  .[, .(ratio, int_value)]

### Cool! 

### Saving these for later:

saveRDS(object = int_res,
        file = "balancing_intercepts.RDS")





































