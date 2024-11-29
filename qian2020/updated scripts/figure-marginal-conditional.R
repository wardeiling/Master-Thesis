# This script generates a figure that shows how marginal and conditional effects may not 
# be equivalent under certain conditions (GM3) in a multilevel linear model.

library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(gee)

# set the seed
set.seed(123)

# Create a function to generate the data
generate_data <- function(n_i = 5000, sigma_u = 1, sigma_Z1 = 1, sigma_e = 1, beta_1 = 0.8, gamma_00 = 2) {
  data_list <- list()
  
  # Simulate data for each individual
  for (i in 1:n_i) {
    u_0i <- rnorm(1, 0, sigma_u)
    Z_i2_lag <- Z_i1 <- rnorm(1, 0, sigma_Z1)
    Y_i2 <- rnorm(1, gamma_00 + beta_1 * Z_i1 + u_0i, sigma_e)
    Z_i3_lag <- Z_i2 <- Y_i2
    Y_i3 <- rnorm(1, gamma_00 + beta_1 * Z_i2 + u_0i, sigma_e)
    
    # Store the data in a list
    subject_data <- data.frame(
      id = i,
      time = 1:2,
      Y = c(Y_i2, Y_i3),
      Z_lag1 = c(Z_i2_lag, Z_i3_lag)
    )
    
    data_list[[i]] <- subject_data
  }
  
  # Combine all subjects' data into a single data frame
  data_long <- do.call(rbind, data_list)
  return(data_long)
}

# Create a function that runs the simulations
run_simulations <- function(n_sim, n_i, sigma_u, sigma_Z1, sigma_e, beta_1, gamma_00) {
  
  # Initialize a list to store results
  all_estimates <- list()
  
  # Simulation loop
  for (sim in 1:n_sim) {
    # Generate the data
    data_sim <- generate_data(n_i, sigma_u, sigma_Z1, sigma_e, beta_1, gamma_00)
    
    # Fit the models
    mlm_mle <- lmer(Y ~ Z_lag1 + (1 | id), data = data_sim, REML = FALSE)
    gee_exch <- gee(Y ~ Z_lag1, id = id, data = data_sim, family = gaussian, corstr = "exchangeable")
    gee_ind <- gee(Y ~ Z_lag1, id = id, data = data_sim, family = gaussian, corstr = "independence")
    gee_ar1 <- gee(Y ~ Z_lag1, id = id, data = data_sim, family = gaussian, corstr = "AR-M", Mv = 1)
    glm <- glm(Y ~ Z_lag1, data = data_sim, family = gaussian)
    # gls_symm <- nlme::gls(Y ~ Z_lag1, data = data_sim, correlation = corSymm(form = ~ 1 | id), method = "ML")
    
    # Extract the fixed effect estimates for each model
    estimates <- data.frame(
      sim = sim,
      MLM_mle_intercept = fixef(mlm_mle)[1],
      MLM_mle_slope = fixef(mlm_mle)[2],
      GEE_exch_intercept = coef(gee_exch)[1],
      GEE_exch_slope = coef(gee_exch)[2],
      GEE_ind_intercept = coef(gee_ind)[1],
      GEE_ind_slope = coef(gee_ind)[2],
      GEE_ar1_intercept = coef(gee_ar1)[1],
      GEE_ar1_slope = coef(gee_ar1)[2],
      GLM_intercept = coef(glm)[1],
      GLM_slope = coef(glm)[2] #,
      #GLS_symm_intercept = coef(gls_symm)[1],
      #GLS_symm_slope = coef(gls_symm)[2]
    )
    
    # Store the estimates in the list
    all_estimates[[sim]] <- estimates
  }
  
  # Combine results into a single data frame
  results <- do.call(rbind, all_estimates)
  return(results)
}

# Set the parameters for the simulation
n_sim <- 1  # Number of simulations
n_i <- 100000   # Number of individuals per simulation
sigma_u <- 0  # Variance of random intercept
sigma_Z1 <- 1 # Variance of Z1
sigma_e <- 1 # Residual variance
beta_1 <- 0.8  # Slope
gamma_00 <- 2  # Intercept

# Run the simulation and store the results
simulation_results <- run_simulations(n_sim, n_i, sigma_u, sigma_Z1, sigma_e, beta_1, gamma_00)

# Calculate mean and standard deviation of estimates
section2.2_summary_stats <- data.frame(
  row.names = c("Intercept", "Z_lag1"),
  MLM_mle = c(mean(simulation_results$MLM_mle_intercept), mean(simulation_results$MLM_mle_slope)),
  GEE_exch = c(mean(simulation_results$GEE_exch_intercept), mean(simulation_results$GEE_exch_slope)),
  GEE_ind = c(mean(simulation_results$GEE_ind_intercept), mean(simulation_results$GEE_ind_slope)),
  GEE_ar1 = c(mean(simulation_results$GEE_ar1_intercept), mean(simulation_results$GEE_ar1_slope)),
  GLM = c(mean(simulation_results$GLM_intercept), mean(simulation_results$GLM_slope))) 

### Step 1. Obtain data for 10 different individuals with increasing random effects

p1 <- generate_data(n_i = 10, sigma_u = 1, sigma_Z1 = 1, sigma_e = 1, beta_1 = 0.8, gamma_00 = 2)

# Create a data frame with the results
results_long <- gather(simulation_results, parameter, estimate, -sim)

### Step 2. Obtain the conditional average: the subject with b_i = 0, the mean individual

# when we fill in sigma_u = 0, the random effect is 0 for all individuals.
# The estimated effect should therefore represent the conditional average.

cond_slope <- 0.8
cond_intercept <- 2

### Step 3. Obtain the population average: the average Y for every level of X over all individuals

marg_slope <- 0.98
marg_intercept <- 1.82



