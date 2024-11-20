# Created: 2024.11.19, Ward B. Eiling
start_time <- Sys.time() # time the script

library(foreach)
library(doParallel)
library(data.table)

# Load the functions
source("simulation/scripts/functions/generative_models_function.R")
source("simulation/scripts/functions/model_fitting_function.R")

# Parameters for the simulation
num_reps <- 1000  # Adjust as needed
set.seed(120)

# Set up parallel backend
cores <- detectCores()
cl <- makeCluster(cores - 1)
registerDoParallel(cl)

# Create combinations of settings
settings <- expand.grid(sample_size = c(30, 100, 200), total_T = c(10, 30), dgm_type = 1:3)

# Run simulation and model fitting
sim_results <- foreach(setting = iter(settings, by = "row"), .combine = "rbind", .packages = c("lme4", "geepack", "data.table")) %dopar% {
  sample_size <- setting$sample_size
  total_T <- setting$total_T
  dgm_type <- setting$dgm_type
  
  # Simulate data
  simulated_data <- replicate(num_reps, dgm_with_treatment(sample_size, total_T, dgm_type), simplify = FALSE)
  
  # Fit models to each dataset
  lapply(simulated_data, function(data) {
    fit_models(data, gm_type = dgm_type)
  })
}

# Stop the cluster
stopCluster(cl)

# Save results
saveRDS(sim_results, "simulation/output/simulation_results_1000reps.rds")

### Data Extraction ###

sim_results <- readRDS("simulation/output/simulation_results_1000reps.rds")

# Create function that retrieves the estimates from the model results
num_reps <- 20

# Initialize a list to store the averages for each condition
averages <- list()

conditions <- c("result.1", "result.2", "result.3", "result.4", "result.5", "result.6", "result.7", "result.8", "result.9", 
                "result.10", "result.11", "result.12", "result.13", "result.14", "result.15", "result.16", "result.17", "result.18")
models <- c("mlm", "gee_ind", "gee_exch", "gee_ar1")

# Loop over all conditions
for (condition in conditions) {
  # Extract values for the current condition from all iterations
  cond_values <- sim_results[condition,]
  unl_cond_values <- unlist(cond_values)
  
  # select elements for each model
  mlm_est <- unl_cond_values[grepl("mlm", names(unl_cond_values))]
  gee_ind_est <- unl_cond_values[grepl("gee_ind", names(unl_cond_values))]
  gee_exch_est <- unl_cond_values[grepl("gee_exch", names(unl_cond_values))]
  gee_ar1_est <- unl_cond_values[grepl("gee_ar1", names(unl_cond_values))]
  
  # Compute the averages for each model
  averages[[condition]] <- c(mlm = mean(mlm_est, na.rm = TRUE),
                             gee_ind = mean(gee_ind_est, na.rm = TRUE),
                             gee_exch = mean(gee_exch_est, na.rm = TRUE),
                             gee_ar1 = mean(gee_ar1_est, na.rm = TRUE))
}

# View the averages
# averages

# Convert the averages to a dataframe
results_df <- t(as.data.frame(averages))

# Add setting information to the dataframe "setting"
results_df <- cbind(settings, results_df)
results_df <- results_df[, c(3, 2, 1, 4, 5, 6, 7)]
colnames(results_df) <- c("GM", "T", "N", "mlm", "gee_ind", "gee_exch", "gee_ar1")

# save as RDS file
saveRDS(results_df, "simulation/output/simulation_results_table_1000reps.rds")

# compute bias (treatment effect beta_20 = 1)
results_df_bias <- results_df
results_df_bias$mlm <- results_df$mlm - 1
results_df_bias$gee_ind <- results_df$gee_ind - 1
results_df_bias$gee_exch <- results_df$gee_exch - 1
results_df_bias$gee_ar1 <- results_df$gee_ar1 - 1

saveRDS(results_df_bias, "simulation/output/simulation_results_table_bias_1000reps.rds")
  
### Create LaTeX Table ###

library(xtable)

# Create a table with the results of raw estimates
table_results <- xtable(results_df, caption = "Simulation Results of treatment effect ($beta_{20}$) estimate for different estimation methods, 1000 replications", 
                        label = "tab:sim_results", digits = c(0, 0, 0, 0, 3, 3, 3, 3))
print(table_results, include.rownames = FALSE, hline.after = c(-1, 0, seq(from = 6, to = nrow(table_results), by = 6)), 
                                                               file = "simulation/output/simulation_results_table_1000reps.tex")

# Create a table with the results of bias estimates
table_results_bias <- xtable(results_df_bias, caption = "Simulation Results of bias in treatment effect ($beta_{20}$) estimate for different estimation methods, 1000 replications", 
                        label = "tab:sim_results_bias", digits = c(0, 0, 0, 0, 3, 3, 3, 3))
print(table_results_bias, include.rownames = FALSE, hline.after = c(-1, 0, seq(from = 6, to = nrow(table_results_bias), by = 6)), 
                                                               file = "simulation/output/simulation_results_table_bias_1000reps.tex")

end_time <- Sys.time() # end time
end_time - start_time # ~ 55 minutes with 11 cores
