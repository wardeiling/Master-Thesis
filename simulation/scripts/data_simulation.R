# Created: 2024.11.19, Ward B. Eiling

library(foreach)
library(doParallel)

# Load the generative model function from the file
source("simulation/scripts/functions/generative_models_function.R")
source("simulation/scripts/functions/model_fitting_function.R")

# Parameters for the simulation
sample_sizes <- c(30, 100, 200)
num_obs <- c(10, 30)
dgm_types <- 1:3
num_reps <- 20 # increase later
set.seed(120)

# Set up parallel backend
cores <- detectCores()
cl <- makeCluster(cores - 1)  # Leave one core free
registerDoParallel(cl)

# Create combinations of settings
settings <- expand.grid(sample_size = sample_sizes, total_T = num_obs, dgm_type = dgm_types)

# Initialize a list to store results
sim_data <- foreach(setting = iter(settings, by = "row"), .combine = "rbind", .packages = "data.table") %dopar% {
  sample_size <- setting$sample_size
  total_T <- setting$total_T
  dgm_type <- setting$dgm_type
  
  # Run the simulation for this setting
  replicate(num_reps, dgm_with_treatment(sample_size, total_T, dgm_type), simplify = FALSE)
}

# Stop the cluster
stopCluster(cl)

# Optionally save the results
saveRDS(sim_data, "simulation/output/simulation_data_only_TE.rds")

### Data Extraction

# Load the results
sim_data <- readRDS("simulation/output/simulation_data_only_TE.rds")

# Extract the data
example <- sim_data[[1]]
