######################################
#### MODULARIZED SIMULATION SCRIPT ###
######################################

rm(list = ls()) # clear workspace
set.seed(123) # set global seed
runname <- "run1" # set a runname

# make a directory in simulation_results based on runname
dir.create(paste0("simulation_results_glmm/", runname), showWarnings = FALSE)

# load libraries
library(lme4)
library(gee)
library(geepack)
library(tidyverse)
library(foreach)
library(doParallel)
library(doRNG)

# load helper functions
source("scripts/hamaker2020_extended/helper-functions/data-generation.R")
source("scripts/hamaker2020_extended/helper-functions/model-fitting.R")
source("scripts/hamaker2020_extended/helper-functions/result-formatting.R")

# Set up parallel backend
cores <- detectCores()
cl <- makeCluster(cores - 1, outfile = paste0("simulation_results_glmm/", runname, "/log.txt"))
registerDoParallel(cl)

# set the number of simulations
nsim <- 1000





# MAIN SIMULATION
run_simulation <- function(nsim = 100, N_total = 5000, T_total = 20, predictor.type = "continuous", outcome.type = "continuous",
                           sdX.within = sqrt(1), sdX.between = sqrt(4), g.00 = 0, g.01 = 2, sd.u0 = 1, 
                           g.10 = 1, sd.u1 = 0, sd.e = 1) {
  
  # PARALLELIZATION SETUP
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  result <- foreach(isim = 1:nsim, .packages = c("geepack", "lme4", "tidyverse"), .options.RNG=120) %dorng% {
                      if (isim %% 10 == 0) {
                        cat(paste("Starting iteration",isim,"\n"))
                      }
    
  data <- glmm_data_generation(N_total, T_total, predictor.type, outcome.type, seed, sdX.within, sdX.between, g.00, g.01, sd.u0, g.10, sd.u1, sd.e)
  models <- glmm_model_fitting(data, outcome.type)
  result_table <- format_results(models)
  
  return(result_table)
}
  
  return(result)
}

run_simulation()



# RUN MULTIPLE SIMULATIONS IN PARALLEL
results <- foreach(i = 1:10, .combine = rbind, .packages = c("lme4", "gee", "geepack", "tidyverse")) %dopar% {
  run_simulation(N = 200, n = 10, seed = 3859 + i)
}

stopCluster(cl)
print(results)


# ### EXAMPLE WITH 1 REPLICATION
# data <- glmm_data_generation(N_total = 5000, T_total = 20, predictor.type = "continuous", outcome.type = "continuous",
#                              sdX.within = sqrt(1), sdX.between = sqrt(4), g.00 = 0, g.01 = 2, sd.u0 = 1, 
#                              g.10 = 1, sd.u1 = 0, sd.e = 1)
# models <- glmm_model_fitting(data, outcome.type = "continuous")
# result_table <- glmm_formating_results(models)
