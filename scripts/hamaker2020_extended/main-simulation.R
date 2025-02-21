######################################
#### MODULARIZED SIMULATION SCRIPT ###
######################################

# Author: Ward B Eiling
# Last edited: 21-02-2025

### Simulation set-up -------------------------------------------------------

rm(list = ls()) # clear workspace
set.seed(123) # set global seed
runname <- "run3" # set a runname

# make a directory in simulation_results based on runname
dir.create(paste0("simulation_results_glmm/", runname), showWarnings = FALSE)

# load libraries
library(lme4) # for generalized linear mixed models
library(geepack) # for generalized estimating equations
# library(tidyverse) # for data manipulation
library(foreach) # for parallelization
library(doParallel) # for parallelization
library(doRNG) # for reproducibility

# load helper functions
source("scripts/hamaker2020_extended/helper-functions/data-generation.R")
source("scripts/hamaker2020_extended/helper-functions/model-fitting.R")
source("scripts/hamaker2020_extended/helper-functions/result-formatting.R")

# Set up parallel backend
cores <- detectCores()
cl <- makeCluster(cores - 1, outfile = paste0("simulation_results_glmm/", runname, "/log.txt"))
registerDoParallel(cl)

# set the number of simulations
nsim <- 10

### Simulation ---------------------------------------------------------------

# run the simulation
parallel_results <- foreach(isim = 1:nsim, .packages = c("geepack", "lme4"), .options.RNG=120, .verbose = TRUE) %dorng% {
  if (isim %% 10 == 0) {
    cat(paste("Starting iteration",isim,"\n"))
  }

  data <- glmm_data_generation(N_total = 200, T_total = 10, predictor.type = "continuous", outcome.type = "continuous",
                               sdX.within = sqrt(1), sdX.between = sqrt(4), g.00 = 0, g.01 = 2, sd.u0 = 1, 
                               g.10 = 1, sd.u1 = 0, sd.e = 1)
  models <- glmm_model_fitting(data, outcome.type = "continuous")
  result_table <- glmm_formating_results(models)
}

# stop the parallel backend
stopCluster(cl)

# save results
saveRDS(parallel_results, file = paste0("simulation_results_glmm/", runname, "/parallel_results.rds"))

### Post-processing of Results ----------------------------------------------

# average results across lists
mean_results <- Reduce("+", parallel_results) / nsim

# save results
saveRDS(mean_results, file = paste0("simulation_results_glmm/", runname, "/mean_results.rds"))

# Compute variance (unbiased sample variance)
var_results <- Reduce("+", lapply(parallel_results, function(x) (x - mean_results)^2)) / (nsim - 1)

# save results
saveRDS(var_results, file = paste0("simulation_results_glmm/", runname, "/var_results.rds"))
