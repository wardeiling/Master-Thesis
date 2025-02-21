rm(list = ls()) # clear workspace

# load libraries
library(lme4) # for generalized linear mixed models
library(geepack) # for generalized estimating equations
library(foreach) # for parallelization
library(doParallel) # for parallelization
library(doRNG) # for reproducibility

# load helper functions
source("scripts/hamaker2020_extended/helper-functions/data-generation.R")
source("scripts/hamaker2020_extended/helper-functions/model-fitting.R")
source("scripts/hamaker2020_extended/helper-functions/result-formatting.R")

### WRAPPER FUNCTION ###
run_simulation <- function(runname = "run1", seed = 4243, nsim = 1000, N_total = 200, T_total = 10, 
                           predictor.type = "continuous", outcome.type = "continuous",
                           sdX.within = sqrt(1), sdX.between = sqrt(4), 
                           g.00 = 0, g.01 = 2, sd.u0 = 1, g.10 = 1, 
                           sd.u1 = 0, sd.e = 1) {
  
  # runname more information
  runname_upd <- paste0(runname, "_pred", predictor.type, "_out", outcome.type, "_Nrep", nsim)
  
  # Create results directory
  dir.create(paste0("simulation_results_glmm/", runname_upd), 
                    showWarnings = FALSE, recursive = TRUE)
  
  # set global seed
  set.seed(123) 
  
  # Set up parallel backend
  cores <- detectCores()
  cl <- makeCluster(cores - 1, outfile = paste0("simulation_results_glmm/", runname_upd, "/log.txt"))
  registerDoParallel(cl)
  
  # Export functions to workers
  clusterExport(cl, varlist = c("glmm_data_generation", "glmm_model_fitting", "glmm_formating_results"))
  
  # Run simulations in parallel
  parallel_results <- foreach(isim = 1:nsim, .packages = c("geepack", "lme4"), .options.RNG=120, .verbose = TRUE) %dorng% {
    if (isim %% 10 == 0) {
      cat(paste("Starting iteration", isim, "\n"))
    }
    
    # Generate data
    data <- glmm_data_generation(N_total = N_total, T_total = T_total, 
                                 predictor.type = predictor.type, outcome.type = outcome.type,
                                 sdX.within = sdX.within, sdX.between = sdX.between, 
                                 g.00 = g.00, g.01 = g.01, sd.u0 = sd.u0, 
                                 g.10 = g.10, sd.u1 = sd.u1, sd.e = sd.e)
    
    # Fit models
    models <- glmm_model_fitting(data, outcome.type = outcome.type)
    
    # Format results
    result_table <- glmm_formating_results(models)
    
    return(result_table)
  }
  
  # Stop parallel backend
  stopCluster(cl)
  
  # Compute mean results
  mean_results <- Reduce("+", parallel_results) / nsim
  
  # Compute variance (unbiased estimator)
  var_results <- Reduce("+", lapply(parallel_results, function(x) (x - mean_results)^2)) / (nsim - 1)
  sd_results <- sqrt(var_results)
  monte_carlo_se <- sd_results / sqrt(nsim)
  
  settings <- list(runname = runname, seed = seed, nsim = nsim, N_total = N_total, T_total = T_total, 
                   predictor.type = predictor.type, outcome.type = outcome.type,
                   sdX.within = sdX.within, sdX.between = sdX.between, 
                   g.00 = g.00, g.01 = g.01, sd.u0 = sd.u0, g.10 = g.10, 
                   sd.u1 = sd.u1, sd.e = sd.e)
  
  output <- list(settings = settings, all_results = parallel_results, mean_results = mean_results, sd_results = sd_results, monte_carlo_se = monte_carlo_se)
  for (name in names(output)) {
    saveRDS(output[[name]], file = paste0("simulation_results_glmm/", runname_upd, "/", name, ".rds"))
  }
  
  return(output)
}

### RUNNING CONDITIONS ###

# stable specifications across simulations
seed = 4243
nsim = 1000
N_total = 200
T_total = 10

sdX.within = sqrt(1)
g.00 = 0
g.01 = 2
sd.u0 = 1
g.10 = 1
sd.u1 = 0

# varying specifications across simulations
sdX.between.continuous = sqrt(4)
sdX.between.binary = sqrt(0.4) # to ensure that we have reallistic values (not some people only having 0 or 1)
sd.e.continuous = 1
sd.e.binary = 0 # irrelevant for binary outcomes

# contxy_sim <- run_simulation(runname = "run4", seed = seed, nsim = nsim, N_total = N_total, T_total = T_total,
#                            predictor.type = "continuous", outcome.type = "continuous",
#                            sdX.within = sdX.within, sdX.between = sdX.between.continuous,
#                            g.00 = g.00, g.01 = g.01, sd.u0 = sd.u0, g.10 = g.10,
#                            sd.u1 = sd.u1, sd.e = sd.e.continuous)
# 
# contxy_sim$mean_results
# contxy_sim$monte_carlo_se
#
# observations
# - once we increase T_total, the total effect is comprised more of the within-person effect, which explains
#   the similarity between the "uninterpretable blend" of raw X and the within-person effects.

# binx_conty_sim <- run_simulation(runname = "run4", seed = seed, nsim = nsim, N_total = N_total, T_total = T_total, 
#                            predictor.type = "binary", outcome.type = "continuous",
#                            sdX.within = sdX.within, sdX.between = sdX.between.binary, 
#                            g.00 = g.00, g.01 = g.01, sd.u0 = sd.u0, g.10 = g.10, 
#                            sd.u1 = sd.u1, sd.e = sd.e.continuous)
# 
# binx_conty_sim$mean_results
# binx_conty_sim$monte_carlo_se

# contx_biny_sim <- run_simulation(runname = "run1",seed = seed, nsim = nsim, N_total = N_total, T_total = T_total, 
#                            predictor.type = "continuous", outcome.type = "binary",
#                            sdX.within = sdX.within, sdX.between = sdX.between.continuous, 
#                            g.00 = g.00, g.01 = g.01, sd.u0 = sd.u0, g.10 = g.10, 
#                            sd.u1 = sd.u1, sd.e = sd.e.binary)
# 
# contx_biny_sim$mean_results
# contx_biny_sim$monte_carlo_se

binaryxy_sim <- run_simulation(runname = "run1", seed = seed, nsim = nsim, N_total = N_total, T_total = T_total, 
                           predictor.type = "binary", outcome.type = "binary",
                           sdX.within = sdX.within, sdX.between = sdX.between.binary, 
                           g.00 = g.00, g.01 = g.01, sd.u0 = sd.u0, g.10 = g.10, 
                           sd.u1 = sd.u1, sd.e = sd.e.binary)

binaryxy_sim$mean_results
binaryxy_sim$monte_carlo_se
