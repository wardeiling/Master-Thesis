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
    
    return(models)
  }
  
  # Stop parallel backend
  stopCluster(cl)
  
  # Save all raw results
  saveRDS(parallel_results, file = paste0("simulation_results_glmm/", runname_upd, "/parallel_results.rds"))
  
  # Apply result formatting **only once**
  formatted_results <- lapply(parallel_results, glmm_formating_results)
  
  # Compute mean results
  mean_results <- Reduce("+", formatted_results) / nsim
  
  # Compute variance (unbiased estimator)
  var_results <- Reduce("+", lapply(formatted_results, function(x) (x - mean_results)^2)) / (nsim - 1)
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
# N_total = 200
# T_total = 10

# contxy_sim <- run_simulation(runname = "update_output_format2", seed = seed, nsim = nsim, N_total = 200, T_total = 10, 
#                              predictor.type = "continuous", outcome.type = "continuous",
#                              sdX.within = 0.25, sdX.between = 0.5, g.00 = 0, g.01 = 1, sd.u0 = 0.5,
#                              g.10 = 0.5, sd.u1 = 0, sd.e = 0.5)
# 
# contxy_sim$mean_results
# contxy_sim$monte_carlo_se
#
# observations
# - once we increase T_total, the total effect is comprised more of the within-person effect, which explains
#   the similarity between the "uninterpretable blend" of raw X and the within-person effects.

# initialize timer
start_time <- Sys.time()
binx_conty_sim <- run_simulation(runname = "speed_old", seed = seed, nsim = nsim, N_total = 200, T_total = 10, 
                                 predictor.type = "binary", outcome.type = "continuous",
                                 sdX.within = NA, sdX.between = 0.5, g.00 = 0, g.01 = 1, sd.u0 = 0.5,
                                 g.10 = 0.5, sd.u1 = 0, sd.e = 0.5)

binx_conty_sim$mean_results
binx_conty_sim$monte_carlo_se

# finalize timer
end_time <- Sys.time()
end_time - start_time

# contx_biny_sim <- run_simulation(runname = "run1",seed = seed, nsim = nsim, N_total = 200, T_total = 10, 
#                                  predictor.type = "continuous", outcome.type = "binary",
#                                  sdX.within = 0.25, sdX.between = 0.5, g.00 = 0, g.01 = 1, sd.u0 = 0.5,
#                                  g.10 = 0.5, sd.u1 = 0, sd.e = NA)
# 
# contx_biny_sim$mean_results
# contx_biny_sim$monte_carlo_se
# 
# binaryxy_sim <- run_simulation(runname = "run1", seed = seed, nsim = nsim, N_total = 200, T_total = 10, 
#                                predictor.type = "binary", outcome.type = "binary",
#                                sdX.within = NA, sdX.between = 0.5, g.00 = -0.25, g.01 = 0.5, sd.u0 = 0.5,
#                                g.10 = 0.5, sd.u1 = 0, sd.e = NA)
# 
# binaryxy_sim$mean_results
# binaryxy_sim$monte_carlo_se

### Retrieve results

# output from GLMM with logit link is in log-odds, so to get the probability we need to transform it
# p = exp(logit) / (1 + exp(logit))

# transform every element in matrix mean_results
# apply(mean_results, c(1, 2), function(x) exp(x) / (1 + exp(x)))
