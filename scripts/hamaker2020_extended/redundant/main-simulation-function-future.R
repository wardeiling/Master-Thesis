### WRAPPER FUNCTION ###
run_simulation <- function(runname = "run1", seed = 4243, nsim = 1000, N_total = 200, T_total = 10, 
                           predictor.type = "continuous", outcome.type = "continuous",
                           sdX.within = sqrt(1), sdX.between = sqrt(4), 
                           g.00 = 0, g.01 = 2, sd.u0 = 1, g.10 = 1, 
                           sd.u1 = 0, sd.e = 1) {
  
  # load libraries
  library(lme4) # for generalized linear mixed models
  library(geepack) # for generalized estimating equations
  library(doFuture) # for parallelization (neatly handles RNG, errors, making objects available for workers)
  
  # load helper functions
  source("scripts/hamaker2020_extended/helper-functions/data-generation-centeredX.R")
  source("scripts/hamaker2020_extended/helper-functions/model-fitting.R")
  source("scripts/hamaker2020_extended/helper-functions/result-formatting.R")
  
  # runname more information
  runname_upd <- paste0(runname, "_pred", predictor.type, "_out", outcome.type, "_Nrep", nsim)
  directory <- paste0("simulation_results_glmm/", runname_upd)
  
  # Create results directory
  dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  
  # set global seed
  set.seed(seed) 
  
  # Set up parallel backend
  future::plan(multisession, workers = parallelly::availableCores())
  
  # Run simulations in parallel
  parallel_results <- foreach(isim = 1:nsim,  .options.future = list(seed = TRUE), .verbose = TRUE) %dofuture% {
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
  
  # Save all raw results
  saveRDS(parallel_results, file = paste0(directory, "/parallel_results.rds"))
  
  # Apply result formatting **only once**
  formatted_results <- lapply(parallel_results, glmm_formating_results)
  
  # Compute mean results
  mean_results <- Reduce("+", formatted_results) / nsim
  
  # Compute variance (unbiased estimator)
  var_results <- Reduce("+", lapply(formatted_results, function(x) (x - mean_results)^2)) / (nsim - 1)
  sd_results <- sqrt(var_results)
  monte_carlo_se <- sd_results / sqrt(nsim)
  
  # record settings
  settings <- list(runname = runname, seed = seed, nsim = nsim, N_total = N_total, T_total = T_total, 
                   predictor.type = predictor.type, outcome.type = outcome.type,
                   sdX.within = sdX.within, sdX.between = sdX.between, 
                   g.00 = g.00, g.01 = g.01, sd.u0 = sd.u0, g.10 = g.10, 
                   sd.u1 = sd.u1, sd.e = sd.e)
  
  # combine all output
  output <- list(settings = settings, all_results = parallel_results, mean_results = mean_results, sd_results = sd_results, monte_carlo_se = monte_carlo_se)
  
  # save all output
  for (name in names(output)) {
    saveRDS(output[[name]], file = paste0(directory, "/", name, ".rds"))
  }
  
  # save warnings to text file
  sink(file = paste0(directory, "/warnings.txt"))
  summary(warnings())
  sink()
  
  return(output)
}

### RUNNING CONDITIONS ###

# stable specifications across simulations
# seed = 4243
# nsim = 1000
# N_total = 200
# T_total = 10

contxy_sim <- run_simulation(runname = "March6_checknewpc_noworker_try3", seed = 4243, nsim = 1000,
                             N_total = 200, T_total = 20, predictor.type = "continuous", outcome.type = "continuous",
                             sdX.within = 0.25, sdX.between = 0.5, g.00 = 0, g.01 = 1, sd.u0 = 0.7,
                             g.10 = 0.5, sd.u1 = 0, sd.e = 0.5)

round(contxy_sim$mean_results, 4)
contxy_sim$monte_carlo_se
summary(warnings())

# observations
# - once we increase T_total, the total effect is comprised more of the within-person effect, which explains
#   the similarity between the "uninterpretable blend" of raw X and the within-person effects.

binx_conty_sim1 <- run_simulation(runname = "postfix_reference", seed = 4243, nsim = 100,
                                 N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "continuous",
                                 sdX.within = NA, sdX.between = 0.5, g.00 = 0, g.01 = 0.5, sd.u0 = 0.7,
                                 g.10 = 0.5, sd.u1 = 0, sd.e = 0.5)

# round(binx_conty_sim$mean_results, 4)
# binx_conty_sim$monte_carlo_se
# summary(warnings())

binx_conty_sim1 <- run_simulation(runname = "postfix_g.01andg.10is2", seed = 4243, nsim = 1000,
                                  N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "continuous",
                                  sdX.within = NA, sdX.between = 0.5, g.00 = 0, g.01 = 2, sd.u0 = 0.7,
                                  g.10 = 2, sd.u1 = 0, sd.e = 0.5)

round(binx_conty_sim1$mean_results, 4)

binx_conty_sim1.1 <- run_simulation(runname = "postfix_g.01is1g.10is-1", seed = 4243, nsim = 1000,
                                  N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "continuous",
                                  sdX.within = NA, sdX.between = 0.5, g.00 = 0, g.01 = 1, sd.u0 = 0.7,
                                  g.10 = -1, sd.u1 = 0, sd.e = 0.5)

round(binx_conty_sim1.1$mean_results, 4)

binx_conty_sim2 <- run_simulation(runname = "postfix_highsdxbetw", seed = 4243, nsim = 100,
                                 N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "continuous",
                                 sdX.within = NA, sdX.between = 3, g.00 = 0, g.01 = 0.5, sd.u0 = 0.7,
                                 g.10 = 0.5, sd.u1 = 0, sd.e = 0.5)

# round(binx_conty_sim$mean_results, 4)
# binx_conty_sim$monte_carlo_se
# summary(warnings())

binx_conty_sim3 <- run_simulation(runname = "postfix_highg.01_lowT_realmean", seed = 4243, nsim = 1000,
                                 N_total = 200, T_total = 5, predictor.type = "binary", outcome.type = "continuous",
                                 sdX.within = NA, sdX.between = 0.5, g.00 = 0, g.01 = 2, sd.u0 = 0.7,
                                 g.10 = 0.5, sd.u1 = 0, sd.e = 0.5)

round(binx_conty_sim3$mean_results, 4)

binx_conty_sim3.3 <- run_simulation(runname = "postfix_highg.01_realmean2", seed = 4243, nsim = 1000,
                                  N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "continuous",
                                  sdX.within = NA, sdX.between = 0.5, g.00 = 0, g.01 = 2, sd.u0 = 0.7,
                                  g.10 = 0.5, sd.u1 = 0, sd.e = 0.5)

round(binx_conty_sim3.3$mean_results, 4)

binx_conty_sim3.1 <- run_simulation(runname = "postfix_highg.01_highT", seed = 4243, nsim = 1000,
                                    N_total = 200, T_total = 30, predictor.type = "binary", outcome.type = "continuous",
                                    sdX.within = NA, sdX.between = 0.5, g.00 = 0, g.01 = 2, sd.u0 = 0.7,
                                    g.10 = 0.5, sd.u1 = 0, sd.e = 0.5)


binx_conty_sim4 <- run_simulation(runname = "postfix_highg.01andg.10", seed = 4243, nsim = 100,
                                  N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "continuous",
                                  sdX.within = NA, sdX.between = 0.5, g.00 = 0, g.01 = 2, sd.u0 = 0.7,
                                  g.10 = 3, sd.u1 = 0, sd.e = 0.5)

# round(binx_conty_sim4$mean_results, 4)
# # binx_conty_sim$monte_carlo_se
# summary(warnings())

binx_conty_sim5 <- run_simulation(runname = "postfix_highg.01andg.10andsdXbet", seed = 4243, nsim = 100,
                                  N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "continuous",
                                  sdX.within = NA, sdX.between = 2, g.00 = 0, g.01 = -1, sd.u0 = 0.7,
                                  g.10 = 3, sd.u1 = 0, sd.e = 0.5)

# round(binx_conty_sim5$mean_results, 4)
# # binx_conty_sim$monte_carlo_se
# summary(warnings())

binx_conty_sim6 <- run_simulation(runname = "newpar_postfix_highg.01andg.10andsdXbet1", seed = 4243, nsim = 100,
                                  N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "continuous",
                                  sdX.within = NA, sdX.between = 1, g.00 = 0, g.01 = 2, sd.u0 = 0.7,
                                  g.10 = 3, sd.u1 = 0, sd.e = 0.5)

# round(binx_conty_sim6$mean_results, 4)
# # binx_conty_sim$monte_carlo_se
# summary(warnings())

# write the summary(warnings()) to file

contx_biny_sim <- run_simulation(runname = "March6", seed = 4243, nsim = 1000,
                                 N_total = 200, T_total = 20, predictor.type = "continuous", outcome.type = "binary",
                                 sdX.within = 0.25, sdX.between = 0.5, g.00 = 0, g.01 = 0.5, sd.u0 = 0.7,
                                 g.10 = 0.5, sd.u1 = 0, sd.e = NA)

round(contx_biny_sim$mean_results, 4)
contx_biny_sim$monte_carlo_se
summary(warnings())

binaryxy_sim1 <- run_simulation(runname = "lower-g.01-vals-fix_solved", seed = 4243, nsim = 1000,
                               N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "binary",
                               sdX.within = NA, sdX.between = 0.5, g.00 = 0, g.01 = 0.5, sd.u0 = 0.7,
                               g.10 = 0.5, sd.u1 = 0, sd.e = NA)

# round(binaryxy_sim$mean_results, 4)
# binaryxy_sim$monte_carlo_se
# summary(warnings())

binaryxy_sim2 <- run_simulation(runname = "higher-g.01-vals-fix_solved", seed = 4243, nsim = 1000,
                               N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "binary",
                               sdX.within = NA, sdX.between = 0.5, g.00 = 0, g.01 = 0.8, sd.u0 = 0.7,
                               g.10 = 0.5, sd.u1 = 0, sd.e = NA)

binaryxy_sim3 <- run_simulation(runname = "newpc", seed = 4243, nsim = 1000,
                                N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "binary",
                                sdX.within = NA, sdX.between = 1, g.00 = 0, g.01 = 0.8, sd.u0 = 0.7,
                                g.10 = 0.5, sd.u1 = 0, sd.e = NA)

binaryxy_sim3 <- run_simulation(runname = "newpc5", seed = 4243, nsim = 1000,
                                N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "binary",
                                sdX.within = NA, sdX.between = 1, g.00 = 0, g.01 = 0.8, sd.u0 = 0.7,
                                g.10 = 0.5, sd.u1 = 0, sd.e = NA)

binaryxy_sim3 <- run_simulation(runname = "newpc5_norandomeff", seed = 4243, nsim = 1000,
                                N_total = 200, T_total = 20, predictor.type = "binary", outcome.type = "binary",
                                sdX.within = NA, sdX.between = 1, g.00 = 0, g.01 = 0.8, sd.u0 = 0,
                                g.10 = 0.5, sd.u1 = 0, sd.e = NA)

# round(binaryxy_sim$mean_results, 4)
# binaryxy_sim$monte_carlo_se
# summary(warnings())

### Retrieve results

# output from GLMM with logit link is in log-odds, so to get the probability we need to transform it
# p = exp(logit) / (1 + exp(logit))

# transform every element in matrix mean_results
# apply(mean_results, c(1, 2), function(x) exp(x) / (1 + exp(x)))
