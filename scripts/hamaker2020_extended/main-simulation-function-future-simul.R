rm(list = ls()) # clear workspace

set.seed(123) # set global seed
runname <- "March24.7" # set a runname
dir.create(paste0("simulation_results_glmm/", runname), showWarnings = FALSE) # create a directory

# load libraries
library(lme4) # for generalized linear mixed models
library(geepack) # for generalized estimating equations (yields more stable results than gee package: https://www2.stat.duke.edu/~fl35/teaching/610-23F/docs/slides/6-1-GEE.pdf)
library(doFuture) # for parallelization (neatly handles RNG, errors, making objects available for workers)
library(parallelly) # for parallelization 
library(foreach) # for parallelization

# load helper functions
source("scripts/hamaker2020_extended/helper-functions/data-generation-centeredX.R")
source("scripts/hamaker2020_extended/helper-functions/model-fitting.R")
source("scripts/hamaker2020_extended/helper-functions/result-formatting.R")

# set the number of simulations
nsim <- 1000

# simulation for Research Report
design <- expand.grid(N_total = 200, T_total = c(10, 30), 
                      predictor.type = "continuous", outcome.type = "continuous",
                      sdX.within = sqrt(1), sdX.between = sqrt(4), 
                      g.00 = 0, g.01 = 2, sd.u0 = 1, g.10 = c(1, 2), 
                      sd.u1 = 0, sd.e = 1)

for (idesign in 1:nrow(design)) {
  
  # extract design parameters
  N_total <- design$N_total[idesign]
  T_total <- design$T_total[idesign]
  predictor.type <- design$predictor.type[idesign]
  outcome.type <- design$outcome.type[idesign]
  sdX.within <- design$sdX.within[idesign]
  sdX.between <- design$sdX.between[idesign]
  g.00 <- design$g.00[idesign]
  g.01 <- design$g.01[idesign]
  sd.u0 <- design$sd.u0[idesign]
  g.10 <- design$g.10[idesign]
  sd.u1 <- design$sd.u1[idesign]
  sd.e <- design$sd.e[idesign]
  
  # Set up parallel backend
  future::plan(multisession, workers = 16)
  
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
  
  saveRDS(parallel_results, file = paste0("simulation_results_glmm/", runname, "/", idesign, ".RDS"))

}

### collect results ---------------------------------------------------------

library(purrr) # for functional programming
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation

for (idesign in 1:nrow(design)) {
  
  parallel_results_setting <- readRDS(paste0("simulation_results_glmm/", runname, "/", idesign, ".RDS"))
  
  # unlist the lists inside the list
  df <- map_dfr(parallel_results_setting, function(rep) {
    map_dfr(rep, ~ as.data.frame(as.list(.x)), .id = "model")
  }, .id = "replication")
  
  # average across replications (not across models)
  df2 <- df %>%
    group_by(model) %>%
    summarise(across(c("X", "X.cent", "X.cluster.means"), mean, na.rm = TRUE))
  
  # change into wide format
  df3 <- pivot_wider(df2, names_from = model, values_from = c("X", "X.cent", "X.cluster.means"))
  
}
  
  
  
  # unlist the lists inside the list
  results <- lapply(parallel_results, unlist)
  
  # note model names
  model_names <- c("l1", "l2", "l3a", "l4", 
                  "g.independence1", "g.exchangeable1", "g.ar11", 
                  "g.independence2", "g.exchangeable2", "g.ar12", 
                  "g.independence3a", "g.exchangeable3a", "g.ar13a", 
                  "g.independence4", "g.exchangeable4", "g.ar14")
  
  # for every model name, extract the results
  for (model_name in model_names) {
    results_list[grep("solution_lmm", names(results_list))]
    results <- lapply(parallel_results, function(x) x[[model_name]])
    saveRDS(results, file = paste0("simulation_results_glmm/", runname, "/", idesign, "_", model_name, ".RDS"))
  }
  
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


