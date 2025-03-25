rm(list = ls()) # clear workspace

seed = 4243 # set global seed
runname <- "March25_design1_maineffects" # set a runname
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
nsim <- 100

# comprehensive design
# design <- expand.grid(N_total = c(100, 200), T_total = c(5, 20), 
#                       predictor.type = "binary", outcome.type = "continuous",
#                       sdX.within = NA, sdX.between = c(0, 0.5, 1.5), 
#                       g.00 = 0, g.01 = c(-1, 0, 1), sd.u0 = c(0, 0.5, 1.5), g.10 = c(0.5, 1.5, 3), 
#                       sd.u1 = c(0, 0.5, 1.5), sd.e = c(0.5, 1.5))

# subdesign 1
design <- expand.grid(N_total = 200, T_total = 20, 
                      predictor.type = "binary", outcome.type = "continuous",
                      sdX.within = NA, sdX.between = c(0, 1), 
                      g.00 = 0, g.01 = c(-1, 0, 1), sd.u0 = c(0, 1), g.10 = c(0.5, 1.5, 3), 
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
  
  saveRDS(parallel_results, file = paste0("simulation_results_glmm/", runname, "/", idesign, ".RDS"))

}

### collect results ---------------------------------------------------------

library(purrr) # for functional programming
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation

design$l1_X <- design$l2_X.cent <- design$l3a_X.cent <- design$l3a_X.cluster.means <- design$l4_X <- design$l4_X.cluster.means

for (idesign in 1:nrow(design)) {
  
  parallel_results_setting <- readRDS(paste0("simulation_results_glmm/", runname, "/", idesign, ".RDS"))
  
  # unlist the lists inside the list
  df <- map_dfr(parallel_results_setting, function(rep) {
    map_dfr(rep, ~ as.data.frame(as.list(.x)), .id = "model")
  }, .id = "replication")
  
  # average across replications
  df_processed <- df %>%
    group_by(model) %>%
    summarise(across(c("X", "X.cent", "X.cluster.means"), mean, na.rm = TRUE)) %>%
    # for now only select the mlm models
    filter(model %in% c("l1", "l2", "l3a", "l4")) %>%
    # change into wide format
    pivot_wider(names_from = model, values_from = c("X", "X.cent", "X.cluster.means"), names_glue = "{model}_{.value}") %>%
    # only select possible columns
    select("l1_X", "l2_X.cent", "l3a_X.cent", "l3a_X.cluster.means", "l4_X", "l4_X.cluster.means")
  
  # add to design
  design[idesign, c("l1_X", "l2_X.cent", "l3a_X.cent", "l3a_X.cluster.means", "l4_X", "l4_X.cluster.means")] <- df_processed
  
}

# save design
saveRDS(design, paste0("simulation_results_glmm/", runname, "/summary-results.RDS"))
