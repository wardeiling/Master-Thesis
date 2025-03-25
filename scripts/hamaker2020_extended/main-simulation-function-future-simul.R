rm(list = ls()) # clear workspace

seed <- 6384
set.seed(seed)
runname <- "March25_design1_maineffects_contextual" # set a runname
dir.create(paste0("simulation_results_glmm/", runname), showWarnings = FALSE) # create a directory

# load libraries
library(lme4) # for generalized linear mixed models
library(geepack) # for generalized estimating equations (yields more stable results than gee package: https://www2.stat.duke.edu/~fl35/teaching/610-23F/docs/slides/6-1-GEE.pdf)
library(doFuture) # for parallelization (neatly handles RNG, errors, making objects available for workers)
library(parallelly) # for parallelization 
library(foreach) # for parallelization

# load helper functions
source("scripts/hamaker2020_extended/helper-functions/data-generation-mundlak.R")
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

# save the empty design and settings to the directory
settings <- list(nsim = nsim, seed = seed, runname = runname, design = design)
saveRDS(settings, paste0("simulation_results_glmm/", runname, "/settings.RDS"))

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
  future::plan(multisession, workers = 8)
  
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

### Version with Absolute Values
design_abs <- design
design_abs$l1_X <- design_abs$l2_X.cent <- design_abs$l3a_X.cent <- design_abs$l3a_X.cluster.means <- design_abs$l4_X <- design_abs$l4_X.cluster.means

for (idesign in 1:nrow(design_abs)) {
  
  # read in the results
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
  design_abs[idesign, c("l1_X", "l2_X.cent", "l3a_X.cent", "l3a_X.cluster.means", "l4_X", "l4_X.cluster.means")] <- df_processed
  
}

# save design
saveRDS(design_abs, paste0("simulation_results_glmm/", runname, "/summary-results-absolute.RDS"))

# selected results (from column 6 onwards)
design_abs_selected <- design_abs[,6:ncol(design_abs)]
round(design_abs_selected, 3)

### Version with Bias in Estimates

design_bias <- design
design_bias$l1_X <- design_bias$l2_X.cent <- design_bias$l3a_X.cent <- design_bias$l3a_X.cluster.means <- design_bias$l4_X <- design_bias$l4_X.cluster.means

for (idesign in 1:nrow(design_bias)) {
  
  # extract parameter values form design
  g.00 <- design_bias$g.00[idesign]
  g.01 <- design_bias$g.01[idesign]
  g.10 <- design_bias$g.10[idesign]
  
  # read in the results
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
    select("l1_X", "l2_X.cent", "l3a_X.cent", "l3a_X.cluster.means", "l4_X", "l4_X.cluster.means") %>%
    # calculate bias (l2, l3a, l4)
    mutate(l2_g.10_bias = l2_X.cent - g.10,
           l3a_g.10_bias = l3a_X.cent - g.10,
           l3a_g.01_bias = l3a_X.cluster.means - g.01,
           l4_g.10_bias = l4_X - g.10,
           l4_g.01_bias = l4_X.cluster.means - g.01
           ) %>%
    # select only bias columns
    select("l1_X", "l2_g.10_bias", "l3a_g.10_bias", "l3a_g.01_bias", "l4_g.10_bias", "l4_g.01_bias")
  
  # add to design
  design_bias[idesign, c("l1_X", "l2_g.10_bias", "l3a_g.10_bias", "l3a_g.01_bias", "l4_g.10_bias", "l4_g.01_bias")] <- df_processed
  
}

# save design
saveRDS(design_bias, paste0("simulation_results_glmm/", runname, "/summary-results-bias.RDS"))

# selected results (from column 6 onwards)
design_bias_selected <- design_bias[,6:ncol(design_bias)]
round(design_bias_selected, 3)