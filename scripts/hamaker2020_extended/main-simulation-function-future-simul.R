######################################
#### MODULARIZED SIMULATION SCRIPT ###
######################################

# Author: Ward B Eiling
# Last edited: 26-03-2025

### Simulation set-up -------------------------------------------------------

rm(list = ls()) # clear workspace

seed <- 6384
set.seed(seed) # set seed for reproducibility
runname <- "March27_design5b_ludtkesbias_contextual_trueclustermeans" # set a runname
parametrization <- "mundlak" # set the parametrization (mundlak or centeredX)
dir.create(paste0("simulation_results_glmm/", runname), showWarnings = FALSE) # create a directory

# load libraries
library(lme4) # for generalized linear mixed models
library(geepack) # for generalized estimating equations (yields more stable results than gee package: https://www2.stat.duke.edu/~fl35/teaching/610-23F/docs/slides/6-1-GEE.pdf)
library(doFuture) # for parallelization (neatly handles RNG, errors, making objects available for workers)
library(parallelly) # for parallelization 
library(foreach) # for parallelization

# load helper functions
if (parametrization == "mundlak") {
  source("scripts/hamaker2020_extended/helper-functions/data-generation-mundlak.R")
} else if (parametrization == "centeredX") {
  source("scripts/hamaker2020_extended/helper-functions/data-generation-centeredX.R")
}
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

# design 1 (test influence of general parameters on bias)
# design <- expand.grid(N_total = 200, T_total = 20, 
#                       predictor.type = "binary", outcome.type = "continuous",
#                       sdX.within = NA, sdX.between = c(0, 1), 
#                       g.00 = 0, g.01 = c(-1, 0, 1), sd.u0 = c(0, 1), g.10 = c(0.5, 1.5, 3), 
#                       sd.u1 = 0, sd.e = 1, true_cluster_means = FALSE)

# design 2 (test influence of true cluster means on bias)
# design <- expand.grid(N_total = 200, T_total = 20,
#                       predictor.type = "binary", outcome.type = "continuous",
#                       sdX.within = NA, sdX.between = c(0, 1),
#                       g.00 = 0, g.01 = c(-1, 0, 1), sd.u0 = 1, g.10 = c(0.5, 1.5, 3),
#                       sd.u1 = 0, sd.e = 1, true_cluster_means = c(FALSE, TRUE))

# # design 2b (true cluster mean)
# design <- expand.grid(N_total = 200, T_total = 20,
#                       predictor.type = "binary", outcome.type = "continuous",
#                       sdX.within = NA, sdX.between = c(0, 1),
#                       g.00 = 0, g.01 = c(-1, 0, 1), sd.u0 = 1, g.10 = c(0.5, 1.5, 3),
#                       sd.u1 = 0, sd.e = 1, true_cluster_means = TRUE)

# # design 3 (test influence T, N, sdX.between on bias)
# design <- expand.grid(N_total = c(100, 200), T_total = c(5, 20),
#                       predictor.type = "binary", outcome.type = "continuous",
#                       sdX.within = NA, sdX.between = c(0, 1, 3),
#                       g.00 = 0, g.01 = 1, sd.u0 = 1, g.10 = c(0.8, 2),
#                       sd.u1 = 0, sd.e = 1, true_cluster_means = c(FALSE, TRUE))

# # design 3b (true cluster mean)
# design <- expand.grid(N_total = c(100, 200), T_total = c(5, 20),
#                       predictor.type = "binary", outcome.type = "continuous",
#                       sdX.within = NA, sdX.between = c(0, 1, 3),
#                       g.00 = 0, g.01 = 1, sd.u0 = 1, g.10 = c(0.8, 2),
#                       sd.u1 = 0, sd.e = 1, true_cluster_means = TRUE)

# design 4 (test influence random intercept and slope on bias)
# design <- expand.grid(N_total = 200, T_total = 20, 
#                       predictor.type = "binary", outcome.type = "continuous",
#                       sdX.within = NA, sdX.between = 1, 
#                       g.00 = c(0, 1), g.01 = c(0, 1.5), sd.u0 = c(0, 1), g.10 = 0.8, 
#                       sd.u1 = c(0, 1), sd.e = 1, true_cluster_means = c(FALSE, TRUE))

# design 4b (true cluster mean)
# design <- expand.grid(N_total = 200, T_total = 20,
#                       predictor.type = "binary", outcome.type = "continuous",
#                       sdX.within = NA, sdX.between = 1,
#                       g.00 = c(0, 1), g.01 = c(0, 1.5), sd.u0 = c(0, 1), g.10 = 0.8,
#                       sd.u1 = c(0, 1), sd.e = 1, true_cluster_means = TRUE)

# design 5a (test Ludtke's bias for all predictor and outcome types)
# design <- expand.grid(N_total = c(100, 200), T_total = c(5, 20),
#                       predictor.type = c("binary", "continuous"),
#                       outcome.type = c("binary", "continuous"),
#                       sdX.within = 1, sdX.between = c(0, 1, 3),
#                       g.00 = 0, g.01 = c(0, 1), g.10 = 1.5,
#                       sd.u0 = c(0, 1), sd.u1 = 0, sd.e = 1,
#                       true_cluster_means = FALSE)
# # remove scenarios with sdX.between == 0 and g.01 != 0
# design <- design[!(design$sdX.between == 0 & design$g.01 != 0),]

# design 5b
design <- expand.grid(N_total = c(100, 200), T_total = c(5, 20),
                      predictor.type = c("binary", "continuous"),
                      outcome.type = c("binary", "continuous"),
                      sdX.within = 1, sdX.between = c(0, 1, 3),
                      g.00 = 0, g.01 = c(0, 1), g.10 = 1.5,
                      sd.u0 = c(0, 1), sd.u1 = 0, sd.e = 1,
                      true_cluster_means = TRUE)
# remove scenarios with sdX.between == 0 and g.01 != 0
design <- design[!(design$sdX.between == 0 & design$g.01 != 0),]

# save the empty design and settings to the directory
settings <- list(nsim = nsim, seed = seed, runname = runname, parametrization = parametrization, design = design)
saveRDS(settings, paste0("simulation_results_glmm/", runname, "/settings.RDS"))

### Simulation ---------------------------------------------------------------

# loop over the design
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
  true_cluster_means <- design$true_cluster_means[idesign]
  
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
                                 g.10 = g.10, sd.u1 = sd.u1, sd.e = sd.e, 
                                 true_cluster_means = true_cluster_means)
    
    # Fit models
    models <- glmm_model_fitting(data, outcome.type = outcome.type)
    return(models)
  }
  saveRDS(parallel_results, file = paste0("simulation_results_glmm/", runname, "/", idesign, ".RDS"))
}

### collect results ---------------------------------------------------------

if(FALSE) {
  
  # design 1
  runname <- "March27_design5_ludtkesbias_contextual_estclustermeans"
  design <- readRDS(paste0("simulation_results_glmm/", runname, "/settings.RDS"))$design
  parametrization <- readRDS(paste0("simulation_results_glmm/", runname, "/settings.RDS"))$parametrization
  
  # design 2
  runname <- "March27_design5b_ludtkesbias_contextual_trueclustermeans"
  design <- readRDS(paste0("simulation_results_glmm/", runname, "/settings.RDS"))$design
  parametrization <- readRDS(paste0("simulation_results_glmm/", runname, "/settings.RDS"))$parametrization
}

# load packages
library(purrr) # for functional programming
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation

# initialize design and add columns for parameter values
design_all <- design %>%
  # add columns for MLM estimation
  mutate(l1_X = NA, l2_X.cent = NA, l3a_X.cent = NA, l3a_X.cluster.means = NA, l4_X = NA, l4_X.cluster.means = NA,
         # add columns for MLM estimation error / bias
         l2_g.10_bias = NA, l3a_g.10_bias = NA, l3a_g.01_bias = NA, l4_g.10_bias = NA, l4_g.01_bias = NA,
         # add columns for GEE estimation
         g.independence1_X = NA, g.exchangeable1_X = NA, g.ar11_X = NA,
         g.independence2_X.cent = NA, g.exchangeable2_X.cent = NA, g.ar12_X.cent = NA,
         g.independence3_X.cent = NA, g.independence3_X.cluster.means = NA,
         g.exchangeable3_X.cent = NA, g.exchangeable3_X.cluster.means = NA,
         g.ar13_X.cent = NA, g.ar13_X.cluster.means = NA,
         g.independence4_X = NA, g.independence4_X.cluster.means = NA,
         g.exchangeable4_X = NA, g.exchangeable4_X.cluster.means = NA,
         g.ar14_X = NA, g.ar14_X.cluster.means = NA,
         # add columns for GEE estimation error / bias
         g.independence2_g.10_bias = NA, g.exchangeable2_g.10_bias = NA, g.ar12_g.10_bias = NA,
         g.independence3_g.10_bias = NA, g.independence3_g.01_bias = NA,
         g.exchangeable3_g.10_bias = NA, g.exchangeable3_g.01_bias = NA,
         g.ar13_g.10_bias = NA, g.ar13_g.01_bias = NA,
         g.independence4_g.10_bias = NA, g.independence4_g.01_bias = NA,
         g.exchangeable4_g.10_bias = NA, g.exchangeable4_g.01_bias = NA,
         g.ar14_g.10_bias = NA, g.ar14_g.01_bias = NA)

# loop over the design
for (idesign in 1:nrow(design_all)) {
  
  # extract parameter values from design
  g.00 <- design_all$g.00[idesign]
  g.01 <- design_all$g.01[idesign]
  g.10 <- design_all$g.10[idesign]
  
  # determine the beta values based on the parametrization
  if (parametrization == "centeredX") {
    
    # if centeredX, g.01 is the between-person effect
    beta_between <- g.01
    beta_within <- g.10
    beta_contextual <- beta_between - beta_within
    
  } else if (parametrization == "mundlak") {
    
    # if mundlak, g.01 is the contextual effect
    beta_contextual <- g.01
    beta_within <- g.10
    beta_between <- beta_within + beta_contextual
  }
  
  # read in the results
  parallel_results_setting <- readRDS(paste0("simulation_results_glmm/", runname, "/", idesign, ".RDS"))
  
  # unlist the lists inside the list
  df <- map_dfr(parallel_results_setting, function(rep) {
    map_dfr(rep, ~ as.data.frame(as.list(.x)), .id = "model")
  }, .id = "replication")
  
  # rename columns (to address problem where X.cluster.means is not estimated)
  colnames(df) <- c("replication", "model", "X.Intercept.", "X", "X.cent", "X.cluster.means")
  
  # average across replications
  df_processed <- df %>%
    group_by(model) %>%
    summarise(across(c("X", "X.cent", "X.cluster.means"), mean, na.rm = TRUE)) %>%
    # for now only select the mlm models
    # filter(model %in% c("l1", "l2", "l3a", "l4")) %>%
    # change into wide format
    pivot_wider(names_from = model, values_from = c("X", "X.cent", "X.cluster.means"), names_glue = "{model}_{.value}") %>%
    # only select possible columns
    select("l1_X", "l2_X.cent", "l3a_X.cent", "l3a_X.cluster.means", "l4_X", "l4_X.cluster.means",
           # also for GEE
           "g.independence1_X", "g.exchangeable1_X", "g.ar11_X",
           "g.independence2_X.cent", "g.exchangeable2_X.cent", "g.ar12_X.cent",
           "g.independence3_X.cent", "g.independence3_X.cluster.means",
           "g.exchangeable3_X.cent", "g.exchangeable3_X.cluster.means",
           "g.ar13_X.cent", "g.ar13_X.cluster.means",
           "g.independence4_X", "g.independence4_X.cluster.means",
           "g.exchangeable4_X", "g.exchangeable4_X.cluster.means",
           "g.ar14_X", "g.ar14_X.cluster.means"
           ) %>%
    # calculate estimation error (l2, l3a, l4)
    # l3a_X.cluster.means is the between-person effect beta(between)
    # l4_X.cluster.means is the contextual effect = beta(between) - beta(within)
    mutate(l2_g.10_bias = l2_X.cent - beta_within,
           l3a_g.10_bias = l3a_X.cent - beta_within,
           l3a_g.01_bias = l3a_X.cluster.means - beta_between,
           l4_g.10_bias = l4_X - beta_within,
           l4_g.01_bias = l4_X.cluster.means - beta_contextual,
           # also for GEE
           g.independence2_g.10_bias = g.independence2_X.cent - beta_within,
           g.exchangeable2_g.10_bias = g.exchangeable2_X.cent - beta_within,
           g.ar12_g.10_bias = g.ar12_X.cent - beta_within,
           g.independence3_g.10_bias = g.independence3_X.cent - beta_within,
           g.independence3_g.01_bias = g.independence3_X.cluster.means - beta_between,
           g.exchangeable3_g.10_bias = g.exchangeable3_X.cent - beta_within,
           g.exchangeable3_g.01_bias = g.exchangeable3_X.cluster.means - beta_between,
           g.ar13_g.10_bias = g.ar13_X.cent - beta_within,
           g.ar13_g.01_bias = g.ar13_X.cluster.means - beta_between,
           g.independence4_g.10_bias = g.independence4_X - beta_within,
           g.independence4_g.01_bias = g.independence4_X.cluster.means - beta_contextual,
           g.exchangeable4_g.10_bias = g.exchangeable4_X - beta_within,
           g.exchangeable4_g.01_bias = g.exchangeable4_X.cluster.means - beta_contextual,
           g.ar14_g.10_bias = g.ar14_X - beta_within,
           g.ar14_g.01_bias = g.ar14_X.cluster.means - beta_contextual
           ) 
  
  # add everything to design, making sure that all columns align
  design_all[idesign, colnames(df_processed)] <- df_processed
  
}

# reorder the columns (first mlm, then GEE exchangeable, GEE AR(1), GEE independence)
design_all_reordered <- design_all %>%
  select(N_total, T_total, predictor.type, outcome.type, sdX.within, sdX.between,
         g.00, g.01, g.10, sd.u0, sd.u1, sd.e, true_cluster_means,
         # first absolute value columns (ordered according to centering method)
         l1_X, g.exchangeable1_X, g.ar11_X, g.independence1_X,
         l2_X.cent, g.exchangeable2_X.cent, g.ar12_X.cent,  g.independence2_X.cent,
         l3a_X.cent, g.exchangeable3_X.cent, g.ar13_X.cent, g.independence3_X.cent, 
         l3a_X.cluster.means, g.exchangeable3_X.cluster.means, g.ar13_X.cluster.means, g.independence3_X.cluster.means, 
         l4_X, g.exchangeable4_X, g.ar14_X, g.independence4_X,
         l4_X.cluster.means, g.exchangeable4_X.cluster.means, g.ar14_X.cluster.means, g.independence4_X.cluster.means,
         # then bias columns (ordered according to model type)
         l2_g.10_bias, l3a_g.10_bias, l3a_g.01_bias, l4_g.10_bias, l4_g.01_bias,
         g.exchangeable2_g.10_bias, g.exchangeable3_g.10_bias, g.exchangeable3_g.01_bias, g.exchangeable4_g.10_bias, g.exchangeable4_g.01_bias,
         g.ar12_g.10_bias, g.ar13_g.10_bias, g.ar13_g.01_bias, g.ar14_g.10_bias, g.ar14_g.01_bias,
         g.independence2_g.10_bias, g.independence3_g.10_bias, g.independence3_g.01_bias, g.independence4_g.10_bias, g.independence4_g.01_bias
         )

# save design
saveRDS(design_all_reordered, paste0("simulation_results_glmm/", runname, "/summary-results-all.RDS"))

# create a rounded version of the design
design_all_rounded <- design_all_reordered %>%
  # round numbers to 3 decimals: everything that startswith "l" or "g."
  mutate(across(starts_with("l"), ~ round(., 3)),
         across(starts_with("g."), ~ round(., 3)))

# save design
saveRDS(design_all_rounded, paste0("simulation_results_glmm/", runname, "/summary-results-all-rounded.RDS"))

# remove all bias columns
design_absolute <- design_all_rounded %>%
  select(-contains("g.10"), -contains("g.01"))

# save design
saveRDS(design_absolute, paste0("simulation_results_glmm/", runname, "/summary-results-absolute.RDS"))
write.csv(design_absolute, paste0("simulation_results_glmm/", runname, "/summary-results-absolute.csv"), row.names = FALSE)

# remove absolute value columns
design_bias <- design_all_rounded %>%
  select(-contains("X", ignore.case = FALSE))

# save design
saveRDS(design_bias, paste0("simulation_results_glmm/", runname, "/summary-results-bias.RDS"))
write.csv(design_bias, paste0("simulation_results_glmm/", runname, "/summary-results-bias.csv"), row.names = FALSE)

# create separate versions containing the within-person and contextual effect
design_bias_g01 <- design_bias %>%
  select(-contains("g.10"))
design_bias_g10 <- design_bias %>%
  select(-contains("g.01"))

# save design
saveRDS(design_bias_g01, paste0("simulation_results_glmm/", runname, "/summary-results-bias-g01.RDS"))
write.csv(design_bias_g01, paste0("simulation_results_glmm/", runname, "/summary-results-bias-g01.csv"), row.names = FALSE)
saveRDS(design_bias_g10, paste0("simulation_results_glmm/", runname, "/summary-results-bias-g10.RDS"))
write.csv(design_bias_g10, paste0("simulation_results_glmm/", runname, "/summary-results-bias-g10.csv"), row.names = FALSE)

### Optional: Retrieve output

if(FALSE) {
  # design 2
  runname <- "March26_design2_maineffects_contextual_bothclustermeans"
  design_bias_temp <- readRDS(paste0("simulation_results_glmm/", runname, "/summary-results-bias.RDS")) %>%
    # round the last six columns to 3 decimals
    mutate(across(starts_with("l"), ~ round(., 3)))
  design_bias_temp_sel <- design_bias_temp[,6:ncol(design_bias_temp)]
  
  # design 5
  runname <- "March27_design5_ludtkesbias_contextual_estclustermeans"
  design_bias_temp <- readRDS(paste0("simulation_results_glmm/", runname, "/summary-results-bias.RDS")) %>%
    # round the last six columns to 3 decimals
    mutate(across(starts_with("l"), ~ round(., 3)))
  design_bias_temp_sel <- design_bias_temp[,6:ncol(design_bias_temp)]
  
  design_bias_temp_sel <- design_bias_temp %>%
    filter(sdX.between == 1 & g.01 == 1 & sd.u0 == 1 & g.10 == 1.5 & sd.u1 == 0 & sd.e == 1) %>%
    # create column that combines the two factor variables predictor.type and outcome.type
    mutate(variabel.types = paste("pred", predictor.type, "out", outcome.type, sep = "_")) 
  
  # create line plots showing bias of g.01 (for model L3a and L4) against number of timepoints (T_total) with lines 
  # for each different model and in a 2x2 grid for predictor and outcome type
  library(ggplot2)
  
  # step 1: change the model type into long-format
  design_bias_temp_sel <- design_bias_temp_sel %>%
    pivot_longer(cols = starts_with("l"), names_to = "model", values_to = "bias") %>%
    mutate(model = factor(model, levels = c("l2_g.10_bias", "l3a_g.10_bias", "l3a_g.01_bias", "l4_g.10_bias", "l4_g.01_bias"))) 
    # only select g.01 from l3a and L4
    # filter(model %in% c("l3a_g.01_bias", "l4_g.01_bias"))
  
  # step 2: create the plot
  ggplot(design_bias_temp_sel, aes(x = T_total, y = bias, color = model)) +
    geom_line() +
    facet_grid(predictor.type ~ outcome.type) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(x = "Number of timepoints (T_total)", y = "Estimation error")
  
}
