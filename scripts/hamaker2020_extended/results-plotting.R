# load packages
library(purrr) # for functional programming
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation

runname <- "April10_fullsimulation"

# retrieve the design and parameterization
design <- readRDS(paste0("simulation_results_glmm/", runname, "/settings.RDS"))$design
parametrization <- readRDS(paste0("simulation_results_glmm/", runname, "/settings.RDS"))$parametrization
nsim <- readRDS(paste0("simulation_results_glmm/", runname, "/settings.RDS"))$nsim

# initialize design and add columns for parameter values
design_all <- design %>%
  # add columns for MLM estimation
  mutate(# add columns for MLM estimation error / bias
         l2_g.10_bias = NA, l3a_g.10_bias = NA, l3a_g.01_bias = NA, l4_g.10_bias = NA, l4_g.01_bias = NA,
         # add columns for GEE estimation error / bias
         g.independence2_g.10_bias = NA, g.exchangeable2_g.10_bias = NA, g.ar12_g.10_bias = NA,
         g.independence3_g.10_bias = NA, g.independence3_g.01_bias = NA,
         g.exchangeable3_g.10_bias = NA, g.exchangeable3_g.01_bias = NA,
         g.ar13_g.10_bias = NA, g.ar13_g.01_bias = NA,
         g.independence4_g.10_bias = NA, g.independence4_g.01_bias = NA,
         g.exchangeable4_g.10_bias = NA, g.exchangeable4_g.01_bias = NA,
         g.ar14_g.10_bias = NA, g.ar14_g.01_bias = NA)

# loop over the design
for (idesign in 1:2) {
  
  ### Extract parameter values from design
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
  
  # # tranpose structure of lists
  # parallel_results_setting_transp <- transpose(parallel_results_setting)
  
  # unlist the lists inside the list
  df <- map_dfr(parallel_results_setting, function(rep) {
    map_dfr(rep, ~ as.data.frame(as.list(.x)), .id = "model")
  }, .id = "replication")
  
  df_wide <- df %>%
    select(-X.Intercept.) %>%                            # (1) Remove intercept
    pivot_wider(id_cols = replication, names_from = model, values_from = c("X", "X.cent", "X.cluster.means"), names_glue = "{model}_{.value}") %>%
    select(where(~ !all(is.na(.))))
  
}






