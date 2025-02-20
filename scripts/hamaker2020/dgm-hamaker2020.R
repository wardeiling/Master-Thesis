######################################
#### MODULARIZED SIMULATION SCRIPT ###
######################################

library(lme4)
library(gee)
library(geepack)
library(tidyverse)
library(foreach)
library(doParallel)


# MODEL FITTING FUNCTION
fit_models <- function(data, outcome.type) {
  models <- list()
  
  if (outcome.type == "continuous") {
    family_arg <- gaussian()
  } else if (outcome.type == "binary") {
    family_arg <- binomial(link = "logit")
  }
  
  models$l1 <- summary(lmer(Y ~ X + (1 | Cluster), data = data, REML = FALSE))$coefficients[, 1]
  models$g.indep1 <- coef(gee(Y ~ X, id = Cluster, data = data, corstr = "independence", family = family_arg))
  models$g.exch1 <- coef(gee(Y ~ X, id = Cluster, data = data, corstr = "exchangeable", family = family_arg))
  models$g.ar1 <- coef(gee(Y ~ X, id = Cluster, data = data, corstr = "AR-M", Mv = 1, family = family_arg))
  
  models$l2 <- summary(lmer(Y ~ X.cent + (1 | Cluster), data = data, REML = FALSE))$coefficients[, 1]
  models$g.indep2 <- coef(gee(Y ~ X.cent, id = Cluster, data = data, corstr = "independence", family = family_arg))
  models$g.exch2 <- coef(gee(Y ~ X.cent, id = Cluster, data = data, corstr = "exchangeable", family = family_arg))
  models$g.ar2 <- coef(gee(Y ~ X.cent, id = Cluster, data = data, corstr = "AR-M", Mv = 1, family = family_arg))
  
  models$l3a <- summary(lmer(Y ~ X.cent + cluster.means + (1 | Cluster), data = data, REML = FALSE))$coefficients[, 1]
  models$g.indep3a <- coef(gee(Y ~ X.cent + cluster.means, id = Cluster, data = data, corstr = "independence", family = family_arg))
  models$g.exch3a <- coef(gee(Y ~ X.cent + cluster.means, id = Cluster, data = data, corstr = "exchangeable", family = family_arg))
  models$g.ar3a <- coef(gee(Y ~ X.cent + cluster.means, id = Cluster, data = data, corstr = "AR-M", Mv = 1, family = family_arg))
  
  models$l4 <- summary(lmer(Y ~ X + cluster.means + (1 | Cluster), data = data, REML = FALSE))$coefficients[, 1]
  models$g.indep4 <- coef(gee(Y ~ X + cluster.means, id = Cluster, data = data, corstr = "independence", family = family_arg))
  models$g.exch4 <- coef(gee(Y ~ X + cluster.means, id = Cluster, data = data, corstr = "exchangeable", family = family_arg))
  models$g.ar4 <- coef(gee(Y ~ X + cluster.means, id = Cluster, data = data, corstr = "AR-M", Mv = 1, family = family_arg))
  
  return(models)
}

# RESULT FORMATTING FUNCTION
format_results <- function(models) {
  table <- matrix(NA, nrow = 16, ncol = 3, 
                  dimnames = list(c("L1", "GEE.indep1", "GEE.exch1", "GEE.ar1",
                                    "L2", "GEE.indep2", "GEE.exch2", "GEE.ar2",
                                    "L3a", "GEE.indep3a", "GEE.exch3a", "GEE.ar3a",
                                    "L4", "GEE.indep4", "GEE.exch4", "GEE.ar4"),
                                  c("X", "X.cent", "cluster.means")))
  
  for (i in 1:16) {
    names <- names(models[[i]])
    if ("X" %in% names) table[i, 1] <- models[[i]]["X"]
    if ("X.cent" %in% names) table[i, 2] <- models[[i]]["X.cent"]
    if ("cluster.means" %in% names) table[i, 3] <- models[[i]]["cluster.means"]
  }
  
  return(table)
}

# MAIN WRAPPER FUNCTION
run_simulation <- function(N = 5000, n = 20, predictor.type = "continuous", outcome.type = "continuous", seed = 3859, 
                           sdX.within = sqrt(1), sdX.between = sqrt(4), g.00 = 0, g.01 = 2, sd.u0 = 1, 
                           g.10 = 1, sd.u1 = 0, sd.e = 1) {
  
  data <- generate_data(N, n, predictor.type, outcome.type, seed, sdX.within, sdX.between, g.00, g.01, sd.u0, g.10, sd.u1, sd.e)
  data <- prepare_data(data)
  models <- fit_models(data, outcome.type)
  result_table <- format_results(models)
  
  return(result_table)
}

# PARALLELIZATION SETUP
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# RUN MULTIPLE SIMULATIONS IN PARALLEL
results <- foreach(i = 1:10, .combine = rbind, .packages = c("lme4", "gee", "geepack", "tidyverse")) %dopar% {
  run_simulation(N = 1000, n = 10, seed = 3859 + i)
}

stopCluster(cl)
print(results)
