# simulation study for the random effects model paper (Statistical Science)

rm(list = ls())
set.seed(123)

library(rootSolve)
library(geepack)
library(lme4)
library(mvtnorm)
library(foreach)
library(doParallel)  # Replace doMC with doParallel

source("simulation/scripts/functions/generative_models_function.R")

# simulation: using lmer() package -----------------------------------------

# Set up parallel backend to use many processors
cores <- detectCores()
cl <- makeCluster(cores - 1)  # Create a cluster
registerDoParallel(cl)        # Register the parallel backend

set.seed(120)

nsim <- 20

# full simulation in Appendix B
# design <- expand.grid(sample_size = c(30, 50, 100, 200), total_T = c(10, 20, 30), dgm_type = 1:3)

# simulation in Section 4
design <- expand.grid(sample_size = c(30, 100, 200), total_T = c(10, 30), dgm_type = 1:3)

design <- design[order(design$dgm_type), ]

for (idesign in 1:nrow(design)) {
  sample_size <- design$sample_size[idesign]
  total_T <- design$total_T[idesign]
  dgm_type <- design$dgm_type[idesign]
  
  # Create template output structure for simulated trials with estimation error
  dta <- dgm_with_treatment(sample_size, total_T, dgm_type = dgm_type)
  if (dgm_type %in% c(1,3,4)) {
    fit <- lme4::lmer(Y ~ Z * X + (1 + X| userid), data = dta)
  } else if (dgm_type == 2) {
    fit <- lme4::lmer(Y ~ Z * X + (Z * X| userid), data = dta)
  }
  coef <- summary(fit)$coefficients
  varcor <- summary(fit)$varcor
  for (irow in 1:nrow(coef)) {
    for (icol in 1:ncol(coef)) {
      coef[irow, icol] <- NA
    }
  }
  for (irow in 1:nrow(varcor$userid)) {
    for (icol in 1:ncol(varcor$userid)) {
      varcor$userid[irow, icol] <- NA
    }
  }
  attr(varcor, "sc") <- NA
  coef_NA_fill <- coef
  varcor_NA_fill <- varcor
  
  # Start parallel jobs
  writeLines(c(""), "log.txt")
  sink("log.txt", append = FALSE)
  set.seed(123)
  
  result <- foreach(isim = 1:nsim, .combine = "c") %dopar% {  # Use %dopar% for parallel execution
    
    library(lme4) # Load the library in the worker process
    
    if (isim %% 10 == 0) {
      cat(paste("Starting iteration", isim, "\n"))
    }
    dta <- dgm_with_treatment(sample_size, total_T, dgm_type = dgm_type)
    
    # solution <- tryCatch(
    #   {
    #     if (dgm_type == 1 | dgm_type == 3) {
    #       fit <- lme4::lmer(Y ~ Z * X + (1 + X| userid), data = dta)
    #     } else if (dgm_type == 2) {
    #       fit <- lme4::lmer(Y ~ Z * X + (Z * X| userid), data = dta)
    #     }
    #     list(coef = summary(fit)$coefficients, varcor = summary(fit)$varcor)
    #   },
    #   error = function(cond) {
    #     message("\nCatched error in lmer():")
    #     message(cond)
    #     return(list(coef = coef_NA_fill, varcor = varcor_NA_fill))
    #   })
    
    solution <- tryCatch(
      {
        if (dgm_type == 1 | dgm_type == 3) {
          fit <- lmer(Y ~ Z * X + (1 + X| userid), data = dta)
        } else if (dgm_type == 2) {
          fit <- lmer(Y ~ Z * X + (Z * X| userid), data = dta)
        }
        list(coef = summary(fit)$coefficients, varcor = summary(fit)$varcor)
      },
      error = function(cond) {
        message("\nCaught error in lmer():")
        message(cond)
        message("Skipping this iteration and moving to the next one.")
        return(list(coef = coef_NA_fill, varcor = varcor_NA_fill))
      })
    
    output <- list(solution)
  }
  
  sink()
  
  dir.create("simulation_result", showWarnings = FALSE)
  saveRDS(result, file = paste0("simulation_result/", idesign, ".RDS"))
}

# Stop the cluster after the parallel tasks are done
stopCluster(cl)
