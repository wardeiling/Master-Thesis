
# 2018.10.22 Tianchen Qian

# simulation study for the random effects model paper (Statistical Science)

# generative model:
# Y_t+1 = alpha_0 + alpha_1 X_t + b_0i + b_1i X_t + A_t (beta_0 + beta_1 X_t + b_2i + b_3i X_t) + epsilon_it

rm(list = ls())

set.seed(123) # set global seed
runname <- "GM123ad-1000reps-researchreport" # set a runname

# make a directory in simulation_results based on runname
dir.create(paste0("simulation_results/", runname), showWarnings = FALSE)

# load packages
# library(rootSolve) # OLD
library(geepack)
library(lme4)
# library(mvtnorm) # OLD
library(foreach)
library(doParallel) # NEW
library(doRNG)

# source the generative model
source("scripts/updated scripts/generative_model_updated_newmodels.R")

### Simulation -----------------------------------------

# NEW: Set up parallel backend
cores <- detectCores()
cl <- makeCluster(cores - 1, outfile = paste0("simulation_results/", runname, "/log.txt"))
registerDoParallel(cl)

# set the number of simulations
nsim <- 1000

# simulation in Section 4
# design <- expand.grid(sample_size = c(30, 100, 200), total_T = c(10, 30), dgm_type = 1:3)

# simulation with larger sample sizes
# design <- expand.grid(sample_size = c(30, 100, 200), total_T = c(10, 30), dgm_type = 2)

# simulation with "a" models
# design <- expand.grid(sample_size = c(30, 100, 200), total_T = c(10, 30), dgm_type = c("1a", "2a", "3a"))

# simulation with "b" models
# design <- expand.grid(sample_size = c(30, 100, 200), total_T = c(10, 30), dgm_type = c("1b", "2b", "3b"))

# simulation for all models with N = 200, 1000, T = 10, 30
# design <- expand.grid(sample_size = c(200, 1000), total_T = c(10, 30), dgm_type = c(1, "1a", "1b", 2, "2a", "2b", 3, "3a", "3b"))

# simulation for Research Report
design <- expand.grid(sample_size = c(30, 100, 200), total_T = c(10, 30), dgm_type = c(1,2,3,"3a","3d"))

# simulation for models 3, 3a, 3b, 3c, 3d with N = 1000, T = 10, 30
# design <- expand.grid(sample_size = c(200), total_T = c(10, 30), dgm_type = c(3, "3c"))

# make sure dgm_type is not a factor
design$dgm_type <- as.character(design$dgm_type)

for (idesign in 1:nrow(design)) {
    sample_size <- design$sample_size[idesign]
    total_T <- design$total_T[idesign]
    dgm_type <- design$dgm_type[idesign]
    
    ### create template output structure
    
    # for models without interaction beta_1
    if (dgm_type %in% c("1b", "2b", "3b", "3d")) {
        coef_NA_fill_mlm <- matrix(NA, ncol = 3, 
                                   nrow = 3, 
                                   dimnames = list(c("(Intercept)", "X", "A"), 
                                                   c("Estimate", "Std. Error", "t value")))
        coef_NA_fill_gee <- matrix(NA, ncol = 4, nrow = 3, 
                                   dimnames = list(c("(Intercept)", "X", "A"), 
                                                   c("Estimate", "Std.err", "Wald", "Pr(>|W|)")))
        
    # for model with interaction but without alpha_1
    } else if (dgm_type %in% c("3h")) {
        coef_NA_fill_mlm <- matrix(NA, ncol = 3, 
                                   nrow = 3, 
                                   dimnames = list(c("(Intercept)", "A", "A:X"), 
                                                   c("Estimate", "Std. Error", "t value")))
        coef_NA_fill_gee <- matrix(NA, ncol = 4, nrow = 3, 
                                   dimnames = list(c("(Intercept)", "A", "A:X"), 
                                                   c("Estimate", "Std.err", "Wald", "Pr(>|W|)")))
    
    # for all other models with all fixed effects
    } else {
        row_mlm <- c("(Intercept)", "X", "A", "X:A")
        col_mlm <- c("Estimate", "Std. Error", "t value")
        coef_NA_fill_mlm <- matrix(NA, ncol = 3, 
                                   nrow = 4, 
                                   dimnames = list(c("(Intercept)", "X", "A", "X:A"), 
                                                   c("Estimate", "Std. Error", "t value")))
        coef_NA_fill_gee <- matrix(NA, ncol = 4, nrow = 4, 
                                   dimnames = list(c("(Intercept)", "X", "A", "X:A"), 
                                                   c("Estimate", "Std.err", "Wald", "Pr(>|W|)")))
    }
    
    # set.seed(123) # original study had this here, but having a loop within the for loop
    # potentially makes it so that we get the same dataset for each of the parallel processes, 
    # which is undesirable as it is not really random.
    
    result <- foreach(isim = 1:nsim, .combine = "c", .errorhandling = "remove", 
                      .packages = c("geepack", "lme4"), .options.RNG=120) %dorng% {
        if (isim %% 10 == 0) {
            cat(paste("Starting iteration",isim, ", GM", dgm_type, ", N =", sample_size, ", T =", total_T, "\n"))
        }
        dta <- dgm_with_treatment(sample_size, total_T, dgm_type = dgm_type)
        
        ### Data Analysis
        
        # for mlm
        solution_lmm <- tryCatch(
          {
            if (dgm_type %in% c(1,3,"3c","3e","3f","3g")) {
              mlm_fit <- lmer(Y ~ X * A + (1 + A| userid), data = dta)
            } else if (dgm_type == 2) {
              mlm_fit <- lmer(Y ~ X * A + (X * A| userid), data = dta)
            } else if (dgm_type %in% c("1a", "2a", "3a")) {
              mlm_fit <- lmer(Y ~ X * A + (1 | userid), data = dta)
            } else if (dgm_type %in% c("1b", "2b", "3b")) {
              mlm_fit <- lmer(Y ~ X + A + (1 | userid), data = dta)
            } else if (dgm_type == "3d") {
              mlm_fit <- lmer(Y ~ X + A + (1 + A | userid), data = dta)
            } else if (dgm_type == "3h") {
              mlm_fit <- lmer(Y ~ A + A:X + (1 + A | userid), data = dta)
            } 
            
            # Check for singular fits
            if (isSingular(mlm_fit, tol = 1e-05)) {
              message("\nModel is singular.")
              return(list(coef = coef_NA_fill_mlm))
            }
            
            list(coef = summary(mlm_fit)$coefficients) #, varcor = summary(fit)$varcor)
            
          }, warning = function(w) {
            # log the warning message
            message("\nCaught warning in lmer():")
            message(w)
            
            # handle specific warnings
            if (grepl("Model failed to converge", w$message)) {
              return(list(coef = coef_NA_fill_mlm))
            }
            
            # For other warnings, return NA coefficients
            return(list(coef = coef_NA_fill_mlm))
            
          }, error = function(cond) {
            message("\nCaught error in lmer():")
            message(cond)
            return(list(coef = coef_NA_fill_mlm)) #, varcor = varcor_NA_fill))
          })
        
        # for gee with independence
        solution_gee_ind <- tryCatch(
            {
              if (dgm_type %in% c("1b", "2b", "3b", "3d")) {
                gee_ind_fit <- geeglm(Y ~ X + A, id = userid, family = gaussian, 
                                      corstr = "independence", data = dta)
              } else if (dgm_type == "3h") {
                gee_ind_fit <- geeglm(Y ~ A + A:X, id = userid, family = gaussian, 
                                      corstr = "independence", data = dta)
              } else {
                gee_ind_fit <- geeglm(Y ~ X * A, id = userid, family = gaussian, 
                                      corstr = "independence", data = dta)
              }
              
              list(coef = summary(gee_ind_fit)$coefficients)
            },
            error = function(cond) {
              message(paste("\nCaught error in geeglm() with independence:"))
              message(cond)
              return(list(coef = coef_NA_fill_gee)) 
            }
          )
        
        # for gee with exchangeable
        solution_gee_ex <- tryCatch(
          {
            if (dgm_type %in% c("1b", "2b", "3b", "3d")) {
              gee_ex_fit <- geeglm(Y ~ X + A, id = userid, family = gaussian, 
                                    corstr = "exchangeable", data = dta)
            } else if (dgm_type == "3h"){
              gee_ex_fit <- geeglm(Y ~ A + A:X, id = userid, family = gaussian, 
                                    corstr = "exchangeable", data = dta)
            } else {
              gee_ex_fit <- geeglm(Y ~ X * A, id = userid, family = gaussian, 
                                    corstr = "exchangeable", data = dta)
            }
            
            list(coef = summary(gee_ex_fit)$coefficients)
          },
          error = function(cond) {
            message(paste("\nCaught error in geeglm() with exchangeable:"))
            message(cond)
            return(list(coef = coef_NA_fill_gee)) 
          }
        )
        
        # for gee with AR(1)
        solution_gee_ar1 <- tryCatch(
          {
            if (dgm_type %in% c("1b", "2b", "3b", "3d")) {
              gee_ar1_fit <- geeglm(Y ~ X + A, id = userid, family = gaussian, 
                                   corstr = "ar1", data = dta)
            } else if (dgm_type == "3h") {
              gee_ar1_fit <- geeglm(Y ~ A + A:X, id = userid, family = gaussian, 
                                   corstr = "ar1", data = dta)
            } else {
              gee_ar1_fit <- geeglm(Y ~ X * A, id = userid, family = gaussian, 
                                   corstr = "ar1", data = dta)
            }

            list(coef = summary(gee_ar1_fit)$coefficients)
          },
          error = function(cond) {
            message(paste("\nCaught error in geeglm() with ar1:"))
            message(cond)
            return(list(coef = coef_NA_fill_gee)) 
          }
        )
        
        # store results
        output <- list(solution_lmm = solution_lmm, solution_gee_ind = solution_gee_ind,
                       solution_gee_ex = solution_gee_ex, solution_gee_ar1 = solution_gee_ar1)
      
    }
    
    saveRDS(result, file = paste0("simulation_results/", runname, "/", idesign, ".RDS"))
}

# Stop the cluster
# stopCluster(cl)

### collect results ---------------------------------------------------------

# create dataframe to store results
design$mlm_alpha_0_bias <- design$mlm_alpha_0_sd <- 
    design$mlm_alpha_1_bias <- design$mlm_alpha_1_sd <- 
    design$mlm_beta_0_bias <- design$mlm_beta_0_sd <- 
    design$mlm_beta_1_bias <- design$mlm_beta_1_sd <- 
design$gee_ind_alpha_0_bias <- design$gee_ind_alpha_0_sd <- 
    design$gee_ind_alpha_1_bias <- design$gee_ind_alpha_1_sd <- 
    design$gee_ind_beta_0_bias <- design$gee_ind_beta_0_sd <- 
    design$gee_ind_beta_1_bias <- design$gee_ind_beta_1_sd <-
design$gee_ex_alpha_0_bias <- design$gee_ex_alpha_0_sd <-
    design$gee_ex_alpha_1_bias <- design$gee_ex_alpha_1_sd <- 
    design$gee_ex_beta_0_bias <- design$gee_ex_beta_0_sd <- 
    design$gee_ex_beta_1_bias <- design$gee_ex_beta_1_sd <- 
design$gee_ar1_alpha_0_bias <- design$gee_ar1_alpha_0_sd <-
    design$gee_ar1_alpha_1_bias <- design$gee_ar1_alpha_1_sd <- 
    design$gee_ar1_beta_0_bias <- design$gee_ar1_beta_0_sd <- 
    design$gee_ar1_beta_1_bias <- design$gee_ar1_beta_1_sd <- 
design$mlm_success <- design$gee_ind_success <- design$gee_ex_success <- design$gee_ar1_success <- NA

for (idesign in 1:nrow(design)) {
  
    results_list <- readRDS(paste0("simulation_results/", runname, "/", idesign, ".RDS"))
    
    dgm_type <- design$dgm_type[idesign]
    
    ### Restate the true values of the parameters ###
    
    # dgm_type = 1 or 3
    alpha_0_true <- - 2 # was originally -1 in the code but is -2 in the paper
    alpha_1_true <- - 0.3
    beta_0_true <- 1 # was originallly 0.5 in the code but is 1 in the paper
    beta_1_true <- 0.3 # was originally 0.1 in the code but is 0.3 in the paper
    
    # # set the interactions to 0 for "b" models
    # if (dgm_type %in% c("1b", "2b", "3b", "3d")) {
    #   beta_1_true <- 0
    # }
    
    ### Divide the list into sublists for each model type ###
    
    result_lmm <- results_list[grep("solution_lmm", names(results_list))]
    result_gee_ind <- results_list[grep("solution_gee_ind", names(results_list))]
    result_gee_ex <- results_list[grep("solution_gee_ex", names(results_list))]
    result_gee_ar1 <- results_list[grep("solution_gee_ar1", names(results_list))]
    
    ### Extract the Results ###
    
    ### Linear Mixed Model
    
    # mlm_alpha_0 <- sapply(result_lmm, function(l) l$coef["(Intercept)", "Estimate"])
    # mlm_alpha_0_sd <- sapply(result_lmm, function(l) l$coef["(Intercept)", "Std. Error"])
    # mlm_alpha_1 <- sapply(result_lmm, function(l) l$coef["X", "Estimate"])
    # mlm_alpha_1_sd <- sapply(result_lmm, function(l) l$coef["X", "Std. Error"])
    mlm_beta_0 <- sapply(result_lmm, function(l) l$coef["A", "Estimate"])
    mlm_beta_0_sd <- sapply(result_lmm, function(l) l$coef["A", "Std. Error"])
    
    # if (!dgm_type %in% c("1b", "2b", "3b", "3d")) {
    #   mlm_beta_1 <- sapply(result_lmm, function(l) l$coef["X:A", "Estimate"])
    #   mlm_beta_1_sd <- sapply(result_lmm, function(l) l$coef["X:A", "Std. Error"])
    # }

    # design$mlm_alpha_0_bias[idesign] <- mean(mlm_alpha_0, na.rm = T) - alpha_0_true
    # design$mlm_alpha_0_sd[idesign] <- sd(mlm_alpha_0, na.rm = T)
    # design$mlm_alpha_1_bias[idesign] <- mean(mlm_alpha_1, na.rm = T) - alpha_1_true
    # design$mlm_alpha_1_sd[idesign] <- sd(mlm_alpha_1, na.rm = T)
    design$mlm_beta_0_bias[idesign] <- mean(mlm_beta_0, na.rm = T) - beta_0_true
    design$mlm_beta_0_sd[idesign] <- sd(mlm_beta_0, na.rm = T)
    
    # # interaction beta1 is not present in the "b" models, so only calculate if not b
    # if (!dgm_type %in% c("1b", "2b", "3b", "3d")) {
    #     design$mlm_beta_1_bias[idesign] <- mean(mlm_beta_1, na.rm = T) - beta_1_true
    #     design$mlm_beta_1_sd[idesign] <- sd(mlm_beta_1, na.rm = T)
    # }
    
    # check if models converged by calculating success proportion for beta_0
    # so the value should not be NA or 0
    design$mlm_success[idesign] <- 1 - sum(is.na(mlm_beta_0)) / nsim
    
    ### GEE with independence
    
    # gee_ind_alpha_0 <- sapply(result_gee_ind, function(l) l$coef["(Intercept)", "Estimate"])
    # gee_ind_alpha_0_sd <- sapply(result_gee_ind, function(l) l$coef["(Intercept)", "Std.err"])
    # gee_ind_alpha_1 <- sapply(result_gee_ind, function(l) l$coef["X", "Estimate"])
    # gee_ind_alpha_1_sd <- sapply(result_gee_ind, function(l) l$coef["X", "Std.err"])
    gee_ind_beta_0 <- sapply(result_gee_ind, function(l) l$coef["A", "Estimate"])
    gee_ind_beta_0_sd <- sapply(result_gee_ind, function(l) l$coef["A", "Std.err"])
    
    # if (!dgm_type %in% c("1b", "2b", "3b", "3d")) {
    #   gee_ind_beta_1 <- sapply(result_gee_ind, function(l) l$coef["X:A", "Estimate"])
    #   gee_ind_beta_1_sd <- sapply(result_gee_ind, function(l) l$coef["X:A", "Std.err"])
    # }
    
    # design$gee_ind_alpha_0_bias[idesign] <- mean(gee_ind_alpha_0, na.rm = T) - alpha_0_true
    # design$gee_ind_alpha_0_sd[idesign] <- sd(gee_ind_alpha_0, na.rm = T)
    # design$gee_ind_alpha_1_bias[idesign] <- mean(gee_ind_alpha_1, na.rm = T) - alpha_1_true
    # design$gee_ind_alpha_1_sd[idesign] <- sd(gee_ind_alpha_1, na.rm = T)
    design$gee_ind_beta_0_bias[idesign] <- mean(gee_ind_beta_0, na.rm = T) - beta_0_true
    design$gee_ind_beta_0_sd[idesign] <- sd(gee_ind_beta_0, na.rm = T)
    
    # if (!dgm_type %in% c("1b", "2b", "3b", "3d")) {
    #   design$gee_ind_beta_1_bias[idesign] <- mean(gee_ind_beta_1, na.rm = T) - beta_1_true
    #   design$gee_ind_beta_1_sd[idesign] <- sd(gee_ind_beta_1, na.rm = T)
    # }
    
    # check if models converged by calculating success proportion for beta_0
    design$gee_ind_success[idesign] <- 1 - sum(is.na(gee_ind_beta_0)) / nsim
    
    ### GEE with exchangeable
    
    # gee_ex_alpha_0 <- sapply(result_gee_ex, function(l) l$coef["(Intercept)", "Estimate"])
    # gee_ex_alpha_0_sd <- sapply(result_gee_ex, function(l) l$coef["(Intercept)", "Std.err"])
    # gee_ex_alpha_1 <- sapply(result_gee_ex, function(l) l$coef["X", "Estimate"])
    # gee_ex_alpha_1_sd <- sapply(result_gee_ex, function(l) l$coef["X", "Std.err"])
    gee_ex_beta_0 <- sapply(result_gee_ex, function(l) l$coef["A", "Estimate"])
    gee_ex_beta_0_sd <- sapply(result_gee_ex, function(l) l$coef["A", "Std.err"])
    
    # if (!dgm_type %in% c("1b", "2b", "3b", "3d")) {
    #   gee_ex_beta_1 <- sapply(result_gee_ex, function(l) l$coef["X:A", "Estimate"])
    #   gee_ex_beta_1_sd <- sapply(result_gee_ex, function(l) l$coef["X:A", "Std.err"])
    # }
    
    # design$gee_ex_alpha_0_bias[idesign] <- mean(gee_ex_alpha_0, na.rm = T) - alpha_0_true
    # design$gee_ex_alpha_0_sd[idesign] <- sd(gee_ex_alpha_0, na.rm = T)
    # design$gee_ex_alpha_1_bias[idesign] <- mean(gee_ex_alpha_1, na.rm = T) - alpha_1_true
    # design$gee_ex_alpha_1_sd[idesign] <- sd(gee_ex_alpha_1, na.rm = T)
    design$gee_ex_beta_0_bias[idesign] <- mean(gee_ex_beta_0, na.rm = T) - beta_0_true
    design$gee_ex_beta_0_sd[idesign] <- sd(gee_ex_beta_0, na.rm = T)
    
    # if (!dgm_type %in% c("1b", "2b", "3b", "3d")) {
    #   design$gee_ex_beta_1_bias[idesign] <- mean(gee_ex_beta_1, na.rm = T) - beta_1_true
    #   design$gee_ex_beta_1_sd[idesign] <- sd(gee_ex_beta_1, na.rm = T)
    # }
    
    # check if models converged by calculating success proportion for beta_0
    design$gee_ex_success[idesign] <- 1 - sum(is.na(gee_ex_beta_0)) / nsim
    
    ### GEE with AR(1)
    
    # gee_ar1_alpha_0 <- sapply(result_gee_ar1, function(l) l$coef["(Intercept)", "Estimate"])
    # gee_ar1_alpha_0_sd <- sapply(result_gee_ar1, function(l) l$coef["(Intercept)", "Std.err"])
    # gee_ar1_alpha_1 <- sapply(result_gee_ar1, function(l) l$coef["X", "Estimate"])
    # gee_ar1_alpha_1_sd <- sapply(result_gee_ar1, function(l) l$coef["X", "Std.err"])
    gee_ar1_beta_0 <- sapply(result_gee_ar1, function(l) l$coef["A", "Estimate"])
    gee_ar1_beta_0_sd <- sapply(result_gee_ar1, function(l) l$coef["A", "Std.err"])
    
    # if (!dgm_type %in% c("1b", "2b", "3b", "3d")) {
    #   gee_ar1_beta_1 <- sapply(result_gee_ar1, function(l) l$coef["X:A", "Estimate"])
    #   gee_ar1_beta_1_sd <- sapply(result_gee_ar1, function(l) l$coef["X:A", "Std.err"])
    # }
    
    # design$gee_ar1_alpha_0_bias[idesign] <- mean(gee_ar1_alpha_0, na.rm = T) - alpha_0_true
    # design$gee_ar1_alpha_0_sd[idesign] <- sd(gee_ar1_alpha_0, na.rm = T)
    # design$gee_ar1_alpha_1_bias[idesign] <- mean(gee_ar1_alpha_1, na.rm = T) - alpha_1_true
    # design$gee_ar1_alpha_1_sd[idesign] <- sd(gee_ar1_alpha_1, na.rm = T)
    design$gee_ar1_beta_0_bias[idesign] <- mean(gee_ar1_beta_0, na.rm = T) - beta_0_true
    design$gee_ar1_beta_0_sd[idesign] <- sd(gee_ar1_beta_0, na.rm = T)
    
    # if (!dgm_type %in% c("1b", "2b", "3b", "3d")) {
    #   design$gee_ar1_beta_1_bias[idesign] <- mean(gee_ar1_beta_1, na.rm = T) - beta_1_true
    #   design$gee_ar1_beta_1_sd[idesign] <- sd(gee_ar1_beta_1, na.rm = T)
    # }
    
    # check if models converged by calculating success proportion for beta_0
    design$gee_ar1_success[idesign] <- 1 - sum(is.na(gee_ar1_beta_0)) / nsim
    
}

### Post-processing of Results ----------------------------------------------

# reorder the columns dmg_type, total_T, sample_size to be first
design2 <- design[, c(3, 2, 1, 4:ncol(design))]
# change dmg_type to GM, total_T to T, sample_size to N
colnames(design2) <- c("GM", "T", "N", colnames(design2)[4:ncol(design2)])
colnames(design2) <- gsub("beta_0", "beta0", colnames(design2))
# colnames(design2) <- gsub("beta_1", "beta1", colnames(design2))
# colnames(design2) <- gsub("alpha_0", "alpha0", colnames(design2))
# colnames(design2) <- gsub("alpha_1", "alpha1", colnames(design2))
saveRDS(design2, paste0("simulation_results/", runname, "/results_all.RDS"))

# exclude any column that doesn't include "beta_0", the first three columns and success rates
design3a <- design2[, c(1:3, grep("beta0", colnames(design2)), grep("success", colnames(design2)))]
design3a <- design3a[, c("GM", "T", "N", "mlm_beta0_bias", "mlm_beta0_sd", "gee_ex_beta0_bias", "gee_ex_beta0_sd", 
                         "gee_ar1_beta0_bias", "gee_ar1_beta0_sd", "gee_ind_beta0_bias", "gee_ind_beta0_sd", 
                         "mlm_success", "gee_ex_success", "gee_ar1_success", "gee_ind_success")]
saveRDS(design3a, paste0("simulation_results/", runname, "/results_beta0_bias_sd_success.RDS"))

# exclude any column that doesn't include "beta_0", the first three columns
design3b <- design2[, c(1:3, grep("beta0", colnames(design2)))]
design3c <- design3b[, c("GM", "T", "N", "mlm_beta0_bias", "mlm_beta0_sd", "gee_ex_beta0_bias", "gee_ex_beta0_sd", 
                         "gee_ar1_beta0_bias", "gee_ar1_beta0_sd", "gee_ind_beta0_bias", "gee_ind_beta0_sd")]
saveRDS(design3c, paste0("simulation_results/", runname, "/results_beta0_bias_sd.RDS"))

# exclude any column that doesn't include "beta_0_bias" and the first three columns
design4a <- design2[, c(1:3, grep("beta0_bias", colnames(design2)))]
design4b <- design4a[, c("GM", "T", "N", "mlm_beta0_bias", "gee_ex_beta0_bias", "gee_ar1_beta0_bias", "gee_ind_beta0_bias")]
saveRDS(design4b, paste0("simulation_results/", runname, "/results_beta0_bias.RDS"))

# make LaTeX tables --------------------------------------------------------

library(xtable)

# make table for beta0 bias with Standarddeviation and success rate
colnames(design3a) <- c("GM", "T", "N", "MLM_bias", "MLM_sd", "GEE-Ex_bias", "GEE-Ex_sd", 
                        "GEE-AR1_bias", "GEE-AR1_sd", "GEE-Ind_bias", "GEE-Ind_sd", "MLM_success", "GEE-Ex_success", "GEE-AR1_success", "GEE-Ind_success")
print(xtable(design3a, digits = c(0, 0, 0, 0, rep(3, 8), rep(3, 4)), 
             caption = paste0("Results for beta0 bias with Standarddeviation and success rate, ", nsim, " replications, run: ", runname), label = "tab:beta0_bias_sd_success"), 
      include.rownames = FALSE, file = paste0("simulation_results/", runname, "/results_beta0_bias_sd_success.tex"))

# make table for beta0 bias with Standarddeviation
design3d <- design3c
colnames(design3d) <- c("GM", "T", "N", "MLM_bias", "MLM_sd", "GEE-Ex_bias", "GEE-Ex_sd", 
                        "GEE-AR1_bias", "GEE-AR1_sd", "GEE-Ind_bias", "GEE-Ind_sd")
print(xtable(design3d, digits = c(0, 0, 0, 0, rep(3, 8)), 
             caption = paste0("Results for beta0 bias with Standarddeviation, ", nsim, " replications, run: ", runname), label = "tab:beta0_bias_sd"), 
      include.rownames = FALSE, file = paste0("simulation_results/", runname, "/results_beta0_bias_sd.tex"))

# make table for beta0 bias
design4c <- design4b
colnames(design4c) <- c("GM", "T", "N", "MLM", "GEE-Ex", "GEE-AR1", "GEE-Ind")

print(xtable(design4c, digits = c(0, 0, 0, 0, rep(3, 4)), 
             caption = paste0("Results for beta0 bias, ", nsim, " replications, run: ", runname), label = "tab:beta0_bias"),
      include.rownames = FALSE, file = paste0("simulation_results/", runname, "/results_beta0_bias.tex"))
