
# 2018.10.22 Tianchen Qian

# simulation study for the random effects model paper (Statistical Science)

# generative model:
# Y_t+1 = alpha_0 + alpha_1 X_t + b_0i + b_1i X_t + A_t (beta_0 + beta_1 X_t + b_2i + b_3i X_t) + epsilon_it

rm(list = ls())

set.seed(123) # set global seed
runname <- "GM123ad-1000reps-researchreport-cleanscript_nrep10000" # set a runname

# make a directory in simulation_results based on runname
dir.create(paste0("simulation_results/", runname), showWarnings = FALSE)

# load packages
library(geepack)
library(lme4)
library(foreach)
library(doParallel) # NEW
library(doRNG)

# source the generative model
source("scripts/updated scripts/simplified-research-report/generative_model_updated_newmodels_RR.R")

### Simulation -----------------------------------------

# NEW: Set up parallel backend
cores <- detectCores()
cl <- makeCluster(cores - 1, outfile = paste0("simulation_results/", runname, "/log.txt"))
registerDoParallel(cl)

# set the number of simulations
nsim <- 10000

# simulation for Research Report
design <- expand.grid(sample_size = c(30, 100, 200), total_T = c(10, 30), dgm_type = c(1,2,3,"3a","3d"))

# make sure dgm_type is not a factor
design$dgm_type <- as.character(design$dgm_type)

for (idesign in 1:nrow(design)) {
    sample_size <- design$sample_size[idesign]
    total_T <- design$total_T[idesign]
    dgm_type <- design$dgm_type[idesign]
    
    ### create template output structure
    
    # for models without interaction beta_1
    if (dgm_type %in% c("3d")) {
        coef_NA_fill_mlm <- matrix(NA, ncol = 3, 
                                   nrow = 3, 
                                   dimnames = list(c("(Intercept)", "X", "A"), 
                                                   c("Estimate", "Std. Error", "t value")))
        coef_NA_fill_gee <- matrix(NA, ncol = 4, nrow = 3, 
                                   dimnames = list(c("(Intercept)", "X", "A"), 
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
            if (dgm_type %in% c(1,3)) {
              mlm_fit <- lmer(Y ~ X * A + (1 + A| userid), data = dta)
            } else if (dgm_type == 2) {
              mlm_fit <- lmer(Y ~ X * A + (X * A| userid), data = dta)
            } else if (dgm_type == "3a") {
              mlm_fit <- lmer(Y ~ X * A + (1 | userid), data = dta)
            } else if (dgm_type == "3d") {
              mlm_fit <- lmer(Y ~ X + A + (1 + A | userid), data = dta)
            }
            
            # Check for singular fits
            if (isSingular(mlm_fit, tol = 1e-05)) {
              message("\nModel is singular.")
              return(list(coef = coef_NA_fill_mlm))
            }
            
            list(coef = summary(mlm_fit)$coefficients)
            
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
            return(list(coef = coef_NA_fill_mlm))
          })
        
        # for gee with independence
        solution_gee_ind <- tryCatch(
            {
              if (dgm_type == "3d") {
                gee_ind_fit <- geeglm(Y ~ X + A, id = userid, family = gaussian, 
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
            if (dgm_type == "3d") {
              gee_ex_fit <- geeglm(Y ~ X + A, id = userid, family = gaussian, 
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
            if (dgm_type == "3d") {
              gee_ar1_fit <- geeglm(Y ~ X + A, id = userid, family = gaussian, 
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
design$mlm_beta_0_bias <- design$gee_ind_beta_0_bias <- design$gee_ex_beta_0_bias <- design$gee_ar1_beta_0_bias <-
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
    
    ### Divide the list into sublists for each model type ###
    
    result_lmm <- results_list[grep("solution_lmm", names(results_list))]
    result_gee_ind <- results_list[grep("solution_gee_ind", names(results_list))]
    result_gee_ex <- results_list[grep("solution_gee_ex", names(results_list))]
    result_gee_ar1 <- results_list[grep("solution_gee_ar1", names(results_list))]
    
    ### Extract the Results ###
    
    ### Linear Mixed Model
    
    mlm_beta_0 <- sapply(result_lmm, function(l) l$coef["A", "Estimate"])
    design$mlm_beta_0_bias[idesign] <- mean(mlm_beta_0, na.rm = T) - beta_0_true
    design$mlm_beta_0_sd[idesign] <- sd(mlm_beta_0, na.rm = T)
    
    # check if models converged by calculating success proportion for beta_0
    # so the value should not be NA or 0
    design$mlm_success[idesign] <- 1 - sum(is.na(mlm_beta_0)) / nsim
    
    ### GEE with independence
    
    gee_ind_beta_0 <- sapply(result_gee_ind, function(l) l$coef["A", "Estimate"])
    design$gee_ind_beta_0_bias[idesign] <- mean(gee_ind_beta_0, na.rm = T) - beta_0_true
    design$gee_ind_beta_0_sd[idesign] <- sd(gee_ind_beta_0, na.rm = T)
    
    # check if models converged by calculating success proportion for beta_0
    design$gee_ind_success[idesign] <- 1 - sum(is.na(gee_ind_beta_0)) / nsim
    
    ### GEE with exchangeable
    
    gee_ex_beta_0 <- sapply(result_gee_ex, function(l) l$coef["A", "Estimate"])
    design$gee_ex_beta_0_bias[idesign] <- mean(gee_ex_beta_0, na.rm = T) - beta_0_true
    design$gee_ex_beta_0_sd[idesign] <- sd(gee_ex_beta_0, na.rm = T)
    
    # check if models converged by calculating success proportion for beta_0
    design$gee_ex_success[idesign] <- 1 - sum(is.na(gee_ex_beta_0)) / nsim
    
    ### GEE with AR(1)

    gee_ar1_beta_0 <- sapply(result_gee_ar1, function(l) l$coef["A", "Estimate"])
    design$gee_ar1_beta_0_bias[idesign] <- mean(gee_ar1_beta_0, na.rm = T) - beta_0_true
    design$gee_ar1_beta_0_sd[idesign] <- sd(gee_ar1_beta_0, na.rm = T)
    
    # check if models converged by calculating success proportion for beta_0
    design$gee_ar1_success[idesign] <- 1 - sum(is.na(gee_ar1_beta_0)) / nsim
    
}

### Post-processing of Results ----------------------------------------------

# reorder the columns dmg_type, total_T, sample_size to be first
design2 <- design[, c(3, 2, 1, 4:ncol(design))]
# change dmg_type to GM, total_T to T, sample_size to N
colnames(design2) <- c("GM", "T", "N", colnames(design2)[4:ncol(design2)])
colnames(design2) <- gsub("beta_0", "beta0", colnames(design2))
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
