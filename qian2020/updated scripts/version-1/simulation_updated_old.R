
# 2018.10.22 Tianchen Qian

# simulation study for the random effects model paper (Statistical Science)

# generative model:
# Y_t+1 = alpha_0 + alpha_1 X_t + b_0i + b_1i X_t + A_t (beta_0 + beta_1 X_t + b_2i + b_3i X_t) + epsilon_it


rm(list = ls())
set.seed(123)

library(rootSolve)
library(geepack)
library(lme4)
library(mvtnorm)
library(foreach)
# library(doMC)
library(doParallel) # NEW
library(doRNG)


source("qian2020/updated scripts/generative_model_updated.R")

### simulation: using lmer() packag -----------------------------------------


# if(0) {
    
# max_cores <- 16
# registerDoMC(min(detectCores() - 1, max_cores))

# NEW: Set up parallel backend
cores <- detectCores()
cl <- makeCluster(cores - 1, outfile = "log.txt")
registerDoParallel(cl)

# NEW: Set up parallel RNG
# doRNG::registerDoRNG(seed = 123)
set.seed(120)


nsim <- 100

# full simulation in Appendix B
# design <- expand.grid(sample_size = c(30, 50, 100, 200), total_T = c(10, 20, 30), dgm_type = 1:3)

# simulation in Section 4
design <- expand.grid(sample_size = c(30, 100, 200), total_T = c(10, 30), dgm_type = 1:3)

# only GM3, N = 200, T = 10
# design <- expand.grid(sample_size = 200, total_T = 10, dgm_type = 3) # scenario 15 of section 4

design <- design[order(design$dgm_type), ]

for (idesign in 1:nrow(design)) {
    sample_size <- design$sample_size[idesign]
    total_T <- design$total_T[idesign]
    dgm_type <- design$dgm_type[idesign]
    
    # create template output structure for simulated trials with estimation error
    
    # for mlm
    row_mlm <- c("(Intercept)", "X", "A", "X:A")
    col_mlm <- c("Estimate", "Std. Error", "t value")
    coef_NA_fill_mlm <- matrix(NA, ncol = length(col_mlm), nrow = length(row_mlm), dimnames = list(row_mlm, col_mlm))
    
    # for gee
    row_gee <- c("(Intercept)", "X", "A", "X:A")
    col_gee <- c("Estimate", "Std.err", "Wald", "Pr(>|W|)")
    coef_NA_fill_gee <- matrix(NA, ncol = length(col_gee), nrow = length(row_gee), dimnames = list(row_gee, col_gee))
    
    # dta <- dgm_with_treatment(sample_size, total_T, dgm_type = dgm_type)
    # 
    # coef_NA_fill <- matrix(NA, ncol = 3, nrow = 4)
    # if (dgm_type %in% c(1,3,4)) {
    #     fit <- lmer(Y ~ X * A + (1 + A| userid), data = dta)
    # } else if (dgm_type == 2) {
    #     fit <- lmer(Y ~ X * A + (X * A| userid), data = dta)
    # }
    # coef <- summary(fit)$coefficients
    # for (irow in 1:nrow(coef)) {
    #     for (icol in 1:ncol(coef)) {
    #         coef[irow, icol] <- NA
    #     }
    # }
    # 
    # coef_NA_fill <- coef
    
    # start parallel jobs
    # writeLines(c(""), "log.txt")
    # sink("log.txt", append=FALSE)
    set.seed(123)
    result <- foreach(isim = 1:nsim, .combine = "c", .errorhandling = "remove") %dorng% {
        if (isim %% 10 == 0) {
            cat(paste("Starting iteration",isim,"\n"))
        }
        dta <- dgm_with_treatment(sample_size, total_T, dgm_type = dgm_type)
        
        library(geepack) # NEW
        library(lme4) # NEW
        
        # Mixed model estimation
        solution_lmm <- tryCatch(
          {
            if (dgm_type == 1 | dgm_type == 3) {
              mlm_fit <- lmer(Y ~ X * A + (1 + A| userid), data = dta)
            } else if (dgm_type == 2) {
              mlm_fit <- lmer(Y ~ X * A + (X * A| userid), data = dta)
            }
            list(coef = summary(mlm_fit)$coefficients) #, varcor = summary(fit)$varcor)
          },
          error = function(cond) {
            message("\nCaught error in lmer():")
            message(cond)
            return(list(coef = coef_NA_fill_mlm)) #, varcor = varcor_NA_fill))
          })
        
        # GEE estimation
        solution_gee_ind <- tryCatch(
            {
              gee_ind_fit <- geeglm(Y ~ X * A, id = userid, family = gaussian, 
                                corstr = "independence", data = dta)
              list(coef = summary(gee_ind_fit)$coefficients)
            },
            error = function(cond) {
              message(paste("\nCaught error in geeglm() with independence:"))
              message(cond)
              return(list(coef = coef_NA_fill_gee)) 
            }
          )
        
        solution_gee_ex <- tryCatch(
          {
            gee_ex_fit <- geeglm(Y ~ X * A, id = userid, family = gaussian, 
                                  corstr = "exchangeable", data = dta)
            list(coef = summary(gee_ex_fit)$coefficients)
          },
          error = function(cond) {
            message(paste("\nCaught error in geeglm() with exchangeable:"))
            message(cond)
            return(list(coef = coef_NA_fill_gee)) 
          }
        )
        
        solution_gee_ar1 <- tryCatch(
          {
            gee_ar1_fit <- geeglm(Y ~ X * A, id = userid, family = gaussian, 
                                 corstr = "ar1", data = dta)
            list(coef = summary(gee_ar1_fit)$coefficients)
          },
          error = function(cond) {
            message(paste("\nCaught error in geeglm() with ar1:"))
            message(cond)
            return(list(coef = coef_NA_fill_gee)) 
          }
        )
        
        output <- list(solution_lmm = solution_lmm, solution_gee_ind = solution_gee_ind,
                       solution_gee_ex = solution_gee_ex, solution_gee_ar1 = solution_gee_ar1)
        
        # solution <- tryCatch(
        #     {
        #         if (dgm_type == 1 | dgm_type == 3) {
        #             fit <- lmer(Y ~ X * A + (1 + A| userid), data = dta)
        #         } else if (dgm_type == 2) {
        #             fit <- lmer(Y ~ X * A + (X * A| userid), data = dta)
        #         }
        #         list(coef = summary(fit)$coefficients) #varcor = summary(fit)$varcor)
        #     },
        #     error = function(cond) {
        #         message("\nCatched error in lmer():")
        #         message(cond)
        #         return(list(coef = coef_NA_fill)) # varcor = varcor_NA_fill))
        #     })
        # 
        # output <- list(solution)
    }
    # sink()
    
    dir.create("simulation_result", showWarnings = FALSE)
    saveRDS(result, file = paste0("simulation_result/", idesign, ".RDS"))
}

# Stop the cluster
stopCluster(cl)


### collect results ---------------------------------------------------------

# # NEW: extract result 15
# setting15 <- readRDS("qian2020/output/setting15.RDS")
# beta_0 <- sapply(setting15, function(l) l$coef["A", "Estimate"])
# beta_0_mean <- mean(beta_0)
# beta_0_bias <- beta_0_mean - 1
# beta_0_sd <- sd(beta_0)
# 
# results_list <- setting15 <- readRDS("simulation_result/15.RDS")
# beta_0 <- sapply(setting15, function(l) l$coef["A", "Estimate"])
# beta_0_mean <- mean(beta_0)
# beta_0_bias <- beta_0_mean - 1
# beta_0_sd <- sd(beta_0)
# 
# # Extract Linear Mixed Model solution
# lmm_solution <- setting15[["solution_lmm"]]
# 
# # divide beta_0 into each of the four estimation methods 
# # by the name of the list
# mlm_results <- 

# varcomp_result <- design
# 
# varcomp_result$sigma_b0_bias <- varcomp_result$sigma_b0_sd <- 
#     varcomp_result$sigma_b1_bias <- varcomp_result$sigma_b1_sd <- 
#     varcomp_result$sigma_b2_bias <- varcomp_result$sigma_b2_sd <- 
#     varcomp_result$sigma_b3_bias <- varcomp_result$sigma_b3_sd <- 
#     varcomp_result$sigma_eps_bias <- varcomp_result$sigma_eps_sd <- NA

# design$alpha_0_bias <- design$alpha_0_sd <- design$alpha_0_cp <- 
#     design$alpha_1_bias <- design$alpha_1_sd <- design$alpha_1_cp <- 
#     design$beta_0_bias <- design$beta_0_sd <- design$beta_0_cp <- 
#     design$beta_1_bias <- design$beta_1_sd <- design$beta_1_cp <- NA

# make new one now with gee too
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
    design$gee_ar1_beta_1_bias <- design$gee_ar1_beta_1_sd <- NA

for (idesign in 1:nrow(design)) {
  
    results_list <- readRDS(paste0("simulation_result/", idesign, ".RDS"))
    
    dgm_type <- design$dgm_type[idesign]
    
    ### Restate the true values of the parameters ###
    
    # dgm_type = 1 or 3
    alpha_0_true <- - 2 # was originally -1 in the code but is -2 in the paper
    alpha_1_true <- - 0.3
    beta_0_true <- 1 # was originallly 0.5 in the code but is 1 in the paper
    beta_1_true <- 0.3 # was originally 0.1 in the code but is 0.3 in the paper
    # sigma_b0_true <- 2
    # sigma_b1_true <- 0
    # sigma_b2_true <- 1
    # sigma_b3_true <- 0
    # sigma_eps_true <- 1
    
    if (dgm_type == 2) {
        sigma_b1_true <- sigma_b3_true <- 0.5
    }
    if (dgm_type == 4) {
        sigma_b2_true <- 0
    }
    
    ### Divide the list into sublists for each model type ###
    
    result_lmm <- results_list[grep("solution_lmm", names(results_list))]
    result_gee_ind <- results_list[grep("solution_gee_ind", names(results_list))]
    result_gee_ex <- results_list[grep("solution_gee_ex", names(results_list))]
    result_gee_ar1 <- results_list[grep("solution_gee_ar1", names(results_list))]
    
    ### Extract the Results ###
    
    ### Linear Mixed Model
    
    mlm_alpha_0 <- sapply(result_lmm, function(l) l$coef["(Intercept)", "Estimate"])
    mlm_alpha_0_sd <- sapply(result_lmm, function(l) l$coef["(Intercept)", "Std. Error"])
    # mlm_alpha_0_df <- sapply(result_lmm, function(l) l$coef["(Intercept)", "df"])
    mlm_alpha_1 <- sapply(result_lmm, function(l) l$coef["X", "Estimate"])
    mlm_alpha_1_sd <- sapply(result_lmm, function(l) l$coef["X", "Std. Error"])
    # mlm_alpha_1_df <- sapply(result_lmm, function(l) l$coef["X", "df"])
    mlm_beta_0 <- sapply(result_lmm, function(l) l$coef["A", "Estimate"])
    mlm_beta_0_sd <- sapply(result_lmm, function(l) l$coef["A", "Std. Error"])
    # mlm_beta_0_df <- sapply(result_lmm, function(l) l$coef["A", "df"])
    mlm_beta_1 <- sapply(result_lmm, function(l) l$coef["X:A", "Estimate"])
    mlm_beta_1_sd <- sapply(result_lmm, function(l) l$coef["X:A", "Std. Error"])
    # mlm_beta_1_df <- sapply(result_lmm, function(l) l$coef["X:A", "df"])
    
    # varcor <- lapply(result_lmm, function(l) l$varcor$userid) # variance matrix for random effects; this is "G"
    
    # sigma_eps <- sapply(result_lmm, function(l) attr(l$varcor, "sc"))
    
    # quantiles <- qt(0.975, alpha_0_df)
    design$mlm_alpha_0_bias[idesign] <- mean(mlm_alpha_0) - alpha_0_true
    design$mlm_alpha_0_sd[idesign] <- sd(mlm_alpha_0)
    # design$alpha_0_cp[idesign] <- mean((alpha_0_true < alpha_0 + quantiles * alpha_0_sd) & (alpha_0_true > alpha_0 - quantiles * alpha_0_sd))
    # design$alpha_0_cp[idesign] <- mean((alpha_0_true < alpha_0 + 1.96 * design$alpha_0_sd[idesign]) & (alpha_0_true > alpha_0 - 1.96 * design$alpha_0_sd[idesign]))
    
    # quantiles <- qt(0.975, alpha_1_df)
    design$mlm_alpha_1_bias[idesign] <- mean(mlm_alpha_1) - alpha_1_true
    design$mlm_alpha_1_sd[idesign] <- sd(mlm_alpha_1)
    # design$alpha_1_cp[idesign] <- mean((alpha_1_true < alpha_1 + quantiles * alpha_1_sd) & (alpha_1_true > alpha_1 - quantiles * alpha_1_sd))
    # design$alpha_1_cp[idesign] <- mean((alpha_1_true < alpha_1 + 1.96 * design$alpha_1_sd[idesign]) & (alpha_1_true > alpha_1 - 1.96 * design$alpha_1_sd[idesign]))
    
    # quantiles <- qt(0.975, beta_0_df)
    design$mlm_beta_0_bias[idesign] <- mean(mlm_beta_0) - beta_0_true
    design$mlm_beta_0_sd[idesign] <- sd(mlm_beta_0)
    # design$beta_0_cp[idesign] <- mean((beta_0_true < beta_0 + quantiles * beta_0_sd) & (beta_0_true > beta_0 - quantiles * beta_0_sd))
    # design$beta_0_cp[idesign] <- mean((beta_0_true < beta_0 + 1.96 * design$beta_0_sd[idesign]) & (beta_0_true > beta_0 - 1.96 * design$beta_0_sd[idesign]))
    
    # quantiles <- qt(0.975, beta_1_df)
    design$mlm_beta_1_bias[idesign] <- mean(mlm_beta_1) - beta_1_true
    design$mlm_beta_1_sd[idesign] <- sd(mlm_beta_1)
    # design$beta_1_cp[idesign] <- mean((beta_1_true < beta_1 + quantiles * beta_1_sd) & (beta_1_true > beta_1 - quantiles * beta_1_sd))
    # design$beta_1_cp[idesign] <- mean((beta_1_true < beta_1 + 1.96 * design$beta_1_sd[idesign]) & (beta_1_true > beta_1 - 1.96 * design$beta_1_sd[idesign]))
    
    # varcor_mean <- apply(simplify2array(varcor), 1:2, mean)
    # varcor_sd <- apply(simplify2array(varcor), 1:2, sd)
    # if (dgm_type == 1 | dgm_type == 3) {
    #     varcomp_result_lmm$sigma_b0_bias[idesign] <- varcor_mean[1,1] - sigma_b0_true^2
    #     varcomp_result_lmm$sigma_b0_sd[idesign] <- varcor_sd[1,1]
    #     varcomp_result_lmm$sigma_b2_bias[idesign] <- varcor_mean[2,2] - sigma_b2_true^2
    #     varcomp_result_lmm$sigma_b2_sd[idesign] <- varcor_sd[2,2]
    # } else if (dgm_type == 2) {
    #     varcomp_result_lmm$sigma_b0_bias[idesign] <- varcor_mean[1,1] - sigma_b0_true^2
    #     varcomp_result_lmm$sigma_b0_sd[idesign] <- varcor_sd[1,1]
    #     varcomp_result_lmm$sigma_b1_bias[idesign] <- varcor_mean[2,2] - sigma_b1_true^2
    #     varcomp_result_lmm$sigma_b1_sd[idesign] <- varcor_sd[2,2]
    #     varcomp_result_lmm$sigma_b2_bias[idesign] <- varcor_mean[3,3] - sigma_b2_true^2
    #     varcomp_result_lmm$sigma_b2_sd[idesign] <- varcor_sd[3,3]
    #     varcomp_result_lmm$sigma_b3_bias[idesign] <- varcor_mean[4,4] - sigma_b3_true^2
    #     varcomp_result_lmm$sigma_b3_sd[idesign] <- varcor_sd[4,4]
    # }
    # varcomp_result_lmm$sigma_eps_bias[idesign] <- mean(sigma_eps) - sigma_eps_true^2
    # varcomp_result_lmm$sigma_eps_sd[idesign] <- sd(sigma_eps)
    
    ### GEE with independence
    
    gee_ind_alpha_0 <- sapply(result_gee_ind, function(l) l$coef["(Intercept)", "Estimate"])
    gee_ind_alpha_0_sd <- sapply(result_gee_ind, function(l) l$coef["(Intercept)", "Std.err"])
    gee_ind_alpha_1 <- sapply(result_gee_ind, function(l) l$coef["X", "Estimate"])
    gee_ind_alpha_1_sd <- sapply(result_gee_ind, function(l) l$coef["X", "Std.err"])
    gee_ind_beta_0 <- sapply(result_gee_ind, function(l) l$coef["A", "Estimate"])
    gee_ind_beta_0_sd <- sapply(result_gee_ind, function(l) l$coef["A", "Std.err"])
    gee_ind_beta_1 <- sapply(result_gee_ind, function(l) l$coef["X:A", "Estimate"])
    gee_ind_beta_1_sd <- sapply(result_gee_ind, function(l) l$coef["X:A", "Std.err"])
    
    design$gee_ind_alpha_0_bias[idesign] <- mean(gee_ind_alpha_0) - alpha_0_true
    design$gee_ind_alpha_0_sd[idesign] <- sd(gee_ind_alpha_0)
    design$gee_ind_alpha_1_bias[idesign] <- mean(gee_ind_alpha_1) - alpha_1_true
    design$gee_ind_alpha_1_sd[idesign] <- sd(gee_ind_alpha_1)
    design$gee_ind_beta_0_bias[idesign] <- mean(gee_ind_beta_0) - beta_0_true
    design$gee_ind_beta_0_sd[idesign] <- sd(gee_ind_beta_0)
    design$gee_ind_beta_1_bias[idesign] <- mean(gee_ind_beta_1) - beta_1_true
    design$gee_ind_beta_1_sd[idesign] <- sd(gee_ind_beta_1)
    
    ### GEE with exchangeable
    
    gee_ex_alpha_0 <- sapply(result_gee_ex, function(l) l$coef["(Intercept)", "Estimate"])
    gee_ex_alpha_0_sd <- sapply(result_gee_ex, function(l) l$coef["(Intercept)", "Std.err"])
    gee_ex_alpha_1 <- sapply(result_gee_ex, function(l) l$coef["X", "Estimate"])
    gee_ex_alpha_1_sd <- sapply(result_gee_ex, function(l) l$coef["X", "Std.err"])
    gee_ex_beta_0 <- sapply(result_gee_ex, function(l) l$coef["A", "Estimate"])
    gee_ex_beta_0_sd <- sapply(result_gee_ex, function(l) l$coef["A", "Std.err"])
    gee_ex_beta_1 <- sapply(result_gee_ex, function(l) l$coef["X:A", "Estimate"])
    gee_ex_beta_1_sd <- sapply(result_gee_ex, function(l) l$coef["X:A", "Std.err"])
    
    design$gee_ex_alpha_0_bias[idesign] <- mean(gee_ex_alpha_0) - alpha_0_true
    design$gee_ex_alpha_0_sd[idesign] <- sd(gee_ex_alpha_0)
    design$gee_ex_alpha_1_bias[idesign] <- mean(gee_ex_alpha_1) - alpha_1_true
    design$gee_ex_alpha_1_sd[idesign] <- sd(gee_ex_alpha_1)
    design$gee_ex_beta_0_bias[idesign] <- mean(gee_ex_beta_0) - beta_0_true
    design$gee_ex_beta_0_sd[idesign] <- sd(gee_ex_beta_0)
    design$gee_ex_beta_1_bias[idesign] <- mean(gee_ex_beta_1) - beta_1_true
    design$gee_ex_beta_1_sd[idesign] <- sd(gee_ex_beta_1)
    
    ### GEE with AR(1)
    
    gee_ar1_alpha_0 <- sapply(result_gee_ar1, function(l) l$coef["(Intercept)", "Estimate"])
    gee_ar1_alpha_0_sd <- sapply(result_gee_ar1, function(l) l$coef["(Intercept)", "Std.err"])
    gee_ar1_alpha_1 <- sapply(result_gee_ar1, function(l) l$coef["X", "Estimate"])
    gee_ar1_alpha_1_sd <- sapply(result_gee_ar1, function(l) l$coef["X", "Std.err"])
    gee_ar1_beta_0 <- sapply(result_gee_ar1, function(l) l$coef["A", "Estimate"])
    gee_ar1_beta_0_sd <- sapply(result_gee_ar1, function(l) l$coef["A", "Std.err"])
    gee_ar1_beta_1 <- sapply(result_gee_ar1, function(l) l$coef["X:A", "Estimate"])
    gee_ar1_beta_1_sd <- sapply(result_gee_ar1, function(l) l$coef["X:A", "Std.err"])
    
    design$gee_ar1_alpha_0_bias[idesign] <- mean(gee_ar1_alpha_0) - alpha_0_true
    design$gee_ar1_alpha_0_sd[idesign] <- sd(gee_ar1_alpha_0)
    design$gee_ar1_alpha_1_bias[idesign] <- mean(gee_ar1_alpha_1) - alpha_1_true
    design$gee_ar1_alpha_1_sd[idesign] <- sd(gee_ar1_alpha_1)
    design$gee_ar1_beta_0_bias[idesign] <- mean(gee_ar1_beta_0) - beta_0_true
    design$gee_ar1_beta_0_sd[idesign] <- sd(gee_ar1_beta_0)
    design$gee_ar1_beta_1_bias[idesign] <- mean(gee_ar1_beta_1) - beta_1_true
    design$gee_ar1_beta_1_sd[idesign] <- sd(gee_ar1_beta_1)
    
}

### Post-processing of Results ----------------------------------------------

# reorder the columns dmg_type, total_T, sample_size to be first
design2 <- design[, c(3, 2, 1, 4:ncol(design))]
# change dmg_type to GM, total_T to T, sample_size to N
colnames(design2) <- c("GM", "T", "N", colnames(design2)[4:ncol(design2)])
saveRDS(design2, "qian2020/output/results_all.RDS")

# exclude any column that doesn't include "beta_0" and the first three columns
design3a <- design2[, c(1:3, grep("beta_0", colnames(design2)))]
design3a <- design3a[, c("GM", "T", "N", "mlm_beta_0_bias", "mlm_beta_0_sd", "gee_ex_beta_0_bias", "gee_ex_beta_0_sd", 
                         "gee_ar1_beta_0_bias", "gee_ar1_beta_0_sd", "gee_ind_beta_0_bias", "gee_ind_beta_0_sd")]
saveRDS(design3a, "qian2020/output/results_beta0_bias_sd.RDS")

# exclude any column that doesn't include "beta_0_bias" and the first three columns
design3b <- design2[, c(1:3, grep("beta_0_bias", colnames(design2)))]
design4 <- design3b[, c("GM", "T", "N", "mlm_beta_0_bias", "gee_ex_beta_0_bias", "gee_ar1_beta_0_bias", "gee_ind_beta_0_bias")]
saveRDS(design4, "qian2020/output/results_beta0_bias.RDS")

# varcomp_result_lmm <- varcomp_result_lmm[order(varcomp_result_lmm$dgm_type, varcomp_result_lmm$total_T, varcomp_result_lmm$sample_size), ]

# make LaTeX tables --------------------------------------------------------

library(xtable)

# reorder_col <- c(3, 2, 1, 9, 8)
# design_reorder <- design[, reorder_col]
# colnames(design_reorder) <- c("GM", "T", "N", "beta0_bias", "beta0_sd")
print(xtable(design4, digits = c(0, 0, 0, 0, rep(3, 4)), caption = "Results for beta0 bias, 100 replications", label = "tab:beta0_bias"),
      include.rownames = FALSE, hline.after = c(-1, 0, seq(from = 6, to = nrow(design4), by = 6)))

print(xtable(design3a, digits = c(0, 0, 0, 0, rep(3, 8)), caption = "Results for beta0 bias with Standarddeviation, 100 replications", label = "tab:beta0_bias"),
      include.rownames = FALSE, hline.after = c(-1, 0, seq(from = 6, to = nrow(design4), by = 6)))


# reorder_col <- c(1:3, 9:8, 7:6, 13:12, 11:10, 5:4)
# varcomp_result_lmm_reorder <- varcomp_result_lmm[, reorder_col]
# print(xtable(varcomp_result_lmm_reorder, digits = c(0, 0, 0, 0, rep(3, 10))),
#       include.rownames = FALSE, hline.after = c(-1, 0, seq(from = 4, to = nrow(varcomp_result_lmm_reorder), by = 4)))


