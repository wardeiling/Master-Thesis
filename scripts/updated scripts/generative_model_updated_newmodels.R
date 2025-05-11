# 2018.10.22 Tianchen Qian

# generative model:
# Y_t+1 = alpha_0 + alpha_1 X_t + b_0i + b_1i X_t + A_t (beta_0 + beta_1 X_t + b_2i + b_3i X_t) + epsilon_it

### Naming Scheme 1 ###

# dmg_type with an "a" suffix means that there are no random slopes
# dgm_type with a "b" suffix means that there are no random slopes nor interactions between the treatment and covariate
# dgm_type == 3c: GM3, except higher value for the random slope b_i2 variance (treatment effect)
# dgm_type == 3d: GM3, except no interaction effect beta_1
# dgm_type == 3e: GM3, except higher value for interaction beta_1
# dgm_type == 3f: GM3, except higher value for fixed slope alpha_1
# dgm_type == 3g: GM3, except lower value for fixed slope alpha_1
# dgm_type == 3h: GM3, except no fixed slope alpha_1
# dgm_type == 3i: GM3, except no random slope b_i2 nor fixed slope alpha_1
# dgm_type == 3j: GM3, except no interaction effect beta_1 nor fixed slope alpha_1

dgm_with_treatment <- function(sample_size, total_T, dgm_type) {
    
    # dgm_type is in c(1,2,3,4)
    stopifnot(dgm_type %in% c(1,2,3,"1a", "2a", "3a", "1b", "2b", "3b", "3c", "3d", "3e", "3f", "3g", "3h", "3i", "3j"))
    
    # dgm_type = 1 or 3
    alpha_0 <- - 2 # overall intercept, was originally -1 in the code but is -2 in the paper
    alpha_1 <- - 0.3 # overall slope of covariate, same as in the paper
    beta_0 <- 1 # treatment effect/slope, was originallly 0.5 in the code but is 1 in the paper
    beta_1 <- 0.3 # interaction effect, was originally 0.1 in the code but is 0.3 in the paper
    sigma_b0 <- 2 # sqrt(4)
    sigma_b1 <- 0 
    sigma_b2 <- 1 # sqrt(1)
    sigma_b3 <- 0
    sigma_eps <- 1 # sqrt(1)
    
    # activate the random slopes for model 2
    if (dgm_type == 2) {
      sigma_b1 <- sigma_b3 <- 0.5 # sqrt(0.25)
      # set the random slopes to 0 for "a" models
    } else if (dgm_type %in% c("1a", "2a", "3a")) {
      sigma_b2 <- 0
      # set the interactions to 0 for "b" models
    } else if (dgm_type %in% c("1b", "2b", "3b")) {
      sigma_b2 <- 0
      beta_1 <- 0
      # increase the variance of random slope b_2i for model 3c
    } else if (dgm_type == "3c") {
      sigma_b2 <- 3
    } else if (dgm_type == "3d") {
      beta_1 <- 0
    } else if (dgm_type == "3e") {
      beta_1 <- 0.6
    } else if (dgm_type == "3f") {
      alpha_1 <- -0.6
    } else if (dgm_type == "3g") {
      alpha_1 <- -0.15
    } else if (dgm_type == "3h") {
      alpha_1 <- 0
    } else if (dgm_type == "3i") {
      sigma_b2 <- 0
      alpha_1 <- 0
    } else if (dgm_type == "3j") {
      beta_1 <- 0
      alpha_1 <- 0
    }

    prob_a <- 0.5
    
    df_names <- c("userid", "day", "X", "prob_A", "A", "Y", "b0", "b1", "b2", "b3", "eps", "delta")
    
    dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$userid <- rep(1:sample_size, each = total_T)
    dta$day <- rep(1:total_T, times = sample_size)
    
    # uncorrelated random effects
    b_0i <- rnorm(sample_size, mean = 0, sd = sigma_b0)
    b_1i <- rnorm(sample_size, mean = 0, sd = sigma_b1)
    b_2i <- rnorm(sample_size, mean = 0, sd = sigma_b2)
    b_3i <- rnorm(sample_size, mean = 0, sd = sigma_b3)
    
    # b_1i[b_1i > 2] <- 2
    # b_1i[b_1i < -2] <- -2
    # 
    # b_3i[b_3i > 2] <- 2
    # b_3i[b_3i < -2] <- -2
    
    for (t in 1:total_T) {
        # row index for the rows corresponding to day t for every subject
        row_index <- seq(from = t, by = total_T, length = sample_size)
        
        if (dgm_type %in% c(1, "1a", "1b")) {
            if (t == 1) {
                dta$X[row_index] <- rnorm(sample_size)
            } else {
                dta$X[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size)
            }
            dta$prob_A[row_index] <- rep(prob_a, sample_size)
            
        } else if (dgm_type %in% c(2, "2a", "2b")) {
            if (t == 1) {
                dta$X[row_index] <- rnorm(sample_size)
            } else {
                dta$X[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size)
            }
            dta$prob_A[row_index] <- ifelse(dta$X[row_index] > - 1.27, 0.7, 0.3)
            
        } else if (dgm_type %in% c(3, "3a", "3b", "3c", "3d", "3e", "3f", "3g", "3h", "3i", "3j")) {
            if (t == 1) {
                dta$X[row_index] <- rnorm(sample_size, mean = b_0i) # X involves b_i!!
            } else {
                dta$X[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size, mean = b_0i) # X involves b_i!!
            }
            dta$prob_A[row_index] <- rep(prob_a, sample_size)
        }
        
        dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index])
        
        dta$eps[row_index] <- rnorm(sample_size, mean = 0, sd = sigma_eps)
        
        dta$delta[row_index] <- beta_0 + beta_1 * dta$X[row_index] + b_2i + b_3i * dta$X[row_index]
        
        dta$Y[row_index] <- alpha_0 + alpha_1 * dta$X[row_index] + b_0i + b_1i * dta$X[row_index] + 
            dta$A[row_index] * dta$delta[row_index] + dta$eps[row_index]
        
        dta$b0[row_index] <- b_0i
        dta$b1[row_index] <- b_1i
        dta$b2[row_index] <- b_2i
        dta$b3[row_index] <- b_3i
        
        row_index_lag1 <- row_index
    }
    
    return(dta)
}


##### example: use of lmer() #####
if( 0 ){
  
    library(lme4)
    library(geepack)
    library(gee)
  
    dta <- dgm_with_treatment(sample_size = 200, total_T = 30, dgm_type = "1")
    # hist(dta$Y)
    summary(dta)
    # dta$A <- dta$A - dta$prob_A # action centering doesn't matter when prob_A is constant
    mlm_fit <- lmer(Y ~ A + A:X + (1 + A | userid), data = dta)
    summary(mlm_fit)
    (mlm1 <- summary(lmer(Y ~ X + A + (1  | userid), data = dta))$coefficients)
    (mlm2 <- lmer(Y ~ X * A + (X * A | userid), data = dta))
    isSingular(mlm2, tol = 1e-5)
    (gee_ex <- summary(geepack::geeglm(Y ~ A + A:X, id = userid, data = dta, family = gaussian, corstr = "exchangeable"))$coefficients)
    (gee_in <- summary(geepack::geeglm(Y ~ X * A, id = userid, data = dta, family = gaussian, corstr = "independence"))$coefficients)
    (gee_ar <- summary(geepack::geeglm(Y ~ X * A, id = userid, data = dta, family = gaussian, corstr = "ar1"))$coefficients)
    (gee_in2 <- summary(gee(Y ~ X * A, id = userid, data = dta, family = gaussian, corstr = "independence"))$coefficients)
    (glm <- summary(glm(Y ~ X * A, data = dta, family = gaussian))$coefficients)
    # gee_ar1 <- gee::gee(Y ~ X * A, id = userid, data = dta, family = gaussian, corstr = "AR1")
    
    col <- c("Estimate", "Std.err", "Wald", "Pr(>|W|)")
    row <- c("(Intercept)", "X", "A", "X:A")
    col == colnames(gee_ex)
    row == rownames(gee_ex)
    
    fit <- lmer(Y ~ X * A + (X * A | userid), data = dta)
    fit
    
    summary(fit)$coefficients
    
    attr(summary(fit)$varcor$userid, "stddev") # estimated standard deviation of random effect
    
    ### Check of output ###
    
    library(geepack)
    
    # Simulate data with perfect multicollinearity
    set.seed(123)
    n <- 50
    id <- rep(1:10, each = 5)  # 10 clusters of size 5
    x <- rnorm(n)             # Predictor
    y <- rbinom(n, 1, prob = 0.5)  # Binary outcome
    x_collinear <- x * 2      # Perfectly collinear variable
    
    # Combine into a data frame
    data <- data.frame(id = id, y = y, x = x, x_collinear = x_collinear)
    
    # Fit a GEE model with collinear predictors
    tryCatch({
      model <- geeglm(y ~ x + x_collinear, family = binomial, id = id, data = data)
      summary(model)
    }, warning = function(w) {
      message("Warning occurred: ", conditionMessage(w))
    }, error = function(e) {
      message("Error occurred: ", conditionMessage(e))
    })
}



