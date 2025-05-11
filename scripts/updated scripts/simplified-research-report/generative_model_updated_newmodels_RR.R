# 2018.10.22 Tianchen Qian

# generative model:
# Y_t+1 = alpha_0 + alpha_1 X_t + b_0i + b_1i X_t + A_t (beta_0 + beta_1 X_t + b_2i + b_3i X_t) + epsilon_it

### Naming Scheme 1 ###

# dgm_type == 3a: GM3, except no random slope b_2i
# dgm_type == 3d: GM3, except no interaction effect beta_1

dgm_with_treatment <- function(sample_size, total_T, dgm_type) {
    
    # dgm_type is in c(1,2,3,4)
    stopifnot(dgm_type %in% c(1,2,3, "3a", "3d"))
    
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
    } else if (dgm_type %in% c("3a")) {
      sigma_b2 <- 0
    } else if (dgm_type == "3d") {
      beta_1 <- 0
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
        
        if (dgm_type == 1) {
            if (t == 1) {
                dta$X[row_index] <- rnorm(sample_size)
            } else {
                dta$X[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size)
            }
            dta$prob_A[row_index] <- rep(prob_a, sample_size)
            
        } else if (dgm_type == 2) {
            if (t == 1) {
                dta$X[row_index] <- rnorm(sample_size)
            } else {
                dta$X[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size)
            }
            dta$prob_A[row_index] <- ifelse(dta$X[row_index] > - 1.27, 0.7, 0.3)
            
        } else if (dgm_type %in% c(3, "3a", "3d")) {
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

