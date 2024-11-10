# generative model:
# Y_t+1 = alpha_0 + alpha_1 X_t + b_0i + A_t (beta_0 + beta_1 X_t + b_2i) + epsilon_it
# Y_{it+1} = \gamma_{00} + \beta_2 Z_{it} + u_{0i} + X_{it} (\gamma_{10} + \beta_3 Z_{it} + u_{1i}) + e_{it+1}

# Data generating models

GM1 <- function(sample_size, total_T) {
    
    # dgm_type = 1 or 3
    alpha_0 <- - 1
    alpha_1 <- - 0.3
    beta_0 <- 0.5
    beta_1 <- 0.1
    sigma_b0 <- 2
    sigma_b2 <- 1
    sigma_eps <- 1
    
    prob_a <- 0.5
    
    df_names <- c("userid", "day", "X", "prob_A", "A", "Y", "b0", "b1", "b2", "b3", "eps", "delta")
    
    dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$userid <- rep(1:sample_size, each = total_T)
    dta$day <- rep(1:total_T, times = sample_size)
    
    # uncorrelated random effects
    b_0i <- rnorm(sample_size, mean = 0, sd = sigma_b0)
    b_2i <- rnorm(sample_size, mean = 0, sd = sigma_b2)
    
    for (t in 1:total_T) {
        # row index for the rows corresponding to day t for every subject
        row_index <- seq(from = t, by = total_T, length = sample_size)
        
        if (t == 1) {
                dta$X[row_index] <- rnorm(sample_size)
            } else {
                dta$X[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size)
            }
            
        dta$prob_A[row_index] <- rep(prob_a, sample_size)
        
        dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index])
        
        dta$eps[row_index] <- rnorm(sample_size, mean = 0, sd = sigma_eps)
        
        dta$delta[row_index] <- beta_0 + beta_1 * dta$X[row_index] + b_2i
        
        dta$Y[row_index] <- alpha_0 + alpha_1 * dta$X[row_index] + b_0i + 
                            dta$A[row_index] * dta$delta[row_index] + dta$eps[row_index]
        
        dta$b0[row_index] <- b_0i
        dta$b2[row_index] <- b_2i
        
        row_index_lag1 <- row_index
    }
    
    return(dta)
}