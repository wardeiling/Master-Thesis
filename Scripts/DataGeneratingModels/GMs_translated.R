# Created by: Tianchen Qian
# Edited by: Ward B. Eiling

### Notation
# Original      Y_t+1 = alpha_0 + alpha_1 X_t + b_0i + b_1i X_t + A_t (beta_0 + beta_1 X_t + b_2i + b_3i X_t) + epsilon_it
# New           Y_{(t+1)i} = \beta_{00} + \beta_{10} Z_{ti} + u_{0i} + u_{1i} Z_{ti} + X_{ti} (\beta_{20} + \beta_{30} Z_{ti} 
#               + u_{2i} + u_{3i} Z_{ti}) + e_{(t+1)i}        

# generative model in Schoot et al. (2017) notation:
dgm_with_treatment_translated <- function(sample_size, total_T, dgm_type) {
  
  # dgm_type is in c(1,2,3,4)
  stopifnot(dgm_type %in% c(1,2,3,4))
  
  # Replacing parameters with equivalent Schoot et al. (2017) values
  beta_00 <- - 1
  beta_10 <- - 0.3
  beta_20 <- 1 # note: Qian script uses 0.5, but paper uses 1
  beta_30 <- 0.1
  sigma_u0 <- 2
  sigma_u1 <- 0
  sigma_u2 <- 1
  sigma_u3 <- 0
  sigma_e <- 1
  
  if (dgm_type == 2) {
    sigma_u1 <- sigma_u3 <- 0.5
  }
  if (dgm_type == 4) {
    sigma_u2 <- 0
  }
  
  prob_x <- 0.5
  
  df_names <- c("userid", "day", "Z", "prob_X", "X", "Y", "u0", "u1", "u2", "u3", "e", "delta")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$day <- rep(1:total_T, times = sample_size)
  
  # uncorrelated random effects
  u_0i <- rnorm(sample_size, mean = 0, sd = sigma_u0)
  u_1i <- rnorm(sample_size, mean = 0, sd = sigma_u1)
  u_2i <- rnorm(sample_size, mean = 0, sd = sigma_u2)
  u_3i <- rnorm(sample_size, mean = 0, sd = sigma_u3)
  
  for (t in 1:total_T) {
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    
    if (dgm_type == 1) {
      if (t == 1) {
        dta$Z[row_index] <- rnorm(sample_size)
      } else {
        dta$Z[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size)
      }
      dta$prob_X[row_index] <- rep(prob_x, sample_size)
    } else if (dgm_type == 2) {
      if (t == 1) {
        dta$Z[row_index] <- rnorm(sample_size)
      } else {
        dta$Z[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size)
      }
      dta$prob_X[row_index] <- ifelse(dta$Z[row_index] > - 1.27, 0.7, 0.3)
    } else if (dgm_type %in% c(3,4)) {
      if (t == 1) {
        dta$Z[row_index] <- rnorm(sample_size, mean = u_0i) # Z involves u_i!!
      } else {
        dta$Z[row_index] <- dta$Y[row_index_lag1] + rnorm(sample_size, mean = u_0i) # Z involves u_i!!
      }
      dta$prob_X[row_index] <- rep(prob_x, sample_size)
    }
    
    dta$X[row_index] <- rbinom(sample_size, 1, dta$prob_X[row_index])
    
    dta$e[row_index] <- rnorm(sample_size, mean = 0, sd = sigma_e)
    
    dta$delta[row_index] <- beta_20 + beta_30 * dta$Z[row_index] + u_2i + u_3i * dta$Z[row_index]
    
    dta$Y[row_index] <- beta_00 + beta_10 * dta$Z[row_index] + u_0i + u_1i * dta$Z[row_index] + 
      dta$X[row_index] * dta$delta[row_index] + dta$e[row_index]
    
    dta$u0[row_index] <- u_0i
    dta$u1[row_index] <- u_1i
    dta$u2[row_index] <- u_2i
    dta$u3[row_index] <- u_3i
    
    row_index_lag1 <- row_index
  }
  
  return(dta)
}

##### example: use of lmer() #####
set.seed(321)
sample_size <- 100000
total_T <- 10

dta <- dgm_with_treatment_translated(sample_size, total_T, dgm_type = 1)
summary(dta)
# dta$X <- dta$X - dta$prob_X # action centering doesn't matter when prob_A is constant

fit <- lmer(Y ~ Z * X + (1 + X | userid), data = dta)
summary(fit)$coefficients
fit <- lmer(Y ~ Z * X + (Z * X | userid), data = dta)
summary(fit)$coefficients

attr(summary(fit)$varcor$userid, "stddev") # estimated standard deviation of random effect
