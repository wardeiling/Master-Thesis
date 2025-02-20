library(tidyverse)

generate_data2 <- function(N_total, T_total, predictor.type, outcome.type, 
                          sdX.within, sdX.between, g.00, g.01, sd.u0, 
                          g.10, sd.u1, sd.e) {
  
  # Precompute Level 2 (Cluster-Level) Variables
  X.mean <- rnorm(N_total, mean = 0, sd = sdX.between)  # Cluster means
  b0 <- g.00 + g.01 * X.mean + rnorm(N_total, 0, sd.u0)  # Random intercepts
  b1 <- g.10 + rnorm(N_total, 0, sd.u1)  # Random slopes
  
  # Initialize matrices
  X <- matrix(NA, N_total, T_total)
  Y <- matrix(NA, N_total, T_total)
  
  # Loop over clusters (Level 2)
  for (j in 1:N_total) {
    # Generate Level 1 (Observation-Level) Predictor
    if (predictor.type == "continuous") {
      X.j <- rnorm(T_total, mean = X.mean[j], sd = sdX.within)
    } else if (predictor.type == "binary") {
      p.j <- plogis(X.mean[j])
      X.j <- rbinom(T_total, 1, p.j)  # Logit link
    }
    
    X[j, ] <- X.j  # Store predictor
    
    # Compute Outcome
    if (outcome.type == "continuous") {
      Y[j, ] <- b0[j] + b1[j] * (X.j - X.mean[j]) + rnorm(T_total, 0, sd.e)
    } else if (outcome.type == "binary") {
      eta <- b0[j] + b1[j] * (X.j - X.mean[j])
      p.y <- plogis(eta)
      Y[j, ] <- rbinom(T_total, 1, p.y)
    }
  }
  
  # Convert to long-format data frame
  data <- data.frame(
    Y = c(t(Y)), 
    X = c(t(X)), 
    Cluster = rep(1:N_total, each = T_total), 
    Time = rep(1:T_total, N_total)
  )
  
  return(data)
}


# DATA GENERATION FUNCTION
generate_data <- function(N, n, predictor.type, outcome.type, 
                          sdX.within, sdX.between, g.00, g.01, sd.u0, 
                          g.10, sd.u1, sd.e) {
  
  # Initialize matrices
  Y <- matrix(NA, N, n)    
  X <- matrix(NA, N, n)  
  
  # Generate level 2 elements
  X.mean <- rnorm(N, mean = 0, sd = sdX.between) # between-subject variance
  u0 <- rnorm(N, mean = 0, sd = sd.u0) 
  b0 <- g.00 + g.01 * X.mean + u0 # level 2 intercept
  u1 <- rnorm(N, mean = 0, sd = sd.u1)
  b1 <- g.10 + u1 # level 2 slope
  
  # loop over all individuals j
  for (j in 1:N) {
    
    X.mean.j <- X.mean[j]
    
    if (predictor.type == "continuous") {
      X.j <- rnorm(n, mean = X.mean.j, sd = sdX.within)
    } else if (predictor.type == "binary") {
      p <- plogis(X.mean.j)
      X.j <- rbinom(n, 1, p)
    }
    
    X[j, ] <- X.j
    b1.j <- b1[j]
    b0.j <- b0[j]
    
    if (outcome.type == "continuous") {
      Y[j, ] <- b0.j + b1.j * (X.j - X.mean.j) + rnorm(n, 0, sd.e)
    } else if (outcome.type == "binary") {
      eta <- b0.j + b1.j * (X.j - X.mean.j)
      p.y <- plogis(eta)
      Y[j, ] <- rbinom(n, 1, p.y)
    }
  }
  
  data <- data.frame(
    Y = c(t(Y)), 
    X = c(t(X)), 
    Cluster = rep(1:N, each = n), 
    Time = rep(1:n, N)
  )
  
  return(data)
}

generate_data3 <- function(N_total, T_total, predictor.type, outcome.type, 
                          sdX.within, sdX.between, g.00, g.01, sd.u0, 
                          g.10, sd.u1, sd.e) {
  
  # Define variable names and initialize data frame
  if (outcome.type == "continuous") {
    df_names <- c("Y", "X", "Cluster", "Time", "X_mean", "b0", "b1", "eta")
  } else if (outcome.type == "binary") {
    df_names <- c("Y", "X", "Cluster", "Time", "X_mean", "b0", "b1", "eta", "p_y")
  }
  
  df_names <- c("Y", "X", "Cluster", "Time", "X_mean", "b0", "b1", "eta")
  dta <- data.frame(matrix(NA, nrow = N_total * T_total, ncol = length(df_names)))
  colnames(dta) <- df_names
  
  # Precompute Level 2 (Cluster-Level) Variables
  X.mean <- rnorm(N_total, mean = 0, sd = sdX.between)  # Cluster means
  b0 <- g.00 + g.01 * X.mean + rnorm(N_total, 0, sd.u0)  # Random intercepts
  b1 <- g.10 + rnorm(N_total, 0, sd.u1)  # Random slopes
  
  # Assign Cluster, Time, and Precomputed Level 2 Variables
  dta$Cluster <- rep(1:N_total, each = T_total)
  dta$Time <- rep(1:T_total, N_total)
  dta$X_mean <- rep(X.mean, each = T_total)  # Store cluster means
  dta$b0 <- rep(b0, each = T_total)  # Store intercepts
  dta$b1 <- rep(b1, each = T_total)  # Store slopes
  
  # Loop over time points (Level 1)
  for (t in 1:T_total) {
    row_index <- seq(from = t, by = T_total, length = N_total)  # Efficient indexing
    
    # Generate Level 1 (Observation-Level) Predictor
    if (predictor.type == "continuous") {
      dta$X[row_index] <- rnorm(N_total, mean = X.mean, sd = sdX.within)
    } else if (predictor.type == "binary") {
      p <- plogis(X.mean.j)
      X.j <- rbinom(n, 1, p)
      dta$X[row_index] <- rbinom(N_total, 1, plogis(X.mean))  # Logit link
    }
    
    # Compute Linear Predictor (eta) Before Applying Link Function
    dta$eta[row_index] <- dta$b0[row_index] + dta$b1[row_index] * (dta$X[row_index] - dta$X_mean[row_index])
    
    # Compute Outcome
    if (outcome.type == "continuous") {
      dta$Y[row_index] <- dta$eta[row_index] + rnorm(N_total, 0, sd.e)
    } else if (outcome.type == "binary") {
      dta$Y[row_index] <- rbinom(N_total, 1, plogis(dta$eta[row_index]))
    }
  }
  
  return(dta)
}


# test functions
data1 <- generate_data(N = 5000, # number of clusters
                       n = 20, # number of observations within a cluster, originally set to 4.
                       predictor.type = "continuous", # type of predictor
                       outcome.type = "continuous", # type of outcome
                       
                       # PREDICTOR
                       sdX.within = sqrt(1),		# within-person variance 
                       sdX.between = sqrt(4),	# between-person variance
                       # if set to zero, the marginal effect will become approximately equal to the conditional effect.
                       
                       # INTERCEPT LEVEL 2
                       g.00 = 0,			# Grand intercept
                       g.01 = 2,			# between-cluster slope
                       sd.u0 = 1,			# SD of residuals intercept at level 2
                       
                       # SLOPE LEVEL 2
                       g.10 = 1,			# fixed within-cluster slope
                       sd.u1 = 0,			# SD of within-cluster slope at level 2
                       
                       # RESIDUALS AT LEVEL 1
                       sd.e = 1			# residual SD at level 1 (only used when outcome.type is continuous)
)

data2 <- generate_data2(N = 5000, # number of clusters
                       n = 20, # number of observations within a cluster, originally set to 4.
                       predictor.type = "continuous", # type of predictor
                       outcome.type = "continuous", # type of outcome
                       
                       # PREDICTOR
                       sdX.within = sqrt(1),		# within-person variance 
                       sdX.between = sqrt(4),	# between-person variance
                       # if set to zero, the marginal effect will become approximately equal to the conditional effect.
                       
                       # INTERCEPT LEVEL 2
                       g.00 = 0,			# Grand intercept
                       g.01 = 2,			# between-cluster slope
                       sd.u0 = 1,			# SD of residuals intercept at level 2
                       
                       # SLOPE LEVEL 2
                       g.10 = 1,			# fixed within-cluster slope
                       sd.u1 = 0,			# SD of within-cluster slope at level 2
                       
                       # RESIDUALS AT LEVEL 1
                       sd.e = 1			# residual SD at level 1 (only used when outcome.type is continuous)
)

summary(data1)
summary(data2)
