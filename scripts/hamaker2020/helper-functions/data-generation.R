library(tidyverse)

generate_data <- function(N_total, T_total, predictor.type, outcome.type, 
                          sdX.within, sdX.between, g.00, g.01, sd.u0, 
                          g.10, sd.u1, sd.e) {
  
  # Precompute Level 2 (Cluster-Level) Variables
  X.mean <- rnorm(N_total, mean = 0, sd = sdX.between)  # Cluster means
  p_X_mean <- if (predictor.type == "binary") plogis(X.mean) else rep(NA, N_total)
  b0 <- g.00 + g.01 * X.mean + rnorm(N_total, 0, sd.u0)  # Random intercepts
  b1 <- g.10 + rnorm(N_total, 0, sd.u1)  # Random slopes
  
  # Initialize storage for long-format data
  dta <- data.frame(
    Cluster = rep(1:N_total, each = T_total),
    Time = rep(1:T_total, times = N_total),
    X_mean = rep(X.mean, each = T_total),
    b0 = rep(b0, each = T_total),
    b1 = rep(b1, each = T_total),
    p_X_mean = rep(p_X_mean, each = T_total),
    X = NA, eta = NA, Y = NA, p_y = NA  # Preallocate necessary columns
  )
  
  # Generate Level 1 data
  for (j in 1:N_total) {
    idx <- which(dta$Cluster == j)  # Rows belonging to cluster j
    
    # Generate predictor
    if (predictor.type == "continuous") {
      X_j <- rnorm(T_total, mean = X.mean[j], sd = sdX.within)
    } else {
      X_j <- rbinom(T_total, 1, p_X_mean[j])
    }
    
    # Compute eta and outcome
    eta_j <- b0[j] + b1[j] * (X_j - X.mean[j])
    
    if (outcome.type == "continuous") {
      Y_j <- eta_j + rnorm(T_total, mean = 0, sd = sd.e)
      p_y_j <- NA  # No probability needed for continuous outcome
    } else {
      p_y_j <- plogis(eta_j)
      Y_j <- rbinom(T_total, 1, p_y_j)
    }
    
    # Store values
    dta$X[idx] <- X_j
    dta$eta[idx] <- eta_j
    dta$Y[idx] <- Y_j
    dta$p_y[idx] <- p_y_j
  }
  
  # Compute cluster-level means and centering
  dta$X_mean_est <- ave(dta$X, dta$Cluster, FUN = mean)
  dta$X_cent <- dta$X - dta$X_mean
  dta$X_cent_est <- dta$X - dta$X_mean_est
  
  return(dta)
}

# test function
data <- generate_data(N_total = 5000, # number of clusters
                      T_total = 20, # number of observations within a cluster, originally set to 4.
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


summary(data)

# generate_data2 <- function(N_total, T_total, predictor.type, outcome.type, 
#                           sdX.within, sdX.between, g.00, g.01, sd.u0, 
#                           g.10, sd.u1, sd.e) {
#   
#   # Precompute Level 2 (Cluster-Level) Variables
#   X.mean <- rnorm(N_total, mean = 0, sd = sdX.between)  # Cluster means
#   b0 <- g.00 + g.01 * X.mean + rnorm(N_total, 0, sd.u0)  # Random intercepts
#   b1 <- g.10 + rnorm(N_total, 0, sd.u1)  # Random slopes
#   
#   # Initialize matrices
#   X <- matrix(NA, N_total, T_total)
#   Y <- matrix(NA, N_total, T_total)
#   
#   # Loop over clusters (Level 2)
#   for (j in 1:N_total) {
#     # Generate Level 1 (Observation-Level) Predictor
#     if (predictor.type == "continuous") {
#       X.j <- rnorm(T_total, mean = X.mean[j], sd = sdX.within)
#     } else if (predictor.type == "binary") {
#       p.j <- plogis(X.mean[j])
#       X.j <- rbinom(T_total, 1, p.j)  # Logit link
#     }
#     
#     X[j, ] <- X.j  # Store predictor
#     
#     # Compute Outcome
#     if (outcome.type == "continuous") {
#       Y[j, ] <- b0[j] + b1[j] * (X.j - X.mean[j]) + rnorm(T_total, 0, sd.e)
#     } else if (outcome.type == "binary") {
#       eta <- b0[j] + b1[j] * (X.j - X.mean[j])
#       p.y <- plogis(eta)
#       Y[j, ] <- rbinom(T_total, 1, p.y)
#     }
#   }
#   
#   # Convert to long-format data frame
#   data <- data.frame(
#     Y = c(t(Y)), 
#     X = c(t(X)), 
#     Cluster = rep(1:N_total, each = T_total), 
#     Time = rep(1:T_total, N_total)
#   )
#   
#   data$X_mean <- rep(X.mean, each = T_total)
#   data$X_mean_est <- ave(data$X, data$Cluster, FUN = mean)
#   data$X_cent <- data$X - data$X_mean
#   data$X_cent_est <- data$X - data$X_mean_est
#   
#   return(data)
# }

generate_data2 <- function(N_total, T_total, predictor.type, outcome.type, 
                           sdX.within, sdX.between, g.00, g.01, sd.u0, 
                           g.10, sd.u1, sd.e) {
  
  # Precompute Level 2 (Cluster-Level) Variables
  X.mean <- rnorm(N_total, mean = 0, sd = sdX.between)  # Cluster means
  p_X_mean <- plogis(X.mean)  # Probability of X when binary
  b0 <- g.00 + g.01 * X.mean + rnorm(N_total, 0, sd.u0)  # Random intercepts
  b1 <- g.10 + rnorm(N_total, 0, sd.u1)  # Random slopes
  
  # Initialize matrices
  X <- matrix(NA, N_total, T_total)
  Y <- matrix(NA, N_total, T_total)
  eta <- matrix(NA, N_total, T_total)
  p_y <- matrix(NA, N_total, T_total)
  
  # Loop over clusters (Level 2)
  for (j in 1:N_total) {
    # Generate Level 1 (Observation-Level) Predictor
    if (predictor.type == "continuous") {
      X.j <- rnorm(T_total, mean = X.mean[j], sd = sdX.within)
    } else if (predictor.type == "binary") {
      X.j <- rbinom(T_total, 1, p_X_mean[j])  # Logit link
    }
    
    X[j, ] <- X.j  # Store predictor
    
    # Compute Linear Predictor (eta)
    eta[j, ] <- b0[j] + b1[j] * (X.j - X.mean[j])
    
    # Compute Outcome
    if (outcome.type == "continuous") {
      Y[j, ] <- eta[j, ] + rnorm(T_total, 0, sd.e)
    } else if (outcome.type == "binary") {
      p_y[j, ] <- plogis(eta[j, ])
      Y[j, ] <- rbinom(T_total, 1, p_y[j, ])
    }
  }
  
  # Convert to long-format data frame
  data <- data.frame(
    Y = c(t(Y)), 
    X = c(t(X)), 
    eta = c(t(eta)), 
    Cluster = rep(1:N_total, each = T_total), 
    Time = rep(1:T_total, N_total),
    b0 = rep(b0, each = T_total),
    b1 = rep(b1, each = T_total),
    X_mean = rep(X.mean, each = T_total),
    p_X_mean = rep(p_X_mean, each = T_total)
  )
  
  # Post-processing
  data$X_mean_est <- ave(data$X, data$Cluster, FUN = mean)
  data$X_cent <- data$X - data$X_mean
  data$X_cent_est <- data$X - data$X_mean_est
  
  if (outcome.type == "binary") {
    data$p_y <- c(t(p_y))  # Store p_y if binary outcome
  }
  
  return(data)
}

# Efficient data generation function
generate_data3 <- function(N_total, T_total, predictor.type, outcome.type, 
                          sdX.within, sdX.between, g.00, g.01, sd.u0, 
                          g.10, sd.u1, sd.e) {
  
  # Define variable names and initialize data frame
  df_names <- c("Cluster", "Time", "Y", "X", "X_cent", "X_cent_est" ,"X_mean", "X_mean_est", "b0", "b1", "eta", "p.X_mean", "p.y")
  dta <- data.frame(matrix(NA, nrow = N_total * T_total, ncol = length(df_names)))
  colnames(dta) <- df_names
  
  # Precompute Level 2 (Cluster-Level) Variables
  X.mean <- rnorm(N_total, mean = 0, sd = sdX.between)  # Cluster means
  b0 <- g.00 + g.01 * X.mean + rnorm(N_total, 0, sd.u0)  # Random intercepts
  b1 <- g.10 + rnorm(N_total, 0, sd.u1)  # Random slopes
  
  # Assign Cluster, Time, and Precomputed Level 2 Variables
  dta$Cluster <- rep(1:N_total, each = T_total)
  dta$Time <- rep(1:T_total, times = N_total)
  dta$X_mean <- rep(X.mean, each = T_total)  # Store cluster means
  dta$b0 <- rep(b0, each = T_total)  # Store intercepts
  dta$b1 <- rep(b1, each = T_total)  # Store slopes
  
  if (predictor.type == "binary"){
    dta$p_X_mean <- plogis(dta$X_mean)  # Store probabilities
  }
  
  # Loop over time points (Level 1)
  for (t in 1:T_total) {
    # row_index <- seq(from = t, by = T_total, length = N_total)  # Efficient indexing
    row_index <- which(dta$Time == t)
    
    # Generate Level 1 (Observation-Level) Predictor
    if (predictor.type == "continuous") {
      dta$X[row_index] <- rnorm(N_total, mean = dta$X_mean[row_index], sd = sdX.within)
    } else if (predictor.type == "binary") {
      dta$X[row_index] <- rbinom(N_total, 1, dta$p_X_mean[row_index])  # Logit link
    }
    
    # Compute Linear Predictor (eta) Before Applying Link Function
    dta$eta[row_index] <- dta$b0[row_index] + dta$b1[row_index] * (dta$X[row_index] - dta$X_mean[row_index])
    
    # Compute Outcome
    if (outcome.type == "continuous") {
      dta$Y[row_index] <- dta$eta[row_index] + rnorm(N_total, mean = 0, sd = sd.e)
    } else if (outcome.type == "binary") {
      dta$p.y[row_index] <- plogis(dta$eta[row_index])
      dta$Y[row_index] <- rbinom(N_total, 1, dta$p.y[row_index])
    }
  }
  
  # Prepare data
  dta$X_mean_est <- ave(dta$X, dta$Cluster, FUN = mean)
  dta$X_cent_est <- dta$X - dta$X_mean_est
  dta$X_cent <- dta$X - dta$X_mean
  
  return(dta)
}




# test functions

data2 <- generate_data2(N_total = 5000, # number of clusters
                       T_total = 20, # number of observations within a cluster, originally set to 4.
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

summary(data2)

library(microbenchmark)

# Define function parameters
N_total <- 5000  # Number of clusters
T_total <- 10   # Number of time points per cluster
predictor.type <- "continuous"  # Can be "continuous" or "binary"
outcome.type <- "continuous"  # Can be "continuous" or "binary"
sdX.within <- 1
sdX.between <- 1
g.00 <- 0
g.01 <- 1
sd.u0 <- 1
g.10 <- 1
sd.u1 <- 1
sd.e <- 1

# Run functions and compare speed
time_comparison <- microbenchmark(
  data3 = generate_data3(N_total, T_total, predictor.type, outcome.type, 
                         sdX.within, sdX.between, g.00, g.01, sd.u0, 
                         g.10, sd.u1, sd.e),
  data2 = generate_data2(N_total, T_total, predictor.type, outcome.type, 
                         sdX.within, sdX.between, g.00, g.01, sd.u0, 
                         g.10, sd.u1, sd.e),
  times = 50
)

print(time_comparison)

set.seed(123)  # Ensure reproducibility

# Define parameters
N_total <- 200
T_total <- 20
predictor.type <- "continuous"
outcome.type <- "continuous"
sdX.within <- 1
sdX.between <- 1
g.00 <- 0
g.01 <- 1
sd.u0 <- 1
g.10 <- 0.5
sd.u1 <- 0.5
sd.e <- 1

# Storage for means and SDs
vars_to_check <- c("Y", "X", "eta", "b0", "b1", "X_mean", "X_cent")

stats_data2 <- matrix(NA, nrow = 10000, ncol = length(vars_to_check) * 2)
stats_data3 <- matrix(NA, nrow = 10000, ncol = length(vars_to_check) * 2)

colnames(stats_data2) <- colnames(stats_data3) <- c(rbind(paste0(vars_to_check, "_mean"), paste0(vars_to_check, "_sd")))

for (i in 1:10000) {
  data2 <- generate_data2(N_total, T_total, predictor.type, outcome.type, 
                          sdX.within, sdX.between, g.00, g.01, sd.u0, 
                          g.10, sd.u1, sd.e)
  
  data3 <- generate_data3(N_total, T_total, predictor.type, outcome.type, 
                          sdX.within, sdX.between, g.00, g.01, sd.u0, 
                          g.10, sd.u1, sd.e)
  
  # Compute means and SDs for each function
  stats_data2[i, ] <- unlist(lapply(data2[vars_to_check], function(x) c(mean(x), sd(x))))
  stats_data3[i, ] <- unlist(lapply(data3[vars_to_check], function(x) c(mean(x), sd(x))))
}

# Compute final mean and SD over all iterations
final_means_data2 <- colMeans(stats_data2)
final_means_data3 <- colMeans(stats_data3)
final_sds_data2 <- apply(stats_data2, 2, sd)
final_sds_data3 <- apply(stats_data3, 2, sd)

# Combine results into a comparison table
comparison_table <- data.frame(
  Variable = rep(vars_to_check, each = 2),
  Statistic = rep(c("Mean", "SD"), times = length(vars_to_check)),
  Data2 = c(final_means_data2, final_sds_data2),
  Data3 = c(final_means_data3, final_sds_data3),
  Difference = c(final_means_data2 - final_means_data3, final_sds_data2 - final_sds_data3)
)

print(comparison_table)

