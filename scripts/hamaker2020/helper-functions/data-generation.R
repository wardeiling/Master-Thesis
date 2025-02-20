library(tidyverse)

glmm_data_generation <- function(N_total, T_total, predictor.type, outcome.type, 
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

### FUNCTION TESTING ###

if(testing == TRUE){
  # generate continuous data
  data_cont <- glmm_data_generation(N_total = 5000, T_total = 20, predictor.type = "continuous", outcome.type = "continuous",
                                    sdX.within = sqrt(1), sdX.between = sqrt(4), g.00 = 0, g.01 = 2, sd.u0 = 1,
                                    g.10 = 1, sd.u1 = 0, sd.e = 1)
  
  ggplot(data_cont, aes(x = Y)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of Y")
  ggplot(data_cont, aes(x = X)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of X")
  
  # generate binary data
  data_binary <- glmm_data_generation(N_total = 5000, T_total = 20, predictor.type = "binary", outcome.type = "binary",
                                      sdX.within = sqrt(1), sdX.between = sqrt(4), g.00 = 0, g.01 = 2, sd.u0 = 1,
                                      g.10 = 1, sd.u1 = 0, sd.e = 1)
  
  ggplot(data_binary, aes(x = Y)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of Y")
  ggplot(data_binary, aes(x = X)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of X")
  
  ### Make Evaluation Plots
  # 1. histograms with density overlay for evaluation of probability distributions
  
  ggplot(data_binary, aes(x = p_X_mean)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of p_X_mean")
  
  ggplot(data_binary, aes(x = p_y)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of p_y")
  
  # when sdX.between is large (e.g, sqrt(4)) the distribution of probabilities 
  # of X and Y are U/bowl shaped (higher density at the edges)
  # when sdX.between is small (e.g., sqrt(0.4)) the distribution of probabilities 
  # of X and Y are inverted U shaped (higher density in the center)
  # when sdX.between is moderate (e.g., sqrt(2)) the distribution of probabilities
  # is more uniform but inverted U shaped for X and slightly U shaped for Y.
  
  # 2. scatterplots for evaluation of relationship between probabilities and originating variables
  
  ggplot(data_binary, aes(x = X_mean, y = p_X_mean)) + geom_point() + labs(title = "Scatterplot of X_mean and p_X_mean")
  # As expected, we see a a sigmoid curve representing the logit relationship between X_mean 
  # and p_X_mean. It is more pronounced (different from linear) when sdX.between is large.
  
  ggplot(data_binary, aes(x = eta, y = p_y)) + geom_point() + labs(title = "Scatterplot of eta and p_y")
  # As expected again, we see a sigmoid curve representing the logit relationship between eta
  # and p_y. It is more pronounced (different from linear) when sdX.between and sd.u0 are large.
}

