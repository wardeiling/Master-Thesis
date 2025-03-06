library(tidyverse)

glmm_data_generation <- function(N_total, T_total, predictor.type, outcome.type, 
                                 sdX.within, sdX.between, g.00, g.01, sd.u0, g.10, sd.u1, sd.e){
  
  # Function to generate data for a generalized linear mixed model
  # ---------------------------------------------------------------------
  # when the outcome is continuous, we simulate the mixed linear model,
  # in which we add residual error to the linear predictor eta.
  # when the outcome is binary, we simulate the mixed logistic model,
  # in which we convert the linear predictor eta to probability p and
  # generate binary outcome Y from the probability.
  # ---------------------------------------------------------------------
  # when the predictor is continuous, we center it around the cluster mean.
  # when the predictor is binary, we don't center it, but we rely on the 
  # (latent) cluster mean in b0.j.
  
  # input element     description
  # ----------------  ---------------------------------------------------
  # N_total           number of clusters
  # T_total           number of observations within a cluster
  # predictor.type    type of predictor
  # outcome.type      type of outcome
  # sdX.within        within-person variance (only used when predictor.type is continuous)
  # sdX.between       between-person variance
  # g.00              grand intercept
  # g.01              between-cluster slope
  # sd.u0             SD of residuals intercept at level 2
  # g.10              fixed within-cluster slope
  # sd.u1             SD of within-cluster slope at level 2
  # sd.e              residual SD at level 1 (only used when outcome.type is continuous)
  
  # Check input
  stopifnot(predictor.type %in% c("continuous", "binary"))
  stopifnot(outcome.type %in% c("continuous", "binary"))
  
  # Precompute Level 2 (Cluster-Level) Variables
  X.mean.j <- rnorm(N_total, mean = 0, sd = sdX.between)  # Cluster means
  p.X.mean.j <- if (predictor.type == "binary") plogis(X.mean.j) else rep(NA, N_total) # convert the log-odds (logit) to probability
  b0.j <- g.00 + g.01 * X.mean.j + rnorm(N_total, 0, sd.u0)  # Random intercepts
  b1.j <- g.10 + rnorm(N_total, 0, sd.u1)  # Random slopes
  
  # Initialize storage for long-format data
  dta <- data.frame(
    Cluster = rep(1:N_total, each = T_total),
    Time = rep(1:T_total, times = N_total),
    X.mean.j = rep(X.mean.j, each = T_total),
    b0 = rep(b0.j, each = T_total),
    b1 = rep(b1.j, each = T_total),
    p.X.mean.j = rep(p.X.mean.j, each = T_total),
    X = NA, eta = NA, Y = NA, p.Y = NA
  )
  
  # Generate Level 1 data
  for (j in 1:N_total) {
    idx <- which(dta$Cluster == j)  # Rows belonging to cluster j
    
    # Generate predictor and linear predictor eta
    if (predictor.type == "continuous") {
      X.jt <- rnorm(T_total, mean = X.mean.j[j], sd = sdX.within)
      eta.jt <- b0.j[j] + b1.j[j] * (X.jt - X.mean.j[j])
    } else {
      X.jt <- rbinom(T_total, 1, p.X.mean.j[j])
      # eta.jt <- b0.j[j] + b1.j[j] * X.jt
      eta.jt <- b0.j[j] + b1.j[j] * (X.jt - X.mean.j[j])
    }
    
    if (outcome.type == "continuous") {
      # add residual error to linear predictor
      Y.jt <- eta.jt + rnorm(T_total, mean = 0, sd = sd.e)
      p.Y.jt <- NA  # No probability needed for continuous outcome
    } else {
      # convert the log-odds (logit) to probability
      p.Y.jt <- plogis(eta.jt)
      # use probability to generate binary outcome
      Y.jt <- rbinom(T_total, 1, p.Y.jt)
    }
    
    # Store values
    dta$X[idx] <- X.jt
    dta$eta[idx] <- eta.jt
    dta$Y[idx] <- Y.jt
    dta$p.Y[idx] <- p.Y.jt
  }
  
  # Compute cluster-level means and centering
  dta$X.cluster.means <- ave(dta$X, dta$Cluster, FUN = mean)
  dta$X.cent <- dta$X - dta$X.cluster.means # using estimated mean
  # dta$X.cent_true <- dta$X - dta$X.mean.j # using true mean
  
  return(dta)
}

### FUNCTION TESTING ###

if(0){
  
  ### Make Evaluation Plots ###
  library(ggplot2)
  
  ### 1 generate continuous X and Y
  # data_cont <- glmm_data_generation(N_total = 5000, T_total = 4, predictor.type = "continuous", outcome.type = "continuous",
  #                                   sdX.within = 1, sdX.between = 2, g.00 = 0, g.01 = 2, sd.u0 = 1,
  #                                   g.10 = 1, sd.u1 = 0, sd.e = 1) # original values
  data_cont <- glmm_data_generation(N_total = 5000, T_total = 10, predictor.type = "continuous", outcome.type = "continuous",
                                    sdX.within = 0.25, sdX.between = 0.5, g.00 = 0, g.01 = 1, sd.u0 = 0.5,
                                    g.10 = 0.5, sd.u1 = 0, sd.e = 0.5)
  # decreased values to ensure that eta values are not too extreme (within -3 to 3 range)
  summary(data_cont)
  
  # 1.1 check distributions of X and Y
  ggplot(data_cont, aes(x = Y)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of Y")
  ggplot(data_cont, aes(x = X)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of X")
  
  # 1.2 check histogram of eta
  ggplot(data_cont, aes(x = eta)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of eta")
  
  ### 2 generate binary X and continuous Y
  data_binx_conty <- glmm_data_generation(N_total = 5000, T_total = 10, predictor.type = "binary", outcome.type = "continuous",
                                    sdX.within = NA, sdX.between = 0.5, g.00 = -0.25, g.01 = 1, sd.u0 = 0.5,
                                    g.10 = 0.5, sd.u1 = 0, sd.e = 0.5)
  # same values as (1) but without sdX.within
  summary(data_binx_conty)
  
  # 2.1 check distributions of X and Y
  ggplot(data_binx_conty, aes(x = Y)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of Y")
  ggplot(data_binx_conty, aes(x = X)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of X")
  
  # 2.2 histograms X.mean.j, p.X.mean.j and eta
  ggplot(data_binx_conty, aes(x = X.mean.j)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of X.mean.j")
  ggplot(data_binx_conty, aes(x = p.X.mean.j)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of p.X.mean.j")
  ggplot(data_binx_conty, aes(x = eta)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of eta")
  # eta is nicely distributed with -3 to 3 range
  
  # when sdX.between is large (e.g, sqrt(4)) the distribution of probabilities 
  # of X and Y are U/bowl shaped (higher density at the edges)
  # when sdX.between is small (e.g., sqrt(0.4)) the distribution of probabilities
  # of X and Y are inverted U shaped (higher density in the center)
  
  # 3 generate continuous X and binary Y
  data_contx_biny <- glmm_data_generation(N_total = 5000, T_total = 10, predictor.type = "continuous", outcome.type = "binary",
                                          sdX.within = 0.25, sdX.between = 0.5, g.00 = 0, g.01 = 1, sd.u0 = 0.5,
                                          g.10 = 0.5, sd.u1 = 0, sd.e = NA)
  summary(data_contx_biny)
  # I lowered the values until the range of eta was approximately -3 to 3, roughly corresponding to probabilities of 0.05 to 0.95.
  # This was done to ensure that the probabilities are not too close to 0 or 1, which would make the logit transformation unstable.
  
  # 3.1 check distributions of X and Y
  ggplot(data_contx_biny, aes(x = Y)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of Y")
  ggplot(data_contx_biny, aes(x = X)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of X")
  
  # 3.2 check distributions of eta and p_y
  ggplot(data_contx_biny, aes(x = eta)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of eta")
  ggplot(data_contx_biny, aes(x = p.Y)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of p_y")
  
  # 4 generate binary X and Y
  data_binary <- glmm_data_generation(N_total = 5000, T_total = 10, predictor.type = "binary", outcome.type = "binary",
                                      sdX.within = NA, sdX.between = 0.5, g.00 = -0.25, g.01 = 0.5, sd.u0 = 0.5,
                                      g.10 = 0.5, sd.u1 = 0, sd.e = NA)
  summary(data_binary)
  # similar to (3), I lowered the values until the range of eta was approximately -3 to 3, roughly corresponding to probabilities of 0.05 to 0.95.
  # In addition, I made sure eta was centered approximately at 0 to ensure balanced probabilities (mean of 0.5)

  # 4.1 check distributions of X and Y
  ggplot(data_binary, aes(x = Y)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of Y")
  ggplot(data_binary, aes(x = X)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of X")
  
  
  # 4.2 histograms with density overlay for evaluation of probability distributions
  ggplot(data_binary, aes(x = eta)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of eta")
  ggplot(data_binary, aes(x = p.X.mean.j)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of p.X.mean.j")
  ggplot(data_binary, aes(x = p.Y)) + geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(alpha = 0.5, fill = "blue") + labs(title = "Histogram of p_y")
  
  # when sdX.between is large (e.g, sqrt(4)) the distribution of probabilities 
  # of X and Y are U/bowl shaped (higher density at the edges)
  # when sdX.between is small (e.g., sqrt(0.4)) the distribution of probabilities 
  # of X and Y are inverted U shaped (higher density in the center)
  # when sdX.between is moderate (e.g., sqrt(2)) the distribution of probabilities
  # is more uniform but inverted U shaped for X and slightly U shaped for Y.
  
  # 4.3 scatterplots for evaluation of relationship between probabilities and originating variables
  
  ggplot(data_binary, aes(x = X.mean.j, y = p.X.mean.j)) + geom_point() + labs(title = "Scatterplot of X.mean.j and p.X.mean.j")
  # As expected, we see a a sigmoid curve representing the logit relationship between X.mean.j 
  # and p.X.mean.j. It is more pronounced (different from linear) when sdX.between is large.
  
  ggplot(data_binary, aes(x = eta, y = p.Y)) + geom_point() + labs(title = "Scatterplot of eta and p_y")
  # As expected again, we see a sigmoid curve representing the logit relationship between eta
  # and p_y. It is more pronounced (different from linear) when sdX.between and sd.u0 are large.
}