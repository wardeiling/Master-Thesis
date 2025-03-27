glmm_data_generation <- function(N_total, T_total, predictor.type, outcome.type, 
                                 sdX.within, sdX.between, g.00, g.01, sd.u0, g.10, sd.u1, sd.e,
                                 true_cluster_means = FALSE){
  
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
  # true_cluster_means whether to use true cluster means for centering
  
  # Check input
  stopifnot(predictor.type %in% c("continuous", "binary"))
  stopifnot(outcome.type %in% c("continuous", "binary"))
  
  # Precompute Level 2 (Cluster-Level) Variables
  X.mean.j <- rnorm(N_total, mean = 0, sd = sdX.between)  # Cluster means
  b1.j <- g.10 + rnorm(N_total, 0, sd.u1)  # Random slopes
  
  if (predictor.type == "continuous"){
    p.X.mean.j <- rep(NA, N_total)  # No probability needed for continuous predictor
    b0.j <- g.00 + g.01 * X.mean.j + rnorm(N_total, 0, sd.u0)  # Random intercepts
  } else if (predictor.type == "binary"){
    p.X.mean.j <- plogis(X.mean.j)  # convert the log-odds (logit) to probability
    b0.j <- g.00 + g.01 * p.X.mean.j + rnorm(N_total, 0, sd.u0)  # Random intercepts
  }
  
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
    if (predictor.type == "continuous"){
      X.jt <- rnorm(T_total, mean = X.mean.j[j], sd = sdX.within)
      eta.jt <- b0.j[j] + b1.j[j] * (X.jt - X.mean.j[j])
    } else if (predictor.type == "binary"){
      X.jt <- rbinom(T_total, 1, p.X.mean.j[j])
      eta.jt <- b0.j[j] + b1.j[j] * (X.jt - p.X.mean.j[j])
    }
    
    if (outcome.type == "continuous"){
      Y.jt <- eta.jt + rnorm(T_total, mean = 0, sd = sd.e) # add residual error to linear predictor
      p.Y.jt <- NA  # No probability needed for continuous outcome
    } else if (outcome.type == "binary"){
      p.Y.jt <- plogis(eta.jt) # convert the log-odds (logit) to probability
      Y.jt <- rbinom(T_total, 1, p.Y.jt) # use probability to generate binary outcome
    }
    
    # Store values
    dta$X[idx] <- X.jt
    dta$eta[idx] <- eta.jt
    dta$Y[idx] <- Y.jt
    dta$p.Y[idx] <- p.Y.jt
  }
  
  if (true_cluster_means == FALSE) {
    # Compute cluster-level means and centering
    dta$X.cluster.means <- ave(dta$X, dta$Cluster, FUN = mean)
    dta$X.cent <- dta$X - dta$X.cluster.means # using estimated mean
  } else { # FOR TESTING PURPOSES ONLY!
    if (predictor.type == "continuous"){
      dta$X.cluster.means <- dta$X.mean.j
    } else if (predictor.type == "binary"){
      dta$X.cluster.means <- dta$p.X.mean.j
    }
    # Use true cluster means for centering
    dta$X.cent <- dta$X - dta$X.cluster.means 
  }
  
  return(dta)
}
