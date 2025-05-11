library(lme4)
library(geepack) # geepack yields more stable results than gee package: https://www2.stat.duke.edu/~fl35/teaching/610-23F/docs/slides/6-1-GEE.pdf

# Model Fitting Function
glmm_model_fitting <- function(data, outcome.type) {
  models <- list()
  
  # Define family argument for GLMM and GEE
  family_arg <- if (outcome.type == "continuous") gaussian(link = "identity") else binomial(link = "logit")
  
  # Define formula for different models
  formulas <- list(
    l1 = Y ~ X + (1 | Cluster),
    l2 = Y ~ X.cent + (1 | Cluster),
    l3a = Y ~ X.cent + X.cluster.means + (1 | Cluster),
    l4 = Y ~ X + X.cluster.means + (1 | Cluster)
  )
  
  gee_formulas <- list(
    g1 = Y ~ X,
    g2 = Y ~ X.cent,
    g3a = Y ~ X.cent + X.cluster.means,
    g4 = Y ~ X + X.cluster.means
  )
  
  # Create NA output template (list with correct names)
  coef_NA_fill <- list(
    l1 = setNames(rep(NA, 2), c("(Intercept)", "X")),
    l2 = setNames(rep(NA, 2), c("(Intercept)", "X.cent")),
    l3a = setNames(rep(NA, 3), c("(Intercept)", "X.cent", "X.cluster.means")),
    l4 = setNames(rep(NA, 3), c("(Intercept)", "X", "X.cluster.means"))
  )
  
  if (outcome.type == "continuous") {
    for (name in names(formulas)) {
      models[[name]] <- tryCatch(
        
        # Use MLE (faster and we are not interested in variance components)
        fixef(lmer(formulas[[name]], data = data, REML = FALSE)),
        
        # Catch warnings and errors in lmer() and return NA (e.g., non-convergence, singularity)
        warning = function(w) {
          message("Warning in lmer: ", conditionMessage(w))
          # Return NA for the coefficients
          return(coef_NA_fill[[name]])
        },
        error = function(e) {
          message("Error in lmer: ", conditionMessage(e))
          # Return NA for the coefficients
          return(coef_NA_fill[[name]])
        }
      )
    }
  } else if (outcome.type == "binary") {
    # Fit GLMM models
    for (name in names(formulas)) {
      models[[name]] <- tryCatch(
        fixef(glmer(formulas[[name]], data = data, family = family_arg)),
        
        # Catch warnings and errors in glmer() and return NA (e.g., non-convergence, singularity)
        warning = function(w) {
          message("Warning in glmer: ", conditionMessage(w))
          # Return NA for the coefficients
          return(coef_NA_fill[[name]])
        },
        error = function(e) {
          message("Error in glmer: ", conditionMessage(e))
          # Return NA for the coefficients
          return(coef_NA_fill[[name]])
        }
      )
    }
  } 
  
  # Fit GEE models with different correlation structures using geepack::geeglm
  cor_structures <- c("independence", "exchangeable", "ar1")
  
  for (i in seq_along(gee_formulas)) {
    for (cor in cor_structures) {
      model_name <- paste0("g.", cor, i)
      models[[model_name]] <- tryCatch(
        coef(geeglm(gee_formulas[[i]], id = Cluster, data = data, corstr = cor, family = family_arg, control = geese.control(maxit = 50))),
        
        # Catch errors in geeglm() and return NA (function does not have warnings!)
        error = function(e) {
          message("Error in geeglm: ", conditionMessage(e))
          # Return NA for the coefficients
          return(coef_NA_fill[[i]])
        }
      )
    }
  }
  
  return(models)
}