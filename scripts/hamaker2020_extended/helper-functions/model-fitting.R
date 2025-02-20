library(lme4)
library(geepack) # geepack yields more stable results than gee package: https://www2.stat.duke.edu/~fl35/teaching/610-23F/docs/slides/6-1-GEE.pdf

# Improved Model Fitting Function
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
  
  # Fit GLMM models
  for (name in names(formulas)) {
    models[[name]] <- tryCatch(
      summary(glmer(formulas[[name]], data = data, family = family_arg))$coefficients[, 1],
      error = function(e) NA
    )
  }
  
  # Fit GEE models with different correlation structures using geepack::geeglm
  cor_structures <- c("independence", "exchangeable", "ar1")
  
  for (i in seq_along(gee_formulas)) {
    for (cor in cor_structures) {
      model_name <- paste0("g.", cor, i)
      models[[model_name]] <- tryCatch(
        coef(geeglm(gee_formulas[[i]], id = Cluster, data = data, corstr = cor, family = family_arg)),
        error = function(e) NA
      )
    }
  }
  
  return(models)
}
