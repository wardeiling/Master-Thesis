library(lme4)
library(geepack) # geepack yields more stable results than gee package: https://www2.stat.duke.edu/~fl35/teaching/610-23F/docs/slides/6-1-GEE.pdf

# Model Fitting Function
glmm_model_fitting <- function(data, outcome.type) {
  
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
  
  ### GENERALIZED LINEAR MIXED MODELS ###
  
  # Fit multilevel linear models
  if (outcome.type == "continuous"){
    
    solution_l1 <- tryCatch(
      fixef(lmer(formulas$l1, data = data, REML = FALSE)),
      error = function(e) NA
    )
    
    solution_l2 <- tryCatch(
      fixef(lmer(formulas$l2, data = data, REML = FALSE)),
      error = function(e) NA
    )
    
    solution_l3a <- tryCatch(
      fixef(lmer(formulas$l3a, data = data, REML = FALSE)),
      error = function(e) NA
    )
    
    solution_l4 <- tryCatch(
      fixef(lmer(formulas$l4, data = data, REML = FALSE)),
      error = function(e) NA
    )
  
  } else if (outcome.type == "binary") {
    
    # Fit GLMM with logit link function
    solution_l1 <- tryCatch(
      fixef(glmer(formulas$l1, data = data, family = family_arg)),
      error = function(e) NA
    )
    
    solution_l2 <- tryCatch(
      fixef(glmer(formulas$l2, data = data, family = family_arg)),
      error = function(e) NA
    )
    
    solution_l3a <- tryCatch(
      fixef(glmer(formulas$l3a, data = data, family = family_arg)),
      error = function(e) NA
    )
    
    solution_l4 <- tryCatch(
      fixef(glmer(formulas$l4, data = data, family = family_arg)),
      error = function(e) NA
    )
  }
  
  ### GENERALIZED ESTIMATING EQUATIONS ###
  
  # Fit GEE models with different correlation structures using geepack::geeglm
  cor_structures <- c("independence", "exchangeable", "ar1")
  
  # models with raw predictor X
  solution_g.independence1 <- tryCatch(
    fixef(geeglm(gee_formulas$g1, data = data, id = Cluster, family = family_arg, corstr = cor_structures[1])),
    error = function(e) NA
  )
  
  solution_g.exchangeable1 <- tryCatch(
    fixef(geeglm(gee_formulas$g1, data = data, id = Cluster, family = family_arg, corstr = cor_structures[2])),
    error = function(e) NA
  )
  
  solution_g.ar11 <- tryCatch(
    fixef(geeglm(gee_formulas$g1, data = data, id = Cluster, family = family_arg, corstr = cor_structures[3])),
    error = function(e) NA
  )
  
  # models with centered predictor X
  solution_g.independence2 <- tryCatch(
    fixef(geeglm(gee_formulas$g2, data = data, id = Cluster, family = family_arg, corstr = cor_structures[1])),
    error = function(e) NA
  )
  
  solution_g.exchangeable2 <- tryCatch(
    fixef(geeglm(gee_formulas$g2, data = data, id = Cluster, family = family_arg, corstr = cor_structures[2])),
    error = function(e) NA
  )
  
  solution_g.ar12 <- tryCatch(
    fixef(geeglm(gee_formulas$g2, data = data, id = Cluster, family = family_arg, corstr = cor_structures[3])),
    error = function(e) NA
  )
  
  # models with centered predictor X and cluster means
  solution_g.independence3a <- tryCatch(
    fixef(geeglm(gee_formulas$g3a, data = data, id = Cluster, family = family_arg, corstr = cor_structures[1])),
    error = function(e) NA
  )
  
  solution_g.exchangeable3a <- tryCatch(
    fixef(geeglm(gee_formulas$g3a, data = data, id = Cluster, family = family_arg, corstr = cor_structures[2])),
    error = function(e) NA
  )
  
  solution_g.ar13a <- tryCatch(
    fixef(geeglm(gee_formulas$g3a, data = data, id = Cluster, family = family_arg, corstr = cor_structures[3])),
    error = function(e) NA
  )
  
  # models with raw predictor X and cluster means
  solution_g.independence4 <- tryCatch(
    fixef(geeglm(gee_formulas$g4, data = data, id = Cluster, family = family_arg, corstr = cor_structures[1])),
    error = function(e) NA
  )
  
  solution_g.exchangeable4 <- tryCatch(
    fixef(geeglm(gee_formulas$g4, data = data, id = Cluster, family = family_arg, corstr = cor_structures[2])),
    error = function(e) NA
  )
  
  solution_g.ar14 <- tryCatch(
    fixef(geeglm(gee_formulas$g4, data = data, id = Cluster, family = family_arg, corstr = cor_structures[3])),
    error = function(e) NA
  )
  
  # store all models in a list
  models <- list(
    l1 = solution_l1,
    l2 = solution_l2,
    l3a = solution_l3a,
    l4 = solution_l4,
    g.independence1 = solution_g.independence1,
    g.exchangeable1 = solution_g.exchangeable1,
    g.ar11 = solution_g.ar11,
    g.independence2 = solution_g.independence2,
    g.exchangeable2 = solution_g.exchangeable2,
    g.ar12 = solution_g.ar12,
    g.independence3a = solution_g.independence3a,
    g.exchangeable3a = solution_g.exchangeable3a,
    g.ar13a = solution_g.ar13a,
    g.independence4 = solution_g.independence4,
    g.exchangeable4 = solution_g.exchangeable4,
    g.ar14 = solution_g.ar14
  )
  
  return(models)
}
