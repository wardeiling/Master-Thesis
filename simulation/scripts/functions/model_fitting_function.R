library(lme4)
library(geepack)

fit_models <- function(data, gm_type) {
  # Initialize result storage
  estimates <- list(
    mlm = NA,
    gee_ind = NA,
    gee_exch = NA,
    gee_ar1 = NA
  )
  
  # Fit MLM
  estimates$mlm <- tryCatch({
    if (gm_type %in% c(1, 3)) {
      coef(summary(lmer(Y ~ Z * X + (1 + X | userid), data = data, REML = TRUE)))["X", "Estimate"]
    } else if (gm_type == 2) {
      coef(summary(lmer(Y ~ Z * X + (Z * X | userid), data = data, REML = TRUE)))["X", "Estimate"]
    }
  }, error = function(e) {
    NA  # Assign NA if an error occurs
  })
  
  # Fit GEE models with different correlation structures
  estimates$gee_ind <- tryCatch({
    coef(summary(geeglm(Y ~ Z * X, id = userid, data = data, corstr = "independence")))["X","Estimate"]
  }, error = function(e) {
    NA
  })
  
  estimates$gee_exch <- tryCatch({
    coef(summary(geeglm(Y ~ Z * X, id = userid, data = data, corstr = "exchangeable")))["X","Estimate"]
  }, error = function(e) {
    NA
  })
  
  estimates$gee_ar1 <- tryCatch({
    coef(summary(geeglm(Y ~ Z * X, id = userid, data = data, corstr = "ar1")))["X","Estimate"]
  }, error = function(e) {
    NA
  })
  
  return(estimates)
}

# fit_models <- function(data, gm_type) {
#   # Fit MLM
#   if (gm_type %in% c(1, 3)) {
#     mlm <- lmer(Y ~ Z * X + (1 + X | userid), data = data, REML = TRUE)
#   } else if (gm_type == 2) {
#     mlm <- lmer(Y ~ Z * X + (Z * X | userid), data = data, REML = TRUE)
#   }
#   
#   # Fit GEE with different correlation structures
#   gee_ind <- geeglm(Y ~ Z * X, id = userid, data = data, corstr = "independence")
#   gee_exch <- geeglm(Y ~ Z * X, id = userid, data = data, corstr = "exchangeable")
#   gee_ar1 <- geeglm(Y ~ Z * X, id = userid, data = data, corstr = "ar1")
#   
#   # Collect estimates for comparison
#   estimates <- list(
#     mlm = coef(summary(mlm)),
#     gee_ind = coef(summary(gee_ind)),
#     gee_exch = coef(summary(gee_exch)),
#     gee_ar1 = coef(summary(gee_ar1))
#   )
#   
#   return(estimates)
# }
