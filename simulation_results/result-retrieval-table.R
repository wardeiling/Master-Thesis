
library(tidyverse)

results <- readRDS("simulation_results/GM123ad-1000reps-researchreport/results_beta0_bias_sd_success.RDS")
# remove columns which contain "ar1" or "ex" in the column name
results2 <- results %>% select(-contains("ar1"), -contains("ex"), -"gee_ind_success")

# create a table for the results
library(xtable)

# make table for beta0 bias with Standarddeviation and success rate
colnames(results2) <- c("GM", "T", "N", "MLM_bias", "MLM_sd", "GEE-Ind_bias", "GEE-Ind_sd", "MLM_success")
print(xtable(results2, digits = c(0, 0, 0, 0, rep(3, 5)), 
             caption = "Results for beta0 bias with Standarddeviation and success rate, 1000 replications", label = "tab:beta0_bias_sd_success"), 
      include.rownames = FALSE, file = "simulation_results/GM123ad-1000reps-researchreport/results_cleaned.tex")
