
library(tidyverse)

gm123 <- readRDS("simulation_results/GM123_trio_1000reps/results_beta0_bias_sd.RDS")
gm13 <- gm123 %>% filter(GM != "2")

# only select columns with "mlm" or "gee_ind"
gm13_mlm_gee <- gm13 %>% select(GM, "T", N, starts_with("mlm"), starts_with("gee_ind"))

library(xtable)

# create a table for the results
tabel <- xtable(gm13_mlm_gee, caption = "Results of the simulation study for GM13", label = "tab:gm13_mlm_gee")
print(tabel, type = "latex", include.rownames = FALSE)
