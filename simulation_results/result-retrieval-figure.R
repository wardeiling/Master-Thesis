library(tidyverse)

gm123 <- readRDS("simulation_results/GM123_trio_1000reps/results_beta0_bias.RDS")
gm123_n200 <- gm123 %>% filter(N == 200)
gm3adh <- readRDS("simulation_results/GM3adh_trio_1000reps/results_beta0_bias.RDS")
# remove first two rows
gm3adh <- gm3adh[-c(1, 2), ]
# combine them row-wise
gm_all <- bind_rows(gm123_n200, gm3adh)

### create a lineplot with the estimation methods having different colors

# "gee_ex_beta0_bias", "gee_ind_beta0_bias", "mlm_beta0_bias" to values
# and the method names to the column "method"

# first reformat the data to long format
gm_all_long <- gm_all %>% pivot_longer(
  cols = c("gee_ex_beta0_bias", "gee_ind_beta0_bias", "gee_ar1_beta0_bias","mlm_beta0_bias"),
  names_to = "method",
  values_to = "values"
)

GM1data <- gm_all_long %>% filter(GM == 1)
GM2data <- gm_all_long %>% filter(GM == 2)
GM3data <- gm_all_long %>% filter(GM == 3)
GM3adata <- gm_all_long %>% filter(GM == "3a")
GM3ddata <- gm_all_long %>% filter(GM == "3d")
GM3hdata <- gm_all_long %>% filter(GM == "3h")

# create plot for each of the GMs with T as x-axis and values as y-axis
# and method as color
GM1plot <- ggplot(GM1data, aes(x = T, y = values, color = method)) +
  geom_line() +
  labs(title = "GM1", x = "T", y = "Bias") +
  theme_minimal()




# only select columns with "mlm" or "gee_ind"
gm13_mlm_gee <- gm13 %>% select(GM, "T", N, starts_with("mlm"), starts_with("gee_ind"))
