library(tidyverse)

### Data Preparation

gm123 <- readRDS("simulation_results/GM123_trio_1000reps/results_beta0_bias.RDS")
gm123_n200 <- gm123 %>% filter(N == 200)
gm3adh <- readRDS("simulation_results/GM3adh_trio_1000reps/results_beta0_bias.RDS")
# remove first two rows
gm3adh <- gm3adh[-c(1, 2), ]
# combine them row-wise
gm_all <- bind_rows(gm123_n200, gm3adh)

# first reformat the data to long format
# "gee_ex_beta0_bias", "gee_ind_beta0_bias", "mlm_beta0_bias" to values
# and the method names to the column "method"
gm_all_long <- gm_all %>% pivot_longer(
  cols = c(
    "gee_ex_beta0_bias",
    "gee_ind_beta0_bias",
    "gee_ar1_beta0_bias",
    "mlm_beta0_bias"
  ),
  names_to = "method",
  values_to = "values"
)
  
# remove GM2
gm_all_long <- gm_all_long %>% filter(GM != 2)
# rename method names
gm_all_long$method <- recode(
  gm_all_long$method,
  "gee_ex_beta0_bias" = "GEE exchangability",
  "gee_ind_beta0_bias" = "GEE independence",
  "gee_ar1_beta0_bias" = "GEE AR(1)",
  "mlm_beta0_bias" = "MLM"
)

# remove gee AR(1) and GEE exhangeability
gm_all_long <- gm_all_long %>% filter(method != "GEE AR(1)", method != "GEE exchangability")



### Make grouped barcharts for each GM

# grouped by estimation method
ggplot(gm_all_long, aes(x = method, y = values, fill = GM)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_minimal()

# grouped by GM
ggplot(gm_all_long, aes(x = GM, y = values, fill = method)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_minimal()

### create a lineplot with the estimation methods having different colors

GM1data <- gm_all_long %>% filter(GM == 1)
GM2data <- gm_all_long %>% filter(GM == 2)
GM3data <- gm_all_long %>% filter(GM == 3)
GM3adata <- gm_all_long %>% filter(GM == "3a")
GM3ddata <- gm_all_long %>% filter(GM == "3d")
GM3hdata <- gm_all_long %>% filter(GM == "3h")

# create plot for each of the GMs with T as x-axis and values as y-axis
# and method as color
GM1plot <- ggplot(GM1data, aes(x = T, y = values, color = method)) +
  geom_point() +
  geom_line() +
  labs(title = "GM1", x = "T", y = "Bias") +
  theme_minimal()


