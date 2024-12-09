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


### Data Preparation

gm_all <- readRDS("simulation_results/GM123ad-1000reps-researchreport/results_beta0_bias.RDS")
# gm123_n200 <- gm123 %>% filter(N == 200)

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

# rename method names
gm_all_long$method <- recode(
  gm_all_long$method,
  "gee_ex_beta0_bias" = "GEE exchangability",
  "gee_ind_beta0_bias" = "GEE independence",
  "gee_ar1_beta0_bias" = "GEE AR(1)",
  "mlm_beta0_bias" = "MLM"
)

# rename T to TT
gm_all_long$TT <- as.factor(gm_all_long$T)
gm_all_long$GM <- as.factor(gm_all_long$GM)
gm_all_long <- gm_all_long %>% filter(GM != 2) # values are way too high

# remove gee AR(1) and GEE exhangeability
# gm_all_long <- gm_all_long %>% filter(method != "GEE AR(1)", method != "GEE exchangability")





### Make grouped barcharts for each GM

# grouped by estimation method
ggplot(gm_all_long, aes(x = method, y = values, fill = method)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_minimal() +
  facet_wrap(~GM)

# grouped by T with different colors for N; separate plots for each GM
ggplot(gm_all_long %>% filter(method == "MLM"), aes(x = TT, y = values, fill = as.factor(N))) +
  geom_col() +
  # geom_bar(position = "dodge", stat = "identity") +
  theme_minimal() +
  facet_grid(.~GM)
  facet_wrap(~GM) 
  

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

# categorize based on estimation method

mlm_data <- gm_all_long %>% filter(method == "MLM")
gee_ind_data <- gm_all_long %>% filter(method == "GEE independence")
gee_ex_data <- gm_all_long %>% filter(method == "GEE exchangability")
gee_ar1_data <- gm_all_long %>% filter(method == "GEE AR(1)")

# create plot for each of the estimation methods with GM as x-axis and values as y-axis
# and GM as color
mlm_plot <- ggplot(mlm_data, aes(x = N, y = values, color = TT)) +
  geom_point() +
  geom_line() +
  labs(title = "MLM", x = "N", y = "Bias") +
  theme_minimal() +
  facet_wrap(~GM)

mlm_plot <- ggplot(gm_all_long %>% filter(GM == 3, method == "MLM"), aes(x = N, y = values, color = TT)) +
  geom_point() +
  geom_line() +
  labs(title = "MLM", x = "N", y = "Bias") +
  theme_minimal() +
  facet_wrap(~GM)

# Create the heatmap
ggplot(gm_all_long %>% filter(method == "MLM"), aes(x = TT, y = as.factor(N), fill = values)) +
  geom_tile(color = "white") + # Add a white border around tiles
  facet_wrap(~ GM) +           # Create a separate heatmap for each GM
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + # Bias color scale
  labs(
    title = "Bias of Treatment Effect by T and N",
    x = "Timepoints (T)",
    y = "Sample Size (N)",
    fill = "Bias"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove gridlines
    strip.text = element_text(size = 12) # Adjust facet label size
  )

# Create the heatmap
ggplot(gm_all_long %>% filter(method == "GEE independence"), aes(x = TT, y = as.factor(N), fill = values)) +
  geom_tile(color = "white") + # Add a white border around tiles
  facet_wrap(~ GM) +           # Create a separate heatmap for each GM
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + # Bias color scale
  labs(
    title = "Difference True Conditional and Estimated Marginal Effect by T and N",
    x = "Timepoints (T)",
    y = "Sample Size (N)",
    fill = "Bias"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove gridlines
    strip.text = element_text(size = 12) # Adjust facet label size
  )

ggplot(gm_all_long %>% filter(method == "GEE independence", GM != 2), aes(x = TT, y = as.factor(N), fill = values)) +
  geom_tile(color = "white") + # Add a white border around tiles
  facet_wrap(~ GM) +           # Create a separate heatmap for each GM
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + # Bias color scale
  labs(
    title = "Difference True Conditional and Estimated Marginal Effect by T and N",
    x = "Timepoints (T)",
    y = "Sample Size (N)",
    fill = "Difference"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove gridlines
    strip.text = element_text(size = 12) # Adjust facet label size
  )

