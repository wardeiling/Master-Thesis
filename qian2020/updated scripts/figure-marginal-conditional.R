# This script generates a figure that shows how marginal and conditional effects may not 
# be equivalent under certain conditions (GM3) in a multilevel linear model.

library(ggplot2)
library(tidyr)
library(dplyr)

# source the generative model
source("qian2020/updated scripts/generative_model_updated_newmodels.R")

# set the seed
set.seed(123)

# simulate data
dta1 <- dgm_with_treatment(sample_size = 10, total_T = 10, dgm_type = "1")
dta3 <- dgm_with_treatment(sample_size = 10, total_T = 10, dgm_type = "3")

# compute the conditional and marginal effects
# conditional: average effect 
# marginal: population average effect

dta1 <- dta1 %>%
  group_by(A) %>%
  mutate(Y_conditional = alpha_0 + alpha_1 * X + b0 + b1 * X + b2 * A + b3 * X * A + eps) %>%
  ungroup() %>%
  group_by(X) %>%
  mutate(Y_marginal = mean(Y_conditional)) %>%
  ungroup()



### plot the data 

# A on X-axis, Y on Y-axis
# multiple grey lines for 10 different individuals with increasing random effects
# black line for conditional average: the subject with b_i = 0, the mean individual
# red line for population average: the average for a given A_it over all individuals


# plot for dgm_type = 1

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Source the generative model
source("qian2020/updated scripts/generative_model_updated_newmodels.R")

# Set the seed for reproducibility
set.seed(123)

# Simulate data for GM1 and GM3
dta1 <- dgm_with_treatment(sample_size = 10, total_T = 10, dgm_type = "1")
dta3 <- dgm_with_treatment(sample_size = 10, total_T = 10, dgm_type = "3")

# Function to compute marginal and conditional averages
compute_effects <- function(data) {
  data %>%
    group_by(A) %>%
    summarize(
      marginal_avg = mean(Y),  # Population average for given A
      conditional_avg = mean(Y[b0 == 0])  # Mean for individuals with b0 = 0
    )
}

# Compute effects for GM1 and GM3
effects1 <- compute_effects(dta1)
effects3 <- compute_effects(dta3)

# Plot function
plot_effects <- function(data, effects, title) {
  ggplot(data, aes(x = A, y = Y, group = userid)) +
    geom_line(aes(color = factor(userid)), size = 0.7, alpha = 0.6, show.legend = FALSE) +
    geom_line(data = effects, aes(x = A, y = conditional_avg), color = "black", size = 1.2) +
    geom_line(data = effects, aes(x = A, y = marginal_avg), color = "red", size = 1.2, linetype = "dashed") +
    labs(
      title = title,
      x = "Treatment (A)",
      y = "Outcome (Y)"
    ) +
    theme_minimal() +
    scale_color_grey()
}

# Create plots
plot1 <- plot_effects(dta1, effects1, "GM1: Marginal and Conditional Effects")
plot3 <- plot_effects(dta3, effects3, "GM3: Marginal and Conditional Effects")

# Display plots side by side
library(gridExtra)
grid.arrange(plot1, plot3, ncol = 2)
