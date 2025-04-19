# load packages
library(purrr) # for functional programming
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation
library(stringr) # for string manipulation
library(ggplot2) # for plotting

### RETRIEVE RESULTS ### ----

### FULL SIMULATION (PART 1) ----
runname <- "April10_fullsimulation"

# retrieve the design and parameterization
settings <- readRDS(paste0("simulation_results_glmm/", runname, "/settings.RDS"))
design <- settings$design
parametrization <- settings$parametrization
nsim <- settings$nsim

# repeat each row of the design nsim times and create replication index
design_long <- design[rep(1:nrow(design), each = nsim), ] %>%
  mutate(replication = rep(1:nsim, times = nrow(design)),
         design_id = rep(1:nrow(design), each = nsim))

# initialize list to store simulation results
results_list <- list()

# loop over the design
for (idesign in 1:nrow(design)) {
  
  ### Extract parameter values from design
  g.00 <- design$g.00[idesign]
  g.01 <- design$g.01[idesign]
  g.10 <- design$g.10[idesign]
  
  # determine the beta values based on the parametrization
  if (parametrization == "centeredX") {
    beta_between <- g.01
    beta_within <- g.10
    beta_contextual <- beta_between - beta_within
  } else if (parametrization == "mundlak") {
    beta_contextual <- g.01
    beta_within <- g.10
    beta_between <- beta_within + beta_contextual
  }
  
  # read in the results
  parallel_results_setting <- readRDS(paste0("simulation_results_glmm/", runname, "/", idesign, ".RDS"))
  
  # unlist the lists inside the list
  df <- map_dfr(parallel_results_setting, function(rep) {
    map_dfr(rep, ~ as.data.frame(as.list(.x)), .id = "model")
  }, .id = "replication")
  
  df_wide <- df %>%
    select(-X.Intercept.) %>%
    pivot_wider(id_cols = replication, names_from = model, values_from = c("X", "X.cent", "X.cluster.means"), names_glue = "{model}_{.value}") %>%
    select(where(~ !all(is.na(.)))) %>%
    mutate(replication = as.integer(replication))  # ensure numeric
  
  # compute bias
  df_wide2 <- df_wide %>%
    mutate(l1_g.10_bias = l1_X - beta_within,
           l2_g.10_bias = l2_X.cent - beta_within,
           l3a_g.10_bias = l3a_X.cent - beta_within,
           l3a_g.01_bias = l3a_X.cluster.means - beta_between,
           l4_g.10_bias = l4_X - beta_within,
           l4_g.01_bias = l4_X.cluster.means - beta_contextual,
           # also for GEE
           # parametrization 1
           g.independence1_g.10_bias = g.independence1_X - beta_within,
           g.exchangeable1_g.10_bias = g.exchangeable1_X - beta_within,
           g.ar11_g.10_bias = g.ar11_X - beta_within,
           # parametrization 2
           g.independence2_g.10_bias = g.independence2_X.cent - beta_within,
           g.exchangeable2_g.10_bias = g.exchangeable2_X.cent - beta_within,
           g.ar12_g.10_bias = g.ar12_X.cent - beta_within,
           # parametrization 3
           g.independence3_g.10_bias = g.independence3_X.cent - beta_within,
           g.independence3_g.01_bias = g.independence3_X.cluster.means - beta_between,
           g.exchangeable3_g.10_bias = g.exchangeable3_X.cent - beta_within,
           g.exchangeable3_g.01_bias = g.exchangeable3_X.cluster.means - beta_between,
           g.ar13_g.10_bias = g.ar13_X.cent - beta_within,
           g.ar13_g.01_bias = g.ar13_X.cluster.means - beta_between,
           # parametrization 4
           g.independence4_g.10_bias = g.independence4_X - beta_within,
           g.independence4_g.01_bias = g.independence4_X.cluster.means - beta_contextual,
           g.exchangeable4_g.10_bias = g.exchangeable4_X - beta_within,
           g.exchangeable4_g.01_bias = g.exchangeable4_X.cluster.means - beta_contextual,
           g.ar14_g.10_bias = g.ar14_X - beta_within,
           g.ar14_g.01_bias = g.ar14_X.cluster.means - beta_contextual,
           # add design_id for merging later
           design_id = idesign) 
  
  # store the result
  results_list[[idesign]] <- df_wide2
}

# combine all results
sim_results_all <- bind_rows(results_list)

# merge with the long design
final_df <- left_join(design_long, sim_results_all, by = c("design_id", "replication"))

# save the final data frame
saveRDS(final_df, paste0("simulation_results_glmm/", runname, "/plotting_bias_df1.RDS"))

### REMAINDER OF SIMULATION (PART 2) ----

runname <- "April17_fullsimulation_contXY"

# retrieve the design and parameterization
settings <- readRDS(paste0("simulation_results_glmm/", runname, "/settings.RDS"))
design <- settings$design
parametrization <- settings$parametrization
nsim <- settings$nsim

# repeat each row of the design nsim times and create replication index
design_long <- design[rep(1:nrow(design), each = nsim), ] %>%
  mutate(replication = rep(1:nsim, times = nrow(design)),
         design_id = rep(1:nrow(design), each = nsim) + 378)

# initialize list to store simulation results
results_list <- list()

# loop over the design
for (idesign in 1:nrow(design)) {
  
  ### Extract parameter values from design
  g.00 <- design$g.00[idesign]
  g.01 <- design$g.01[idesign]
  g.10 <- design$g.10[idesign]
  
  # determine the beta values based on the parametrization
  if (parametrization == "centeredX") {
    beta_between <- g.01
    beta_within <- g.10
    beta_contextual <- beta_between - beta_within
  } else if (parametrization == "mundlak") {
    beta_contextual <- g.01
    beta_within <- g.10
    beta_between <- beta_within + beta_contextual
  }
  
  # read in the results
  parallel_results_setting <- readRDS(paste0("simulation_results_glmm/", runname, "/", idesign, ".RDS"))
  
  # unlist the lists inside the list
  df <- map_dfr(parallel_results_setting, function(rep) {
    map_dfr(rep, ~ as.data.frame(as.list(.x)), .id = "model")
  }, .id = "replication")
  
  df_wide <- df %>%
    select(-X.Intercept.) %>%
    pivot_wider(id_cols = replication, names_from = model, values_from = c("X", "X.cent", "X.cluster.means"), names_glue = "{model}_{.value}") %>%
    select(where(~ !all(is.na(.)))) %>%
    mutate(replication = as.integer(replication))  # ensure numeric
  
  # compute bias
  df_wide2 <- df_wide %>%
    mutate(l1_g.10_bias = l1_X - beta_within,
           l2_g.10_bias = l2_X.cent - beta_within,
           l3a_g.10_bias = l3a_X.cent - beta_within,
           l3a_g.01_bias = l3a_X.cluster.means - beta_between,
           l4_g.10_bias = l4_X - beta_within,
           l4_g.01_bias = l4_X.cluster.means - beta_contextual,
           # also for GEE
           # parametrization 1
           g.independence1_g.10_bias = g.independence1_X - beta_within,
           g.exchangeable1_g.10_bias = g.exchangeable1_X - beta_within,
           g.ar11_g.10_bias = g.ar11_X - beta_within,
           # parametrization 2
           g.independence2_g.10_bias = g.independence2_X.cent - beta_within,
           g.exchangeable2_g.10_bias = g.exchangeable2_X.cent - beta_within,
           g.ar12_g.10_bias = g.ar12_X.cent - beta_within,
           # parametrization 3
           g.independence3_g.10_bias = g.independence3_X.cent - beta_within,
           g.independence3_g.01_bias = g.independence3_X.cluster.means - beta_between,
           g.exchangeable3_g.10_bias = g.exchangeable3_X.cent - beta_within,
           g.exchangeable3_g.01_bias = g.exchangeable3_X.cluster.means - beta_between,
           g.ar13_g.10_bias = g.ar13_X.cent - beta_within,
           g.ar13_g.01_bias = g.ar13_X.cluster.means - beta_between,
           # parametrization 4
           g.independence4_g.10_bias = g.independence4_X - beta_within,
           g.independence4_g.01_bias = g.independence4_X.cluster.means - beta_contextual,
           g.exchangeable4_g.10_bias = g.exchangeable4_X - beta_within,
           g.exchangeable4_g.01_bias = g.exchangeable4_X.cluster.means - beta_contextual,
           g.ar14_g.10_bias = g.ar14_X - beta_within,
           g.ar14_g.01_bias = g.ar14_X.cluster.means - beta_contextual,
           # add design_id for merging later
           design_id = idesign) 
  
  # store the result
  results_list[[idesign]] <- df_wide2
}

# combine all results
sim_results_all <- bind_rows(results_list) %>%
  mutate(design_id = design_id + 378) # adjust design_id to match the first simulation

# merge with the long design
final_df <- left_join(design_long, sim_results_all, by = c("design_id", "replication"))

# save the final data frame
saveRDS(final_df, paste0("simulation_results_glmm/", runname, "/plotting_bias_df2.RDS"))

### COMBINE RESULTS FROM BOTH SIMULATIONS ### ----

# read in the first simulation results
runname1 <- "April10_fullsimulation"
final_df1 <- readRDS(paste0("simulation_results_glmm/", runname1, "/plotting_bias_df1.RDS"))
# read in the second simulation results
runname2 <- "April17_fullsimulation_contXY"
final_df2 <- readRDS(paste0("simulation_results_glmm/", runname2, "/plotting_bias_df2.RDS"))

# combine the two data frames
final_df <- bind_rows(final_df1, final_df2)
# save the final data frame
newrunname <- "April18_fullsimulation_combined"
saveRDS(final_df, paste0("simulation_results_glmm/", newrunname, "/plotting_bias_df.RDS"))

### SELECT RELEVANT CASES FOR PLOTTING ### ----

runname <- "April18_fullsimulation_combined"

# create a matrix with the different combinations of predictor and outcome type
settings <- expand.grid(
  predictor.type = c("binary", "continuous"),
  outcome.type = c("binary", "continuous"),
  stringsAsFactors = FALSE
)
# # remove continuous X and Y setting
# settings <- settings %>%
#   filter(!(predictor.type == "continuous" & outcome.type == "continuous"))

# loop to make all plots
for(i in 1:nrow(settings)) {
  set.predictor.type <- settings$predictor.type[i]
  set.outcome.type <- settings$outcome.type[i]
  
  # create string for the file name
  type <- paste0("pred_", set.predictor.type, "_out_", set.outcome.type, "_")
  
  ### PLOT 1: Grid of sdX.between and g.01 ----
  
  # read in the final data frame
  final_df <- readRDS(paste0("simulation_results_glmm/", runname, "/plotting_bias_df.RDS"))
  
  # select relevant variables and cases (select and filter)
  plot1_df <- final_df %>%
    # select T = 20, N = 200, sd.u0 = 1, predictor.type = "binary" and outcome.type = "continuous"
    filter(T_total == 20, N_total == 200, sd.u0 == 1, predictor.type == set.predictor.type, outcome.type == set.outcome.type) %>%
    select(-c(ends_with("_success"), ends_with("_X"), ends_with("_X.cent"), ends_with("_X.cluster.means"))) 
  
  plot1_df_beta1 <- plot1_df %>%
    select(-ends_with("_g.01_bias")) %>%
    # turn the bias variables into long format, with a new column indicating the model name 
    pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "beta1_bias") %>%
    # remove the "_g.10_bias" suffix from the model names
    mutate(model = str_remove(model, "_g.10_bias")) %>%
    # set factor levels of model to ensure correct order in the plot
    mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2", 
                                            "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
                                            "g.ar14", "g.independence4"))) %>%
    # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
    mutate(sdX.between_str = paste0("sdX.between = ", sdX.between),
           g.01_str = paste0("g.01 = ", g.01)) %>%
    # Set factor levels to ensure correct order in the plot (0, 1, 3)
    mutate(sdX.between_str = factor(sdX.between_str, levels = c("sdX.between = 0", "sdX.between = 1", "sdX.between = 3")),
           g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")))
  
  plot1_df_g01 <- plot1_df %>%
    select(-ends_with("_g.10_bias")) %>%
    # turn the bias variables into long format, with a new column indicating the model name
    pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "g01_bias") %>%
    # remove the "_g.01_bias" suffix from the model names
    mutate(model = str_remove(model, "_g.01_bias")) %>%
    # set factor levels of model to ensure correct order in the plot
    mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2", 
                                            "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
                                            "g.ar14", "g.independence4"))) %>%
    # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
    mutate(sdX.between_str = paste0("sdX.between = ", sdX.between),
           g.01_str = paste0("g.01 = ", g.01)) %>%
    # Set factor levels to ensure correct order in the plot (0, 1, 3)
    mutate(sdX.between_str = factor(sdX.between_str, levels = c("sdX.between = 0", "sdX.between = 1", "sdX.between = 3")),
           g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")))
  
  # For the within-person effect
  ggplot(plot1_df_beta1, aes(x = model, y = beta1_bias, col = model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
    ylim(-1.5, 1.5) +  # Set y-axis limits
    labs(x = "Generative Model", y = "Bias") +
    facet_grid(sdX.between_str ~ g.01_str) + # Show T and N values in labels
    theme_bw() + 
    # remove X axis labels
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  # save
  ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-sd.Xi_within.pdf"), width = 10, height = 8)
  
  # For the between-person effect
  ggplot(plot1_df_g01, aes(x = model, y = g01_bias, col = model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
    ylim(-1.5, 1.5) +  # Set y-axis limits
    labs(x = "Generative Model", y = "Bias") +
    facet_grid(sdX.between_str ~ g.01_str) + # Show T and N values in labels
    theme_bw() + 
    # remove X axis labels
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  # save
  ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-sd.Xi_between.pdf"), width = 10, height = 8)
  
  ### PLOT 2A: Grid of sd.u0 and g.01 ----
  
  # select relevant variables and cases (select and filter)
  plot2_df <- final_df %>%
    # select T = 20, N = 200, sd.u0 = 1, predictor.type = "binary" and outcome.type = "continuous"
    filter(T_total == 20, N_total == 200, sdX.between == 1, predictor.type == set.predictor.type, outcome.type == set.outcome.type) %>%
    select(-c(ends_with("_success"), ends_with("_X"), ends_with("_X.cent"), ends_with("_X.cluster.means"))) 
  
  plot2_df_beta1 <- plot2_df %>%
    select(-ends_with("_g.01_bias")) %>%
    # turn the bias variables into long format, with a new column indicating the model name 
    pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "beta1_bias") %>%
    # remove the "_g.10_bias" suffix from the model names
    mutate(model = str_remove(model, "_g.10_bias")) %>%
    # set factor levels of model to ensure correct order in the plot
    mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2", 
                                            "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
                                            "g.ar14", "g.independence4"))) %>%
    # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
    mutate(g.01_str = paste0("g.01 = ", g.01),
           sd.u0_str = paste0("sd.u0 = ", sd.u0)) %>%
    # Set factor levels to ensure correct order in the plot (0, 1, 3)
    mutate(g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")),
           sd.u0_str = factor(sd.u0_str, levels = c("sd.u0 = 0", "sd.u0 = 1", "sd.u0 = 3")))
  
  plot2_df_g01 <- plot2_df %>%
    select(-ends_with("_g.10_bias")) %>%
    # turn the bias variables into long format, with a new column indicating the model name
    pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "g01_bias") %>%
    # remove the "_g.01_bias" suffix from the model names
    mutate(model = str_remove(model, "_g.01_bias")) %>%
    # set factor levels of model to ensure correct order in the plot
    mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2", 
                                            "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
                                            "g.ar14", "g.independence4"))) %>%
    # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
    mutate(g.01_str = paste0("g.01 = ", g.01),
           sd.u0_str = paste0("sd.u0 = ", sd.u0)) %>%
    # Set factor levels to ensure correct order in the plot (0, 1, 3)
    mutate(g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")),
           sd.u0_str = factor(sd.u0_str, levels = c("sd.u0 = 0", "sd.u0 = 1", "sd.u0 = 3")))
  
  # For the within-person effect
  ggplot(plot2_df_beta1, aes(x = model, y = beta1_bias, col = model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
    ylim(-1.5, 1.5) +  # Set y-axis limits
    labs(x = "Generative Model", y = "Bias") +
    facet_grid(sd.u0_str ~ g.01_str) + # Show T and N values in labels
    theme_bw() + 
    # remove X axis labels
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  # save
  ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-sd.u0_within.pdf"), width = 10, height = 8)
  
  # For the between-person effect
  ggplot(plot2_df_g01, aes(x = model, y = g01_bias, col = model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
    ylim(-1.5, 1.5) +  # Set y-axis limits
    labs(x = "Generative Model", y = "Bias") +
    facet_grid(sd.u0_str ~ g.01_str) + # Show T and N values in labels
    theme_bw() + 
    # remove X axis labels
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  # save
  ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-sd.u0_between.pdf"), width = 10, height = 8)
  
  ### PLOT 3: Grid of T_total and g.01 ----
  
  # select relevant variables and cases (select and filter)
  plot3_df <- final_df %>%
    # select T = 20, N = 200, sd.u0 = 1, predictor.type = "binary" and outcome.type = "continuous"
    filter(N_total == 200, sdX.between == 3, sd.u0 == 1, predictor.type == set.predictor.type, outcome.type == set.outcome.type) %>%
    select(-c(ends_with("_success"), ends_with("_X"), ends_with("_X.cent"), ends_with("_X.cluster.means"))) 
  
  plot3_df_beta1 <- plot3_df %>%
    select(-ends_with("_g.01_bias")) %>%
    # turn the bias variables into long format, with a new column indicating the model name 
    pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "beta1_bias") %>%
    # remove the "_g.10_bias" suffix from the model names
    mutate(model = str_remove(model, "_g.10_bias")) %>%
    # set factor levels of model to ensure correct order in the plot
    mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2", 
                                            "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
                                            "g.ar14", "g.independence4"))) %>%
    # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
    mutate(g.01_str = paste0("g.01 = ", g.01),
           T_total_str = paste0("T = ", T_total)) %>%
    # Set factor levels to ensure correct order in the plot (0, 1, 3)
    mutate(g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")),
           T_total_str = factor(T_total_str, levels = c("T = 5", "T = 10", "T = 20")))
  
  plot3_df_g01 <- plot3_df %>%
    select(-ends_with("_g.10_bias")) %>%
    # turn the bias variables into long format, with a new column indicating the model name
    pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "g01_bias") %>%
    # remove the "_g.01_bias" suffix from the model names
    mutate(model = str_remove(model, "_g.01_bias")) %>%
    # set factor levels of model to ensure correct order in the plot
    mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2", 
                                            "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
                                            "g.ar14", "g.independence4"))) %>%
    # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
    mutate(g.01_str = paste0("g.01 = ", g.01),
           T_total_str = paste0("T = ", T_total)) %>%
    # Set factor levels to ensure correct order in the plot (0, 1, 3)
    mutate(g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")),
           T_total_str = factor(T_total_str, levels = c("T = 5", "T = 10", "T = 20")))
  
  # For the within-person effect
  ggplot(plot3_df_beta1, aes(x = model, y = beta1_bias, col = model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
    ylim(-1.5, 1.5) +  # Set y-axis limits
    labs(x = "Generative Model", y = "Bias") +
    facet_grid(T_total_str ~ g.01_str) + # Show T and N values in labels
    theme_bw() + 
    # remove X axis labels
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  # save
  ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-T_total_within.pdf"), width = 10, height = 8)
  
  # For the between-person effect
  ggplot(plot3_df_g01, aes(x = model, y = g01_bias, col = model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
    ylim(-1.5, 1.5) +  # Set y-axis limits
    labs(x = "Generative Model", y = "Bias") +
    facet_grid(T_total_str ~ g.01_str) + # Show T and N values in labels
    theme_bw() + 
    # remove X axis labels
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  # save
  ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-T_total_between.pdf"), width = 10, height = 8)
  
  
  # ### PLOT 4 ----
  # predictor.type <- "binary"
  # outcome.type <- "continuous"
  # 
  # plot4_df < - final_df %>%
  #   filter(N_total == 200, sdX.between == 3, predictor.type == predictor.type, outcome.type == outcome.type) %>%
  #   select(T_total, predictor.type, outcome.type, ends_with("_g.10_bias")) %>%
  #   pivot_longer(cols = ends_with("_g.10_bias"),
  #                names_to = "model", values_to = "beta1_bias") %>%
  #   mutate(
  #     model = str_remove(model, "_g.10_bias"),
  #     model = factor(model, levels = c("l2", "l3a", "l4", 
  #                                      "g.exchangeable2", "g.ar12", "g.independence2", 
  #                                      "g.exchangeable3", "g.ar13", "g.independence3", 
  #                                      "g.exchangeable4", "g.ar14", "g.independence4")),
  #     predictor.type = factor(predictor.type),
  #     outcome.type = factor(outcome.type)
  #   ) %>%
  #   # remove any bias values exceeding -100 or 100
  #   filter(beta1_bias > -100 & beta1_bias < 100)
  # 
  # # create a function to summarize the data
  # data_summary <- function(data, varname, groupnames){
  #   require(plyr)
  #   summary_func <- function(x, col){
  #     c(mean = mean(x[[col]], na.rm=TRUE),
  #       sd = sd(x[[col]], na.rm=TRUE))
  #   }
  #   data_sum<-ddply(data, groupnames, .fun=summary_func,
  #                   varname)
  #   data_sum <- rename(data_sum, c("mean" = varname))
  #   return(data_sum)
  # }
  # 
  # # employ function to summarize the data
  # plot4_df_summary <- data_summary(plot3_df, varname="beta1_bias", groupnames=c("T_total", "predictor.type", "outcome.type", "model")) %>%
  #   # turn bias values into absolute values
  #   mutate(beta1_bias = abs(beta1_bias))
  # 
  
}

### PLOT 3 ----

runname <- "April18_fullsimulation_combined"

# read in the final data frame
final_df <- readRDS(paste0("simulation_results_glmm/", runname, "/plotting_bias_df.RDS"))

# Prepare data for plotting beta1 bias across timepoints
plot_df_time_beta1 <- final_df %>%
  filter(N_total == 200, sdX.between == 3, sd.u0 == 1, g.01 == 1) %>%
  select(T_total, predictor.type, outcome.type, ends_with("_g.10_bias")) %>%
  pivot_longer(cols = ends_with("_g.10_bias"),
               names_to = "model", values_to = "beta1_bias") %>%
  mutate(
    model = str_remove(model, "_g.10_bias"),
    model = factor(model, levels = c("l2", "l3a", "l4", 
                                     "g.exchangeable2", "g.ar12", "g.independence2", 
                                     "g.exchangeable3", "g.ar13", "g.independence3", 
                                     "g.exchangeable4", "g.ar14", "g.independence4")),
    predictor.type = factor(predictor.type),
    outcome.type = factor(outcome.type)
  ) %>%
  # remove any bias values exceeding -100 or 100
  filter(beta1_bias > -100 & beta1_bias < 100)

# create a function to summarize the data
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# employ function to summarize the data
plot_df_time_beta1_summary <- data_summary(plot_df_time_beta1, varname="beta1_bias", groupnames=c("T_total", "predictor.type", "outcome.type", "model")) %>%
  # turn bias values into absolute values
  mutate(beta1_bias = abs(beta1_bias))

# create lineplot (continous X and Y)
ggplot(plot_df_time_beta1_summary, aes(x = T_total, y = beta1_bias, color = model)) +
  geom_line(aes(colour = model), linewidth = 1) +
  # geom_errorbar(aes(ymin = beta1_bias - sd, ymax = beta1_bias + sd), width = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Number of timepoints", y = "Absolute Mean Estimation Error", color = "Model") +
  facet_grid(outcome.type ~ predictor.type) +  # 2x2 grid
  theme_bw() 



