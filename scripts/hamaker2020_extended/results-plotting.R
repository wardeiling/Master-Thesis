# load packages
library(purrr) # for functional programming
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation
library(stringr) # for string manipulation
library(ggplot2) # for plotting
library(ggpubr)
library(cowplot)
library(latex2exp)

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

# loop to make all plots
for(i in 1:nrow(settings)) {
  set.predictor.type <- settings$predictor.type[i]
  set.outcome.type <- settings$outcome.type[i]

  # create string for the file name
  type <- paste0("pred_", set.predictor.type, "_out_", set.outcome.type, "_")
  
  #### Plot: Grid of sd.u0 and T_total ----
  
  # read in the final data frame
  final_df <- readRDS(paste0("simulation_results_glmm/", runname, "/plotting_bias_df.RDS"))
  
  # select relevant variables and cases (select and filter)
  plot_df <- final_df %>%
    # select T = 20, N = 200, sd.u0 = 1, predictor.type = "binary" and outcome.type = "continuous"
    filter(T_total %in% c(5, 20), N_total == 200, sd.u0 %in% c(1, 3), predictor.type == set.predictor.type, outcome.type == set.outcome.type,
           sdX.between == 3, g.01 == 3) %>%
    select(-c(ends_with("_success"), ends_with("_X"), ends_with("_X.cent"), ends_with("_X.cluster.means")))
  
  plot_df_beta1 <- plot_df %>%
    select(-ends_with("_g.01_bias")) %>%
    # turn the bias variables into long format, with a new column indicating the model name
    pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "beta1_bias") %>%
    # remove the "_g.10_bias" suffix from the model names
    mutate(model = str_remove(model, "_g.10_bias")) %>%
    # remove models with a 3 in the name
    filter(!str_detect(model, "3")) %>%
    # change model names
    mutate(model = recode(model,
                          "l1" = "M1",
                          "l2" = "M2",
                          "l4" = "M3",
                          "g.independence1" = "G1.independence",
                          "g.exchangeable1" = "G1.exchangeable",
                          "g.ar11" = "G1.AR1",
                          "g.independence2" = "G2.independence",
                          "g.exchangeable2" = "G2.exchangeable",
                          "g.ar12" = "G2.AR1",
                          "g.independence4" = "G3.independence",
                          "g.exchangeable4" = "G3.exchangeable",
                          "g.ar14" = "G3.AR1")) %>%
    # set factor levels of model to ensure correct order in the plot
    mutate(model = factor(model, levels = c("M1", "G1.independence", "G1.exchangeable", "G1.AR1",
                                            "M2", "G2.independence", "G2.exchangeable", "G2.AR1",
                                            "M3", "G3.independence", "G3.exchangeable", "G3.AR1"
                                            ))) %>%
    # Turn variables into labels
    mutate(sd.u0_label = factor(sd.u0,
                                levels = c(1, 3),
                                labels = c(expression(sigma[u] == 1), expression(sigma[u] == 3)))) %>%
    mutate(T_total_label = factor(T_total,
                                  levels = c(5, 20),
                                  labels = c("T == 5", "T == 20"))) %>%
    # create new variable indicating method type (so M1 and G1 are "Method 1")
    mutate(method_type = case_when(
      str_detect(model, "M1") ~ "UC",
      str_detect(model, "M2") ~ "CWC",
      str_detect(model, "M3") ~ "MuCo",
      str_detect(model, "G1") ~ "UC",
      str_detect(model, "G2") ~ "CWC",
      str_detect(model, "G3") ~ "MuCo"
    )) %>%
    mutate(estimation_type = case_when(
      str_detect(model, "M1") ~ "GLMM",
      str_detect(model, "M2") ~ "GLMM",
      str_detect(model, "M3") ~ "GLMM",
      str_detect(model, "independence") ~ "GEE-indep",
      str_detect(model, "exchangeable") ~ "GEE-exch",
      str_detect(model, "AR1") ~ "GEE-AR1"
    )) %>%
    # set factor levels of method_type to ensure correct order in the plot
    mutate(method_type = factor(method_type, levels = c("UC", "CWC", "MuCo")),
           estimation_type = factor(estimation_type, levels = c("GLMM", "GEE-indep", "GEE-exch", "GEE-AR1"))) %>%
    # remove all bias values exceeding 100
    mutate(beta1_bias = ifelse(abs(beta1_bias) > 100, NA, beta1_bias))
    
  
  plot_df_g01 <- plot_df %>%
    select(-ends_with("_g.10_bias")) %>%
    # turn the bias variables into long format, with a new column indicating the model name
    pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "g01_bias") %>%
    # remove the "_g.01_bias" suffix from the model names
    mutate(model = str_remove(model, "_g.01_bias")) %>%
    # remove models with a 3 in the name
    filter(model == "l4" | model == "g.independence4" | model == "g.exchangeable4" | model == "g.ar14") %>%
    # change model names
    mutate(model = recode(model,
                          "l4" = "M3",
                          "g.independence4" = "G3.independence",
                          "g.exchangeable4" = "G3.exchangeable",
                          "g.ar14" = "G3.AR1")) %>%
    # set factor levels of model to ensure correct order in the plot
    mutate(model = factor(model, levels = c("M3", "G3.independence", "G3.exchangeable", "G3.AR1"
    ))) %>%
    # # create new variable indicating method type (so M1 and G1 are "Method 1")
    # mutate(method_type = case_when(
    #   str_detect(model, "M3") ~ "MuCo",
    #   str_detect(model, "G3") ~ "MuCo"
    # )) %>%
    mutate(estimation_type = case_when(
      str_detect(model, "M1") ~ "GLMM",
      str_detect(model, "M2") ~ "GLMM",
      str_detect(model, "M3") ~ "GLMM",
      str_detect(model, "independence") ~ "GEE-indep",
      str_detect(model, "exchangeable") ~ "GEE-exch",
      str_detect(model, "AR1") ~ "GEE-AR1"
    )) %>%
    # set factor levels of method_type to ensure correct order in the plot
    mutate(estimation_type = factor(estimation_type, levels = c("GLMM", "GEE-indep", "GEE-exch", "GEE-AR1"))) %>%
    # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
    mutate(sd.u0_label = factor(sd.u0,
                                levels = c(1, 3),
                                labels = c(expression(sigma[u] == 1), expression(sigma[u] == 3)))) %>%
    mutate(T_total_label = factor(T_total,
                                levels = c(5, 20),
                                labels = c("T == 5", "T == 20"))) %>%
    # remove all bias values exceeding 100
    mutate(g01_bias = ifelse(abs(g01_bias) > 100, NA, g01_bias)) 
  
  # The palette with grey:
  # cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

  # For the within-person effect
  ggplot(plot_df_beta1, aes(x = method_type, y = beta1_bias, col = estimation_type)) +
    geom_boxplot() +  # Suppress default outlier points
    # geom_boxplot(position = position_dodge(width = 0.75)) +  # align boxplots
    # geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), 
    #            alpha = 0.01, size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0 +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    scale_y_continuous(breaks = seq(-1.5, 1.5, by = 0.5)) +
    # ylim(-3, 3) +  # Set y-axis limits
    labs(x = "Method", y = "Bias") +
    facet_grid(sd.u0_label ~ T_total_label, labeller = label_parsed) + # Show T and N values in labels
    theme_bw() +
    # scale_x_discrete(breaks = waiver(), labels = new_labels) +  # <<-- overwrite x-axis labels
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # optionally rotate
    # add 2 vertical lines dividing the methods
    # geom_vline(xintercept = c(4.5, 8.5), linetype = "solid", color = "grey") +
    # remove X axis labels
    theme(# remove vertical grid lines
          panel.grid.major.x = element_blank(),
          # increase font size for grid titles
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          # increase font size for X entries
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13),
          # increase legend font size
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 13)
          # legend.position = "none"
          ) +
    # change legend title to "Estimation"
    scale_color_brewer(name = "Estimation", palette = "Spectral") 
  
  # compute mean of GEE independence with method type UC
  # mean_beta1_bias <- plot_df_beta1 %>%
  #   filter(estimation_type == "GEE-indep", method_type == "UC") %>%
  #   group_by(sd.u0, T_total) %>%
  #   summarise(mean_beta1_bias = mean(beta1_bias, na.rm = TRUE)) %>%
  #   ungroup()
  
  # save for test for main direct
  # ggsave("bias_plot_T_total-vs-sd.u0_within.pdf", width = 14, height = 8)
  
  # save
  ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_T_total-vs-sd.u0_within.pdf"), width = 9, height = 7)
  
  # For the contextual effect
  ggplot(plot_df_g01, aes(x = estimation_type, y = g01_bias, col = estimation_type)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
    coord_cartesian(ylim = c(-1.5, 1.5)) +  # Set y-axis limits
    # add tick mark at Y for every 0.5
    scale_y_continuous(breaks = seq(-1.5, 1.5, by = 0.5)) +
    labs(x = "Method", y = "Bias") +
    facet_grid(sd.u0_label ~ T_total_label, labeller = label_parsed) + # Show T and N values in labels
    theme_bw() +
    # remove X axis labels
    theme(# remove vertical grid lines
      panel.grid.major.x = element_blank(),
      # remove X tick marks
      axis.ticks.x = element_blank(),
      # increase font size for grid titles
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      # increase font size for X entries
      axis.text.x = element_text(size = 12, colour = NA),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 13),
      axis.title.x = element_text(size = 13, colour = NA),
      # remove legend
      # legend.text = element_text(size = 11),
      # legend.title = element_text(size = 13)
      legend.position = "none"
          ) +
    # change legend title to "Estimation"
    scale_color_brewer(name = "Estimation", palette = "Spectral") 
    # scale_color_manual(name = "Estimation", values = cbPalette) 
  
  # # save for test for main direct
  # ggsave("bias_plot_T_total-vs-sd.u0_contextual.pdf", width = 5, height = 7)

  # save
  ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_T_total-vs-sd.u0_contextual.pdf"), width = 3, height = 7)
  
  # combine plots native with gridextra
  # p_combined <- ggarrange(p_within, p_contextual, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")
}

  # ### PLOT 1: Grid of sdX.between and g.01 ----
  # 
  # # read in the final data frame
  # final_df <- readRDS(paste0("simulation_results_glmm/", runname, "/plotting_bias_df.RDS"))
  # 
  # # select relevant variables and cases (select and filter)
  # plot1_df <- final_df %>%
  #   # select T = 20, N = 200, sd.u0 = 1, predictor.type = "binary" and outcome.type = "continuous"
  #   filter(T_total == 20, N_total == 200, sd.u0 == 1, predictor.type == set.predictor.type, outcome.type == set.outcome.type) %>%
  #   select(-c(ends_with("_success"), ends_with("_X"), ends_with("_X.cent"), ends_with("_X.cluster.means")))
  # 
  # plot1_df_beta1 <- plot1_df %>%
  #   select(-ends_with("_g.01_bias")) %>%
  #   # turn the bias variables into long format, with a new column indicating the model name
  #   pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "beta1_bias") %>%
  #   # remove the "_g.10_bias" suffix from the model names
  #   mutate(model = str_remove(model, "_g.10_bias")) %>%
  #   # set factor levels of model to ensure correct order in the plot
  #   mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2",
  #                                           "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
  #                                           "g.ar14", "g.independence4"))) %>%
  #   # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
  #   mutate(sdX.between_str = paste0("sdX.between = ", sdX.between),
  #          g.01_str = paste0("g.01 = ", g.01)) %>%
  #   # Set factor levels to ensure correct order in the plot (0, 1, 3)
  #   mutate(sdX.between_str = factor(sdX.between_str, levels = c("sdX.between = 0", "sdX.between = 1", "sdX.between = 3")),
  #          g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")))
  # 
  # plot1_df_g01 <- plot1_df %>%
  #   select(-ends_with("_g.10_bias")) %>%
  #   # turn the bias variables into long format, with a new column indicating the model name
  #   pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "g01_bias") %>%
  #   # remove the "_g.01_bias" suffix from the model names
  #   mutate(model = str_remove(model, "_g.01_bias")) %>%
  #   # set factor levels of model to ensure correct order in the plot
  #   mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2",
  #                                           "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
  #                                           "g.ar14", "g.independence4"))) %>%
  #   # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
  #   mutate(sdX.between_str = paste0("sdX.between = ", sdX.between),
  #          g.01_str = paste0("g.01 = ", g.01)) %>%
  #   # Set factor levels to ensure correct order in the plot (0, 1, 3)
  #   mutate(sdX.between_str = factor(sdX.between_str, levels = c("sdX.between = 0", "sdX.between = 1", "sdX.between = 3")),
  #          g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")))
  # 
  # # For the within-person effect
  # ggplot(plot1_df_beta1, aes(x = model, y = beta1_bias, col = model)) +
  #   geom_boxplot() +
  #   geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
  #   ylim(-1.5, 1.5) +  # Set y-axis limits
  #   labs(x = "Generative Model", y = "Bias") +
  #   facet_grid(sdX.between_str ~ g.01_str) + # Show T and N values in labels
  #   theme_bw() +
  #   # remove X axis labels
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.title.x = element_blank())
  # 
  # # save
  # ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-sd.Xi_within.pdf"), width = 10, height = 8)
  # 
  # # For the between-person effect
  # ggplot(plot1_df_g01, aes(x = model, y = g01_bias, col = model)) +
  #   geom_boxplot() +
  #   geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
  #   ylim(-1.5, 1.5) +  # Set y-axis limits
  #   labs(x = "Generative Model", y = "Bias") +
  #   facet_grid(sdX.between_str ~ g.01_str) + # Show T and N values in labels
  #   theme_bw() +
  #   # remove X axis labels
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.title.x = element_blank())
  # 
  # # save
  # ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-sd.Xi_between.pdf"), width = 10, height = 8)
  # 
  # ### PLOT 2A: Grid of sd.u0 and g.01 ----
  # 
  # # select relevant variables and cases (select and filter)
  # plot2_df <- final_df %>%
  #   # select T = 20, N = 200, sd.u0 = 1, predictor.type = "binary" and outcome.type = "continuous"
  #   filter(T_total == 20, N_total == 200, sdX.between == 1, predictor.type == set.predictor.type, outcome.type == set.outcome.type) %>%
  #   select(-c(ends_with("_success"), ends_with("_X"), ends_with("_X.cent"), ends_with("_X.cluster.means")))
  # 
  # plot2_df_beta1 <- plot2_df %>%
  #   select(-ends_with("_g.01_bias")) %>%
  #   # turn the bias variables into long format, with a new column indicating the model name
  #   pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "beta1_bias") %>%
  #   # remove the "_g.10_bias" suffix from the model names
  #   mutate(model = str_remove(model, "_g.10_bias")) %>%
  #   # set factor levels of model to ensure correct order in the plot
  #   mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2",
  #                                           "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
  #                                           "g.ar14", "g.independence4"))) %>%
  #   # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
  #   mutate(g.01_str = paste0("g.01 = ", g.01),
  #          sd.u0_str = paste0("sd.u0 = ", sd.u0)) %>%
  #   # Set factor levels to ensure correct order in the plot (0, 1, 3)
  #   mutate(g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")),
  #          sd.u0_str = factor(sd.u0_str, levels = c("sd.u0 = 0", "sd.u0 = 1", "sd.u0 = 3")))
  # 
  # plot2_df_g01 <- plot2_df %>%
  #   select(-ends_with("_g.10_bias")) %>%
  #   # turn the bias variables into long format, with a new column indicating the model name
  #   pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "g01_bias") %>%
  #   # remove the "_g.01_bias" suffix from the model names
  #   mutate(model = str_remove(model, "_g.01_bias")) %>%
  #   # set factor levels of model to ensure correct order in the plot
  #   mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2",
  #                                           "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
  #                                           "g.ar14", "g.independence4"))) %>%
  #   # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
  #   mutate(g.01_str = paste0("g.01 = ", g.01),
  #          sd.u0_str = paste0("sd.u0 = ", sd.u0)) %>%
  #   # Set factor levels to ensure correct order in the plot (0, 1, 3)
  #   mutate(g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")),
  #          sd.u0_str = factor(sd.u0_str, levels = c("sd.u0 = 0", "sd.u0 = 1", "sd.u0 = 3")))
  # 
  # # For the within-person effect
  # ggplot(plot2_df_beta1, aes(x = model, y = beta1_bias, col = model)) +
  #   geom_boxplot() +
  #   geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
  #   ylim(-1.5, 1.5) +  # Set y-axis limits
  #   labs(x = "Generative Model", y = "Bias") +
  #   facet_grid(sd.u0_str ~ g.01_str) + # Show T and N values in labels
  #   theme_bw() +
  #   # remove X axis labels
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.title.x = element_blank())
  # 
  # # save
  # ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-sd.u0_within.pdf"), width = 10, height = 8)
  # 
  # # For the between-person effect
  # ggplot(plot2_df_g01, aes(x = model, y = g01_bias, col = model)) +
  #   geom_boxplot() +
  #   geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
  #   ylim(-1.5, 1.5) +  # Set y-axis limits
  #   labs(x = "Generative Model", y = "Bias") +
  #   facet_grid(sd.u0_str ~ g.01_str) + # Show T and N values in labels
  #   theme_bw() +
  #   # remove X axis labels
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.title.x = element_blank())
  # 
  # # save
  # ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-sd.u0_between.pdf"), width = 10, height = 8)
  # 
  # ### PLOT 3: Grid of T_total and g.01 ----
  # 
  # # select relevant variables and cases (select and filter)
  # plot3_df <- final_df %>%
  #   # select T = 20, N = 200, sd.u0 = 1, predictor.type = "binary" and outcome.type = "continuous"
  #   filter(N_total == 200, sdX.between == 3, sd.u0 == 1, predictor.type == set.predictor.type, outcome.type == set.outcome.type) %>%
  #   select(-c(ends_with("_success"), ends_with("_X"), ends_with("_X.cent"), ends_with("_X.cluster.means")))
  # 
  # plot3_df_beta1 <- plot3_df %>%
  #   select(-ends_with("_g.01_bias")) %>%
  #   # turn the bias variables into long format, with a new column indicating the model name
  #   pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "beta1_bias") %>%
  #   # remove the "_g.10_bias" suffix from the model names
  #   mutate(model = str_remove(model, "_g.10_bias")) %>%
  #   # set factor levels of model to ensure correct order in the plot
  #   mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2",
  #                                           "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
  #                                           "g.ar14", "g.independence4"))) %>%
  #   # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
  #   mutate(g.01_str = paste0("g.01 = ", g.01),
  #          T_total_str = paste0("T = ", T_total)) %>%
  #   # Set factor levels to ensure correct order in the plot (0, 1, 3)
  #   mutate(g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")),
  #          T_total_str = factor(T_total_str, levels = c("T = 5", "T = 10", "T = 20")))
  # 
  # plot3_df_g01 <- plot3_df %>%
  #   select(-ends_with("_g.10_bias")) %>%
  #   # turn the bias variables into long format, with a new column indicating the model name
  #   pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "g01_bias") %>%
  #   # remove the "_g.01_bias" suffix from the model names
  #   mutate(model = str_remove(model, "_g.01_bias")) %>%
  #   # set factor levels of model to ensure correct order in the plot
  #   mutate(model = factor(model, levels = c("l2", "l3a", "l4", "g.exchangeable2", "g.ar12", "g.independence2",
  #                                           "g.exchangeable3", "g.ar13", "g.independence3", "g.exchangeable4",
  #                                           "g.ar14", "g.independence4"))) %>%
  #   # Turn label variables (sdX.between, g.01 and sd.u0) into strings with an underscore
  #   mutate(g.01_str = paste0("g.01 = ", g.01),
  #          T_total_str = paste0("T = ", T_total)) %>%
  #   # Set factor levels to ensure correct order in the plot (0, 1, 3)
  #   mutate(g.01_str = factor(g.01_str, levels = c("g.01 = 0", "g.01 = 1", "g.01 = 3")),
  #          T_total_str = factor(T_total_str, levels = c("T = 5", "T = 10", "T = 20")))
  # 
  # # For the within-person effect
  # ggplot(plot3_df_beta1, aes(x = model, y = beta1_bias, col = model)) +
  #   geom_boxplot() +
  #   geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
  #   ylim(-1.5, 1.5) +  # Set y-axis limits
  #   labs(x = "Generative Model", y = "Bias") +
  #   facet_grid(T_total_str ~ g.01_str) + # Show T and N values in labels
  #   theme_bw() +
  #   # remove X axis labels
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.title.x = element_blank())
  # 
  # # save
  # ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-T_total_within.pdf"), width = 10, height = 8)
  # 
  # # For the between-person effect
  # ggplot(plot3_df_g01, aes(x = model, y = g01_bias, col = model)) +
  #   geom_boxplot() +
  #   geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
  #   ylim(-1.5, 1.5) +  # Set y-axis limits
  #   labs(x = "Generative Model", y = "Bias") +
  #   facet_grid(T_total_str ~ g.01_str) + # Show T and N values in labels
  #   theme_bw() +
  #   # remove X axis labels
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.title.x = element_blank())
  # 
  # # save
  # ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_g01-vs-T_total_between.pdf"), width = 10, height = 8)
  # 

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

# }

# ### PLOT 3 ----
# 
# runname <- "April18_fullsimulation_combined"
# 
# # read in the final data frame
# final_df <- readRDS(paste0("simulation_results_glmm/", runname, "/plotting_bias_df.RDS"))
# 
# # Prepare data for plotting beta1 bias across timepoints
# plot_df_time_beta1 <- final_df %>%
#   filter(N_total == 200, sdX.between == 3, sd.u0 == 1, g.01 == 1) %>%
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
# plot_df_time_beta1_summary <- data_summary(plot_df_time_beta1, varname="beta1_bias", groupnames=c("T_total", "predictor.type", "outcome.type", "model")) %>%
#   # turn bias values into absolute values
#   mutate(beta1_bias = abs(beta1_bias))
# 
# # create lineplot (continous X and Y)
# ggplot(plot_df_time_beta1_summary, aes(x = T_total, y = beta1_bias, color = model)) +
#   geom_line(aes(colour = model), linewidth = 1) +
#   # geom_errorbar(aes(ymin = beta1_bias - sd, ymax = beta1_bias + sd), width = 1) +
#   geom_point(size = 2) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Number of timepoints", y = "Absolute Mean Estimation Error", color = "Model") +
#   facet_grid(outcome.type ~ predictor.type) +  # 2x2 grid
#   theme_bw() 

### WITHIN-PERSON EFFECTS PLOTS ### ----

runname <- "April18_fullsimulation_combined"

# read in the final data frame
final_df <- readRDS(paste0("simulation_results_glmm/", runname, "/plotting_bias_df.RDS"))

# scenario 1
general_scenario_df_beta1 <- final_df %>%
  # remove columns ending with "g.01", "X", "X.cent", "X.cluster.means"
  select(-ends_with("_g.01_bias"), -ends_with("_X"), -ends_with("_X.cent"), -ends_with("_X.cluster.means")) %>%
  pivot_longer(cols = ends_with("_g.10_bias"),
               names_to = "model", values_to = "beta1_bias") %>%
  mutate(
    model = str_remove(model, "_g.10_bias"),
    model = factor(model, levels = c("l1", "l2", "l3a", "l4", 
                                     "g.exchangeable1", "g.ar11", "g.independence1",
                                     "g.exchangeable2", "g.ar12", "g.independence2", 
                                     "g.exchangeable3", "g.ar13", "g.independence3", 
                                     "g.exchangeable4", "g.ar14", "g.independence4")),
    
    # Create labels and set factor levels of model to ensure correct order in the plot
    predictor.type_str = paste0(predictor.type, " predictor"),
    predictor.type_str = factor(predictor.type_str, levels = c("continuous predictor", "binary predictor")),
    outcome.type_str = paste0(outcome.type, " outcome"),
    outcome.type_str = factor(outcome.type_str, levels = c("continuous outcome", "binary outcome")),
    T_total_str = paste0("T = ", T_total),
    T_total_str = factor(T_total_str, levels = c("T = 5", "T = 10", "T = 20")),
    sd.u0_str = paste0("sd.u0 = ", sd.u0),
    sd.u0_str = factor(sd.u0_str, levels = c("sd.u0 = 0", "sd.u0 = 1", "sd.u0 = 3"))
  ) %>%
  # remove parametrization 3
  # filter(!model %in% c("l3a", "g.independence3", "g.ar13", "g.exchangeable3")) %>%
  # remove any bias values exceeding -100 or 100 (occuring in the GEE models)
  filter(beta1_bias > -100 & beta1_bias < 100)

# temp
# temp <- general_scenario_df_beta1 %>%
#   filter(model == "g.independence1", N_total == 200, T_total == 10, sdX.between == 3, sd.u0 == 0, g.01 == 3)

# create 2x2 grid of boxplots
boxplot_grid_maker_beta1 <- function(df) {
  ggplot(df, aes(x = model, y = beta1_bias, col = model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
    ylim(-1.5, 1.5) +  # Set y-axis limits
    labs(x = "Generative Model", y = "Bias") +
    facet_grid(predictor.type_str ~ outcome.type_str) + 
    theme_bw() + 
    # remove X axis labels
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
}

# scenario 01: baseline (all should perform very well, as there is no (random nor systematic) heterogeneity in outcome)
scenario01_df_beta1 <- general_scenario_df_beta1 %>% 
  filter(N_total == 200, T_total == 20, sdX.between == 0, sd.u0 == 0, g.01 == 0)
boxplot_grid_maker_beta1(scenario01_df_beta1)
ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_within_scenario01.pdf"), width = 10, height = 5)

# scenario 02: baseline (with between-person differences in X)
scenario02_df_beta1 <- general_scenario_df_beta1 %>%
  filter(N_total == 200, T_total == 20, sdX.between == 1, sd.u0 == 0, g.01 == 0)
boxplot_grid_maker_beta1(scenario02_df_beta1)
ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_within_scenario02.pdf"), width = 10, height = 5)

# scenario 1: challenging for parametrization 1, but no random effect, which is good for GEE
scenario1_df_beta1 <- general_scenario_df_beta1 %>%
  filter(N_total == 200, T_total == 5, sdX.between == 3, sd.u0 == 0, g.01 == 3)
boxplot_grid_maker_beta1(scenario1_df_beta1) # g.independence1 is out of bounds for cont XY
ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_within_scenario1.pdf"), width = 10, height = 5)

# scenario 2: challenging for parametrization 1 and GEE (large random effect)
scenario2_df_beta1 <- general_scenario_df_beta1 %>%
  filter(N_total == 200, T_total == 5, sdX.between == 3, sd.u0 == 3, g.01 == 3)
boxplot_grid_maker_beta1(scenario2_df_beta1) # g.independence1 is out of bounds for cont XY
ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_within_scenario2.pdf"), width = 10, height = 5)

# create 2x2x2 grid of boxplots
boxplot_grid_maker_beta1_3vars <- function(df) {
  ggplot(df, aes(x = model, y = beta1_bias, col = model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
    ylim(-1.5, 1.5) +  # Set y-axis limits
    labs(x = "Generative Model", y = "Bias") +
    facet_grid(predictor.type_str ~ outcome.type_str + T_total_str) + 
    theme_bw() + 
    # remove X axis labels
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
}

# scenario 3: more comprehensive plot with sd.u0
scenario3_df_beta1 <- general_scenario_df_beta1 %>%
  filter(N_total == 200, T_total == 5, sdX.between == 3, sd.u0 %in% c(0, 3), g.01 == 3)

ggplot(scenario3_df_beta1, aes(x = model, y = beta1_bias, col = model)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
  ylim(-1.5, 1.5) +  # Set y-axis limits
  labs(x = "Generative Model", y = "Bias") +
  facet_grid(predictor.type_str ~ outcome.type_str + sd.u0_str) + 
  theme_bw() + 
  # remove X axis labels
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        # put legend under
        legend.position = "bottom")

ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_within_scenario3.pdf"), width = 10, height = 10)

# scenario 4: more comprehensive plot with T_total
scenario4_df_beta1 <- general_scenario_df_beta1 %>%
  filter(N_total == 200, T_total %in% c(5, 20), sdX.between == 3, sd.u0 == 0, g.01 == 3)

ggplot(scenario4_df_beta1, aes(x = model, y = beta1_bias, col = model)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
  ylim(-1.5, 1.5) +  # Set y-axis limits
  labs(x = "Generative Model", y = "Bias") +
  facet_grid(T_total_str + predictor.type_str ~ outcome.type_str) + 
  theme_bw() + 
  # remove X axis labels
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        # put legend under
        legend.position = "bottom")

# ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_within_scenario4.pdf"), width = 10, height = 8)
  
### CONTEXTUAL EFFECTS PLOTS ### ----

runname <- "April18_fullsimulation_combined"

# read in the final data frame
final_df <- readRDS(paste0("simulation_results_glmm/", runname, "/plotting_bias_df.RDS"))

# scenario 1
general_scenario_df_g01 <- final_df %>%
  # remove columns ending with "g.01", "X", "X.cent", "X.cluster.means"
  select(-ends_with("_g.10_bias"), -ends_with("_X"), -ends_with("_X.cent"), -ends_with("_X.cluster.means")) %>%
  pivot_longer(cols = ends_with("_g.01_bias"),
               names_to = "model", values_to = "g01_bias") %>%
  mutate(
    model = str_remove(model, "_g.01_bias"),
    model = factor(model, levels = c("l1", "l2", "l3a", "l4", 
                                     "g.exchangeable1", "g.ar11", "g.independence1",
                                     "g.exchangeable2", "g.ar12", "g.independence2", 
                                     "g.exchangeable3", "g.ar13", "g.independence3", 
                                     "g.exchangeable4", "g.ar14", "g.independence4")),
    # Create labels and set factor levels of model to ensure correct order in the plot
    predictor.type_str = paste0(predictor.type, " predictor"),
    predictor.type_str = factor(predictor.type_str, levels = c("continuous predictor", "binary predictor")),
    outcome.type_str = paste0(outcome.type, " outcome"),
    outcome.type_str = factor(outcome.type_str, levels = c("continuous outcome", "binary outcome")),
    T_total_str = paste0("T = ", T_total),
    T_total_str = factor(T_total_str, levels = c("T = 5", "T = 10", "T = 20")),
    sd.u0_str = paste0("sd.u0 = ", sd.u0),
    sd.u0_str = factor(sd.u0_str, levels = c("sd.u0 = 0", "sd.u0 = 1", "sd.u0 = 3"))
  ) %>%
  # remove parametrization 3
  # filter(!model %in% c("l3a", "g.independence3", "g.ar13", "g.exchangeable3")) %>%
  # remove any bias values exceeding -100 or 100 (occuring in the GEE models)
  filter(g01_bias > -100 & g01_bias < 100)

# create 2x2 grid of boxplots
boxplot_grid_maker_g01 <- function(df) {
  ggplot(df, aes(x = model, y = g01_bias, col = model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
    ylim(-1.5, 1.5) +  # Set y-axis limits
    labs(x = "Generative Model", y = "Bias") +
    facet_grid(predictor.type_str ~ outcome.type_str) + 
    theme_bw() + 
    # remove X axis labels
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
}

# scenario 01: baseline (all should perform very well, as there is no (random nor systematic) heterogeneity in outcome)
scenario01_df_g01 <- general_scenario_df_g01 %>% 
  filter(N_total == 200, T_total == 20, sdX.between == 0, sd.u0 == 0, g.01 == 0)
boxplot_grid_maker_g01(scenario01_df_g01)
ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_contextual_scenario01.pdf"), width = 10, height = 5)

# scenario 02: baseline (with between-person differences in X)
scenario02_df_g01 <- general_scenario_df_g01 %>%
  filter(N_total == 200, T_total == 20, sdX.between == 1, sd.u0 == 0, g.01 == 0)
boxplot_grid_maker_g01(scenario02_df_g01)
ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_contextual_scenario02.pdf"), width = 10, height = 5)

# scenario 1: challenging for parametrization 1, but no random effect, which is good for GEE
scenario1_df_g01 <- general_scenario_df_g01 %>%
  filter(N_total == 200, T_total == 5, sdX.between == 3, sd.u0 == 0, g.01 == 3)
boxplot_grid_maker_g01(scenario1_df_g01) # g.independence1 is out of bounds for cont XY
ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_contextual_scenario1.pdf"), width = 10, height = 5)

# scenario 2: challenging for parametrization 1 and GEE (large random effect)
scenario2_df_g01 <- general_scenario_df_g01 %>%
  filter(N_total == 200, T_total == 5, sdX.between == 3, sd.u0 == 3, g.01 == 3)
boxplot_grid_maker_g01(scenario2_df_g01) # g.independence1 is out of bounds for cont XY
ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_contextual_scenario2.pdf"), width = 10, height = 5)

# # scenario 2b (with large T)
# scenario2b_df_g01 <- general_scenario_df_g01 %>%
#   filter(N_total == 200, T_total == 20, sdX.between == 3, sd.u0 == 0, g.01 == 3)
# boxplot_grid_maker_g01(scenario2b_df_g01) # g.independence1 is out of bounds for cont XY                  

# scenario 3: more comprehensive plot with sd.u0
scenario3_df_g01 <- general_scenario_df_g01 %>%
  filter(N_total == 200, T_total == 5, sdX.between == 3, sd.u0 %in% c(0, 3), g.01 == 3)
ggplot(scenario3_df_g01, aes(x = model, y = g01_bias, col = model)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0
  ylim(-1.5, 1.5) +  # Set y-axis limits
  labs(x = "Generative Model", y = "Bias") +
  facet_grid(predictor.type_str ~ outcome.type_str + sd.u0_str) + 
  theme_bw() + 
  # remove X axis labels
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        # put legend under
        legend.position = "bottom")
ggsave(paste0("simulation_results_glmm/", runname, "/figures/bias_plot_contextual_scenario3.pdf"), width = 10, height = 10)


# checking how many GEEs yielded extreme values
runname <- "April18_fullsimulation_combined"

# read in the final data frame
final_df <- readRDS(paste0("simulation_results_glmm/", runname, "/plotting_bias_df.RDS"))

# scenario 1
GEE_check_df <- final_df %>%
  # remove columns ending with "g.01", "X", "X.cent", "X.cluster.means"
  select(-ends_with("_g.01_bias"), -ends_with("_X"), -ends_with("_X.cent"), -ends_with("_X.cluster.means")) %>%
  pivot_longer(cols = ends_with("_g.10_bias"),
               names_to = "model", values_to = "beta1_bias") %>%
  mutate(
    model = str_remove(model, "_g.10_bias"),
    model = factor(model, levels = c("l1", "l2", "l3a", "l4", 
                                     "g.exchangeable1", "g.ar11", "g.independence1",
                                     "g.exchangeable2", "g.ar12", "g.independence2", 
                                     "g.exchangeable3", "g.ar13", "g.independence3", 
                                     "g.exchangeable4", "g.ar14", "g.independence4"))) %>%
  # select only the GEE models
  filter(model %in% c("g.exchangeable1", "g.ar11", "g.independence1",
                      "g.exchangeable2", "g.ar12", "g.independence2", 
                      "g.exchangeable3", "g.ar13", "g.independence3", 
                      "g.exchangeable4", "g.ar14", "g.independence4")) 


# define a threshold for "extreme" bias
threshold <- 1E+11  # can also be set to 5, doesn't affect the number of extreme values

extreme_prop_df <- GEE_check_df %>%
  group_by(model, design_id) %>%
  summarize(
    proportion_extreme = mean(abs(beta1_bias) > threshold),
    .groups = "drop"
  ) %>%
  # remove rows with 0 proportion extreme
  filter(proportion_extreme > 0)

# count number of unique design_ids
length(unique(extreme_prop_df$design_id))

# compute min and max
min(extreme_prop_df$proportion_extreme)
max(extreme_prop_df$proportion_extreme)

# plot a density plot of extreme values wrapped by model
ggplot(extreme_prop_df, aes(x = proportion_extreme, fill = model)) +
  geom_density(alpha = 0.5) +
  labs(x = "Proportion of extreme values", y = "Density") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~ model, ncol = 3)

# create a table with the number of extreme scenarios per model
extreme_count_df <- GEE_check_df %>%
  group_by(model) %>%
  summarize(
    num_extreme = sum(abs(beta1_bias) > threshold),
    .groups = "drop"
  ) %>%
  arrange(desc(num_extreme))

# take the setting with most frequent extreme values (design_id = 299 and g.ar12) and plot the density of bias
plot1 <- GEE_check_df %>%
  filter(design_id == 299, model == "g.ar12") %>%
  ggplot(aes(x = beta1_bias)) +
  geom_density() +
  labs(x = "Bias", y = "Density") +
  theme_bw() +
  xlim(-5, 5) +
  theme(legend.position = "bottom")

# g.ar12
# 300
plot2 <- GEE_check_df %>%
  filter(design_id == 300, model == "g.ar12") %>%
  ggplot(aes(x = beta1_bias)) +
  geom_density() +
  labs(x = "Bias", y = "Density") +
  theme_bw() +
  xlim(-5, 5) +
  theme(legend.position = "bottom")

# g.ar12
# 353
plot3 <- GEE_check_df %>%
  filter(design_id == 353, model == "g.ar12") %>%
  ggplot(aes(x = beta1_bias)) +
  geom_density() +
  labs(x = "Bias", y = "Density") +
  theme_bw() +
  xlim(-5, 5) +
  theme(legend.position = "bottom")

# g.ar12
# 354
plot4 <- GEE_check_df %>%
  filter(design_id == 354, model == "g.ar12") %>%
  ggplot(aes(x = beta1_bias)) +
  geom_density() +
  labs(x = "Bias", y = "Density") +
  theme_bw() +
  xlim(-5, 5) +
  theme(legend.position = "bottom")

# g.ar13
# 263
plot5 <- GEE_check_df %>%
  filter(design_id == 263, model == "g.ar13") %>%
  ggplot(aes(x = beta1_bias)) +
  geom_density() +
  labs(x = "Bias", y = "Density") +
  theme_bw() +
  xlim(-5, 5) +
  theme(legend.position = "bottom")

# g.ar14
# 263
plot6 <- GEE_check_df %>%
  filter(design_id == 263, model == "g.ar14") %>%
  ggplot(aes(x = beta1_bias)) +
  geom_density() +
  labs(x = "Bias", y = "Density") +
  theme_bw() +
  xlim(-5, 5) +
  theme(legend.position = "bottom")

# combine plots
library(patchwork)
combined_plot <- (plot1 + plot2) / (plot3 + plot4) / (plot5 + plot6)
combined_plot


plot_df_beta1 <- plotting_bias_df %>%
  select(-ends_with("_g.01_bias")) %>%
  # turn the bias variables into long format, with a new column indicating the model name
  pivot_longer(cols = ends_with("_bias"), names_to = "model", values_to = "beta1_bias") %>%
  # remove the "_g.10_bias" suffix from the model names
  mutate(model = str_remove(model, "_g.10_bias")) %>%
  # remove models with a 3 in the name
  filter(!str_detect(model, "3")) %>%
  # change model names
  mutate(model = recode(model,
                        "l1" = "M1",
                        "l2" = "M2",
                        "l4" = "M3",
                        "g.independence1" = "G1.independence",
                        "g.exchangeable1" = "G1.exchangeable",
                        "g.ar11" = "G1.AR1",
                        "g.independence2" = "G2.independence",
                        "g.exchangeable2" = "G2.exchangeable",
                        "g.ar12" = "G2.AR1",
                        "g.independence4" = "G3.independence",
                        "g.exchangeable4" = "G3.exchangeable",
                        "g.ar14" = "G3.AR1")) %>%
  # set factor levels of model to ensure correct order in the plot
  mutate(model = factor(model, levels = c("M1", "G1.independence", "G1.exchangeable", "G1.AR1",
                                          "M2", "G2.independence", "G2.exchangeable", "G2.AR1",
                                          "M3", "G3.independence", "G3.exchangeable", "G3.AR1"
  ))) %>%
  # Turn variables into labels
  mutate(sd.u0_label = factor(sd.u0,
                              levels = c(1, 3),
                              labels = c(expression(sigma[u] == 1), expression(sigma[u] == 3)))) %>%
  mutate(T_total_label = factor(T_total,
                                levels = c(5, 20),
                                labels = c("T == 5", "T == 20"))) %>%
  # create new variable indicating method type (so M1 and G1 are "Method 1")
  mutate(method_type = case_when(
    str_detect(model, "M1") ~ "UC",
    str_detect(model, "M2") ~ "CWC",
    str_detect(model, "M3") ~ "MuCo",
    str_detect(model, "G1") ~ "UC",
    str_detect(model, "G2") ~ "CWC",
    str_detect(model, "G3") ~ "MuCo"
  )) %>%
  mutate(estimation_type = case_when(
    str_detect(model, "M1") ~ "GLMM",
    str_detect(model, "M2") ~ "GLMM",
    str_detect(model, "M3") ~ "GLMM",
    str_detect(model, "independence") ~ "GEE-indep",
    str_detect(model, "exchangeable") ~ "GEE-exch",
    str_detect(model, "AR1") ~ "GEE-AR1"
  )) %>%
  # set factor levels of method_type to ensure correct order in the plot
  mutate(method_type = factor(method_type, levels = c("UC", "CWC", "MuCo")),
         estimation_type = factor(estimation_type, levels = c("GLMM", "GEE-indep", "GEE-exch", "GEE-AR1"))) %>%
  # select cases with sigma_u = 3 nd T = 20
  filter(sd.u0 == 3, T_total == 5, model == "G1.independence")

# create density of the estimates of beta1
ggplot(plot_df_beta1, aes(x = beta1_bias)) +
  geom_density(aes(fill = model), alpha = 0.5) +
  labs(x = "Bias", y = "Density") +
  theme_bw() 

# create normal histogram
ggplot(plot_df_beta1, aes(x = beta1_bias)) +
  geom_histogram(aes(fill = model), alpha = 0.5, bins = 30) +
  labs(x = "Bias", y = "Density") +
  theme_bw() +
  facet_wrap(~ model, ncol = 3)
  
summary(plot_df_beta1$beta1_bias)
