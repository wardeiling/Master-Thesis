# load packages
library(purrr) # for functional programming
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation
library(stringr) # for string manipulation
library(ggplot2) # for plotting

### RETRIEVE RESULTS OF SIMULATIONS ### ----

### PART 1: DGM 2 (binary X and continuous Y), DGM 3 (continuous X and binary Y) and DGM 4 (binary X and Y) ----
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

### PART 2: DGM 1 (continuous X and Y) ----

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
      str_detect(model, "independence") ~ "GEE-independence",
      str_detect(model, "exchangeable") ~ "GEE-exchangeable",
      str_detect(model, "AR1") ~ "GEE-AR(1)"
    )) %>%
    # set factor levels of method_type to ensure correct order in the plot
    mutate(method_type = factor(method_type, levels = c("UC", "CWC", "MuCo")),
           estimation_type = factor(estimation_type, levels = c("GLMM", "GEE-independence", "GEE-exchangeable", "GEE-AR(1)"))) %>%
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
    mutate(estimation_type = case_when(
      str_detect(model, "M1") ~ "GLMM",
      str_detect(model, "M2") ~ "GLMM",
      str_detect(model, "M3") ~ "GLMM",
      str_detect(model, "independence") ~ "GEE-independence",
      str_detect(model, "exchangeable") ~ "GEE-exchangeable",
      str_detect(model, "AR1") ~ "GEE-AR(1)"
    )) %>%
    # set factor levels of method_type to ensure correct order in the plot
    mutate(estimation_type = factor(estimation_type, levels = c("GLMM", "GEE-independence", "GEE-exchangeable", "GEE-AR(1)"))) %>%
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
    geom_hline(yintercept = 0, linetype = "dashed") +  # Dashed horizontal line at 0 +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    scale_y_continuous(breaks = seq(-1.5, 1.5, by = 0.5)) +
    # ylim(-3, 3) +  # Set y-axis limits
    labs(x = "Method", y = "Bias") +
    facet_grid(sd.u0_label ~ T_total_label, labeller = label_parsed) + # Show T and N values in labels
    theme_bw() +
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
      legend.title = element_text(size = 13),
      legend.position = "bottom",
      legend.box.background = element_rect(color="black", size=2)
    ) +
    # change legend title to "Estimation"
    scale_color_brewer(name = "Estimation", palette = "Spectral") 
  
  # save for test for main direct
  ggsave("bias_plot_T_total-vs-sd.u0_within.pdf", width = 9, height = 7)
  
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
      # increase legend font size
      legend.text = element_text(size = 11, colour = NA),
      legend.title = element_text(size = 13, colour = NA),
      legend.position = "bottom",
    ) +
    # change legend title to "Estimation"
    scale_color_brewer(name = "Estimation", palette = "Spectral") +
    guides(color = guide_legend(override.aes = list(color = NA)))
  
  # save
  ggsave(paste0("simulation_results_glmm/", runname, "/figures/", type, "bias_plot_T_total-vs-sd.u0_contextual.pdf"), width = 3, height = 7)

}
