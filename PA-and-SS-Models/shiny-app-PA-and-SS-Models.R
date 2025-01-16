# Load necessary libraries
library(shiny)

# Define UI for application
ui <- fluidPage(
  titlePanel("Population-Averaged vs. Subject-Specific Effects: GLMM and MLM"),
  
  tabsetPanel(
    # Tab for Generalized Mixed Linear Model with Log-Link
    tabPanel("GLMM with Log-Link and Continuous Predictor",
             sidebarLayout(
               sidebarPanel(
                 h4("Adjust Parameters for GLMM"),
                 sliderInput("alpha_0", "Beta_0 (Intercept):", 
                             min = -5, max = 5, value = 0, step = 0.1),
                 sliderInput("alpha_1", "Beta_1 (Slope):", 
                             min = -5, max = 5, value = 0.5, step = 0.1),
                 sliderInput("sigma_b0", "Sigma_b0 (Random Intercept SD):", 
                             min = 0.0, max = 5, value = 0.5, step = 0.1)
               ),
               mainPanel(
                 h3("GLMM with Log-Link and Continuous Predictor"),
                 plotOutput("glmmPlot"),
                 p("This plot illustrates the difference between the marginal (population-averaged) and conditional (subject-specific) relationships in a generalized mixed linear model (GLMM)
                   with a log-link function. Adjust the sliders to see how the parameters influence the curves.")
               )
             )),
    
    # Tab for GLMM with Binary Predictor
    tabPanel("GLMM with Log-Link and Binary Predictor",
             sidebarLayout(
               sidebarPanel(
                 h4("Adjust Parameters for GLMM"),
                 sliderInput("binary_alpha_0", "Beta_0 (Intercept):", 
                             min = -5, max = 5, value = 0, step = 0.1),
                 sliderInput("binary_alpha_1", "Beta_1 (Slope):", 
                             min = -5, max = 5, value = 5, step = 0.1),
                 sliderInput("binary_sigma_b0", "Sigma_b0 (Random Intercept SD):", 
                             min = 0.0, max = 5, value = 2, step = 0.1),
                 checkboxInput("show_interpolation", "Show Logistic Curve Interpolation", TRUE)
               ),
               mainPanel(
                 h3("GLMM with Log-Link and Binary Predictor"),
                 plotOutput("binaryGLMMPlot"),
                 p("This visualization highlights the difference between population-averaged and subject-specific probabilities in a GLMM with a binary predictor. Use the sliders to adjust parameters 
                   and toggle interpolation to observe logistic curve overlays.")
               )
             )),
    
    # Tab for Mixed Linear Model
    tabPanel("Mixed Linear Model (MLM)",
             sidebarLayout(
               sidebarPanel(
                 h4("Under Development"),
                 p("This tab will include an interactive plot and controls for the Mixed Linear Model.")
               ),
               mainPanel(
                 h3("Mixed Linear Model Placeholder"),
                 p("This tab is a placeholder for the Mixed Linear Model (MLM) visualization. You can add content here as needed.")
               )
             ))
  )
)

# Define server logic
server <- function(input, output) {
  
  # Plot for GLMM with Log-Link
  output$glmmPlot <- renderPlot({
    alpha_0 <- input$alpha_0
    alpha_1 <- input$alpha_1
    sigma_b0 <- input$sigma_b0
    
    set.seed(123)
    n_individuals <- 20
    X_grid <- seq(-5, 5, length.out = 10000)
    b0 <- rnorm(n_individuals, mean = 0, sd = sigma_b0)
    
    logit_conditional <- alpha_0 + alpha_1 * X_grid
    pi_conditional <- 1 / (1 + exp(-logit_conditional))
    
    logit_individuals <- sapply(b0, function(b) alpha_0 + alpha_1 * X_grid + b)
    pi_individuals <- 1 / (1 + exp(-logit_individuals))
    pi_population <- rowMeans(pi_individuals)
    
    plot(X_grid, pi_conditional, type = "l", lwd = 4, col = "black",
         ylim = c(0, 1), xlim = c(-5, 5), xlab = expression(X[it]), ylab = expression(pi),
         main = "GLMM with Continuous Predictor")
    for (i in 1:n_individuals) {
      lines(X_grid, pi_individuals[, i], col = "grey", lwd = 1)
    }
    lines(X_grid, pi_population, col = "red", lwd = 4)
    legend("bottomright", legend = c("conditional", "population", "individuals"),
           col = c("black", "red", "grey"), lty = c(1, 1, 1), lwd = c(4, 4, 1))
  })
  
  # Plot for GLMM with Binary Predictor
  output$binaryGLMMPlot <- renderPlot({
    alpha_0 <- input$binary_alpha_0
    alpha_1 <- input$binary_alpha_1
    sigma_b0 <- input$binary_sigma_b0
    show_interpolation <- input$show_interpolation
    
    set.seed(123)
    n_individuals <- 20
    X_grid_cont <- seq(0, 1, length.out = 10000)
    X_grid_bin <- c(0, 1)
    b0 <- rnorm(n_individuals, mean = 0, sd = sigma_b0)
    
    logit_individuals_bin <- sapply(b0, function(b) outer(alpha_0 + alpha_1 * X_grid_bin, b, "+"))
    pi_individuals_bin <- 1 / (1 + exp(-logit_individuals_bin))
    pi_population_bin <- rowMeans(pi_individuals_bin)
    
    plot(X_grid_cont, numeric(length(X_grid_cont)), type = "n", ylim = c(0, 1), xlim = c(0, 1),
         xlab = expression(X[it]), ylab = expression(pi), main = "GLMM with Binary Predictor")
    
    if (show_interpolation) {
      logit_conditional <- alpha_0 + alpha_1 * X_grid_cont
      pi_conditional <- 1 / (1 + exp(-logit_conditional))
      logit_individuals_cont <- sapply(b0, function(b) alpha_0 + alpha_1 * X_grid_cont + b)
      pi_individuals_cont <- 1 / (1 + exp(-logit_individuals_cont))
      pi_population_cont <- rowMeans(pi_individuals_cont)
      
      lines(X_grid_cont, pi_conditional, col = "black", lwd = 4)
      for (i in 1:n_individuals) {
        lines(X_grid_cont, pi_individuals_cont[, i], col = "grey", lwd = 1)
      }
      lines(X_grid_cont, pi_population_cont, col = "red", lwd = 4)
    }
    
    points(X_grid_bin, 1 / (1 + exp(-(alpha_0 + alpha_1 * X_grid_bin))), col = "black", pch = 16, cex = 1.5)
    for (i in 1:n_individuals) {
      points(X_grid_bin, pi_individuals_bin[, i], col = "grey", pch = 16, cex = 0.6)
    }
    points(X_grid_bin, pi_population_bin, col = "red", pch = 16, cex = 1.5)
    
    # Add a legend
    if (show_interpolation) {
      
      legend("bottomright", legend = c("conditional (interpolation)", "population (interpolation)", "individuals (interpolation)", 
                                   "datapoint conditional", "datapoint population", "datapoints individuals"),
             col = c("black", "red", "grey", "black", "red", "grey"), 
             lty = c(1, 1, 1, NA, NA, NA), 
             pch = c(NA, NA, NA, 16, 16, 16), 
             lwd = c(4, 4, 1, NA, NA, NA), pt.cex = c(NA, NA, NA, 1.5, 1.5, 0.6))
    } else {
      legend("bottomright", legend = c("datapoint conditional", "datapoint population", "datapoints individuals"),
             col = c("black", "red", "grey"), 
             lty = c(NA, NA, NA), 
             pch = c(16, 16, 16), 
             lwd = c(NA, NA, NA), pt.cex = c(1.5, 1.5, 0.6))
    }
   
      
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
