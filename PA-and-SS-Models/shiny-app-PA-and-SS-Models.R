# Load necessary libraries
library(shiny)

# Define UI for application
ui <- fluidPage(
  titlePanel("Mixed Models: GLMM and MLM"),
  
  tabsetPanel(
    # Tab for Generalized Mixed Linear Model
    tabPanel("GLMM with Log-Link",
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
                 h3("GLMM: Population-Averaged vs Subject-Specific Effects"),
                 plotOutput("glmmPlot"),
                 p("This plot illustrates the difference between the marginal (population-averaged) and conditional (subject-specific) relationships in a generalized mixed linear model (GLMM)
                   with a log-link function. 
                   
                   In the subject-specific (SS) model the fixed slope Beta_1 represents on average how an individual's probability of a positive response depends on X (i.e., effect for average individual with b_i = 0). 
                   
                   In the population-averaged (PA) model, the fixed slope Beta_1 represents the change on a logit scale in the fraction of positive responses for a one-unit change in X. 
                   
                   We can see that they are identical when there is no heterogeneity (i.e., random effects variances set to 0). Use the sliders on the left to adjust the parameters and observe their effects. ")
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
  
  # Generate plot for GLMM
  output$glmmPlot <- renderPlot({
    # Parameters from sliders
    alpha_0 <- input$alpha_0
    alpha_1 <- input$alpha_1
    sigma_b0 <- input$sigma_b0
    
    set.seed(123) # Ensure reproducibility
    n_individuals <- 20 # Number of individuals to plot
    X_grid <- seq(-5, 5, length.out = 10000) # Values for X_it
    
    # Generate random intercepts
    b0 <- rnorm(n_individuals, mean = 0, sd = sigma_b0)
    
    # Conditional-level logit and probability
    logit_conditional <- alpha_0 + alpha_1 * X_grid + 0
    pi_conditional <- 1 / (1 + exp(-logit_conditional))
    
    # Individual-specific logits and probabilities
    logit_individuals <- sapply(b0, function(b) alpha_0 + alpha_1 * X_grid + b)
    pi_individuals <- 1 / (1 + exp(-logit_individuals))
    
    # Population-level probability
    pi_population <- rowMeans(pi_individuals)
    
    # Plot
    plot(X_grid, pi_conditional, type = "l", lwd = 4, col = "black",
         ylim = c(0, 1), xlim = c(-5, 5), xlab = expression(X[it]), ylab = expression(pi),
         main = "Population-averaged, Conditional and Individual Logit Curves")
    
    # Add individual curves
    for (i in 1:n_individuals) {
      lines(X_grid, pi_individuals[, i], col = "grey", lwd = 1)
    }
    
    # Add population-level curve
    lines(X_grid, pi_population, col = "red", lwd = 4)
    
    # Add legend
    legend("topleft", legend = c("conditional", "population", "individuals"),
           col = c("black", "red", "grey"), lty = c(1, 1, 1), lwd = c(4, 4, 1))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
