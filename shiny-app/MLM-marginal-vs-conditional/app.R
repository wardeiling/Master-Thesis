library(shiny)

ui <- fluidPage(
  titlePanel("Effect of X_t on Y_t+1"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("var_b", "Variance of b:", min = 0, max = 15, value = 3, step = 0.1),
      sliderInput("var_X1", "Variance of X1:", min = 0, max = 15, value = 1, step = 0.1),
      sliderInput("var_epsilon", "Variance of epsilon:", min = 0, max = 15, value = 1, step = 0.1),
      sliderInput("beta_0", "Beta_0:", min = -5, max = 5, value = 0.8, step = 0.1),
      sliderInput("beta_1", "Beta_1:", min = -5, max = 5, value = 2, step = 0.1)
    ),
    mainPanel(
      plotOutput("effectPlot", height = "500px", width = "600px"),
      tableOutput("effectTable"),
      h2("Observations"),
      p("If var(b) >> var(X1): the marginal slope of X2 on Y3 is greater than the conditional slope"),
      p("If 0.1 < var(b) << var(X1): the marginal slope of X2 on Y3 is smaller than the conditional slope"),
      p("If var(b) << var(X1): the marginal slope of X2 on Y3 is smaller than the conditional slope"),
      p("If 10 < var(b) = var(X1): the marginal slope of X2 on Y3 is slightly smaller than the conditional slope"),
      p("As rho and zeta approximate zero, the marginal slope of X2 on Y3 approaches the conditional slope."),
      p("If 2*var(b) = var(X1): the marginal slope of X2 on Y3 is close to the conditional slope."),
      h2("Appendix"),
      withMathJax(
        p("Assume the variables are generated according to the following multilevel linear model (MLM) with a random intercept:"),
        p("$$b_i \\sim N(0, \\sigma_b^2)$$"),
        p("$$X_{i1} \\sim N(0, \\sigma_{X_1}^2) \\text{ independent of } b_i$$"),
        p("$$Y_{i2} \\mid X_{i1}, b_i \\sim N(\\beta_0 + \\beta_1 X_{i1} + b_i, \\sigma_\\epsilon^2)$$"),
        p("$$X_{i2} = Y_{i2}$$"),
        p("$$Y_{i3} \\mid X_{i1}, Y_{i2}, X_{i2}, b_i \\sim N(\\beta_0 + \\beta_1 X_{i2} + b_i, \\sigma_\\epsilon^2)$$"),
        p("This implies the conditional relationship:"),
        p("$$E[Y_{i,t+1} \\mid X_{it}] = \\beta_0 + \\beta_1 X_{it} + b_i$$"),
        p("While the marginal relationship is more complex:"),
        p("$$E[Y_{i2} \\mid X_{i1}] = \\beta_0 + \\beta_1 X_{i1}$$"),
        p("$$E[Y_{i3} \\mid X_{i2}] = (1 - \\rho \\zeta - \\rho) \\beta_0 + [(1 - \\rho \\zeta) \\beta_1 + \\rho] X_{i2}$$"),
        p("where $$\\rho = \\frac{\\sigma_b^2}{\\sigma_b^2 + \\sigma_\\epsilon^2}$$ and $$\\zeta = \\frac{\\beta_1 \\sigma_{X_1}^2}{\\beta_1 \\sigma_{X_1}^2 + \\sigma_b^2 + \\sigma_\\epsilon^2}.$$")
      )
    )
  )
)

server <- function(input, output) {
  output$effectPlot <- renderPlot({
    # Extract input values for readability
    var_b <- input$var_b
    var_X1 <- input$var_X1
    var_epsilon <- input$var_epsilon
    beta_0 <- input$beta_0
    beta_1 <- input$beta_1
    
    # Compute rho and zeta
    rho <- var_b / (var_b + var_epsilon)
    zeta <- beta_1 * var_X1 / (beta_1 * var_X1 + var_b + var_epsilon)
    
    # Compute marginal effects
    marginal_intercept_x1_y2 <- beta_0
    marginal_slope_x1_y2 <- beta_1
    
    marginal_intercept_x2_y3 <- (1 - rho * zeta - rho) * beta_0
    marginal_slope_x2_y3 <- (1 - rho * zeta) * beta_1 + rho
    
    # Plot
    plot(1, type = "n", xlim = c(-15, 15), ylim = c(-40, 40), 
         xlab = expression(X[it]), ylab = expression(Y[it+1]), 
         main = "Effect of X_t on Y_t+1")
    
    abline(a = beta_0, b = beta_1, col = "blue", lty = 1, lwd = 4)
    abline(a = marginal_intercept_x1_y2, b = marginal_slope_x1_y2, col = "red", lwd = 2)
    abline(a = marginal_intercept_x2_y3, b = marginal_slope_x2_y3, col = "orange", lwd = 2)
    
    legend("topleft", legend = c("true conditional effect", "true marginal effect X1 Y2", "true marginal effect X2 Y3"),
           col = c("blue", "red", "orange"), lty = c(1, 1, 1), lwd = c(4, 2, 2))
  })
  
  output$effectTable <- renderTable({
    data.frame(
      Effect = c("Conditional", "Marginal X1 -> Y2", "Marginal X2 -> Y3"),
      Intercept = c(beta_0, marginal_intercept_x1_y2, marginal_intercept_x2_y3),
      Slope = c(beta_1, marginal_slope_x1_y2, marginal_slope_x2_y3)
    )
  })
}

shinyApp(ui = ui, server = server)
