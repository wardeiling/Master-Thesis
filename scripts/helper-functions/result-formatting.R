# RESULT FORMATTING FUNCTION
glmm_formating_results <- function(models) {
  
  table <- matrix(NA, nrow = 16, ncol = 3, 
                  dimnames = list(c("l1", "l2", "l3a", "l4", 
                                    "g.independence1", "g.exchangeable1", "g.ar11", 
                                    "g.independence2", "g.exchangeable2", "g.ar12", 
                                    "g.independence3a", "g.exchangeable3a", "g.ar13a", 
                                    "g.independence4", "g.exchangeable4", "g.ar14"),
                                  c("X", "X.cent", "X.cluster.means")))
  
  for (i in 1:16) {
    names <- names(models[[i]])
    if ("X" %in% names) table[i, 1] <- models[[i]]["X"]
    if ("X.cent" %in% names) table[i, 2] <- models[[i]]["X.cent"]
    if ("X.cluster.means" %in% names) table[i, 3] <- models[[i]]["X.cluster.means"]
  }
  
  # Define the new ordering of rows
  row_order <- c(
    "l1", "g.independence1", "g.exchangeable1", "g.ar11",
    "l2", "g.independence2", "g.exchangeable2", "g.ar12",
    "l3a", "g.independence3a", "g.exchangeable3a", "g.ar13a",
    "l4", "g.independence4", "g.exchangeable4", "g.ar14"
  )
  
  # Reorder the table
  table <- table[row_order, ]
  
  return(table)
}