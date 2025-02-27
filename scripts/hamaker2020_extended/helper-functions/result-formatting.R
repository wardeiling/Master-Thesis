# RESULT FORMATTING FUNCTION
glmm_formating_results <- function(models) {
  
  # Define the new ordering of rows
  row_order <- c(
    "l1", "g.independence1", "g.exchangeable1", "g.ar11",
    "l2", "g.independence2", "g.exchangeable2", "g.ar12",
    "l3a", "g.independence3a", "g.exchangeable3a", "g.ar13a",
    "l4", "g.independence4", "g.exchangeable4", "g.ar14"
  )
  
  # Initialize the table with the new row order
  table <- matrix(NA, nrow = length(row_order), ncol = 3, 
                  dimnames = list(row_order, c("X", "X.cent", "X.cluster.means")))
  
  # Populate the table following the new row order
  for (i in seq_along(row_order)) {
    name <- row_order[i]
    if (name %in% names(models)) {
      if ("X" %in% names(models[[name]])) table[i, 1] <- models[[name]]["X"]
      if ("X.cent" %in% names(models[[name]])) table[i, 2] <- models[[name]]["X.cent"]
      if ("X.cluster.means" %in% names(models[[name]])) table[i, 3] <- models[[name]]["X.cluster.means"]
    }
  }
  
  return(table)
}