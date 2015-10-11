

reflect <- function(theta_new) {
  # Purpose: checks that a MCMC proposal falls into the posterior support
  #          and reflects it back if not
  # Input: double vector theta_new, character scalar model (global)
  # Output: double vector theta_new
  
  if (model == "independent")
    return(theta_new)
  
  results <- logical(3)
  
  while (any(!results)) {
      test <- findInterval(theta_new[21], c(-1, 1))
      if (test == 1)
        results[1] <- TRUE
      else if (test == 0)
        theta_new[21] <- -2 - theta_new[21]
      else if (test == 2)
        theta_new[21] <- 2 - theta_new[21]
      
      test <- findInterval(theta_new[22], c(-1, 1))
      if (test == 1)
        results[2] <- TRUE
      else if (test == 0)
        theta_new[22] <- -2 - theta_new[22]
      else if (test == 2)
        theta_new[22] <- 2 - theta_new[22]
      
      test <- findInterval(theta_new[23], c(-1, 1))
      if (test == 1)
        results[3] <- TRUE
      else if (test == 0)
        theta_new[23] <- -2 - theta_new[23]
      else if (test == 2)
        theta_new[23] <- 2 - theta_new[23]
  }
  theta_new
}
