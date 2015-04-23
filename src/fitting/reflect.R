

reflect <- function(theta_new, model) {
  # Purpose: checks that a DEMCMC proposal falls into the posterior support
  # Input: numerical vector theta_new, character scalar model
  # Output: logical scalar
  
  results <- logical(length(theta_new))
  
  while (any(!results)) {
    if (model == "normal") {
      test <- findInterval(theta_new[18], c(-1, 1))
      if (test == 1)
        results[18] <- TRUE
      else if (test == 0)
        theta_new[18] <- -2 - theta_new[18]
      else if (test == 2)
        theta_new[18] <- 2 - theta_new[18]
      
      test <- findInterval(theta_new[19], c(-1, 1))
      if (test == 1)
        results[19] <- TRUE
      else if (test == 0)
        theta_new[19] <- -2 - theta_new[19]
      else if (test == 2)
        theta_new[19] <- 2 - theta_new[19]
      
      test <- findInterval(theta_new[20], c(-1, 1))
      if (test == 1)
        results[20] <- TRUE
      else if (test == 0)
        theta_new[20] <- -2 - theta_new[20]
      else if (test == 2)
        theta_new[20] <- 2 - theta_new[20]
    }
    
    test <- findInterval(theta_new[1], c(0, .786))
    if (test == 1)
      results[1] <- TRUE
    else if (test == 0)
      theta_new[1] <-  -theta_new[1]
    else if (test == 2)
      theta_new[1] <- 2 * .786 - theta_new[1]
    
    test <- findInterval(theta_new[2], c(0, 1.172))
    if (test == 1)
      results[2] <- TRUE
    else if (test == 0)
      theta_new[2] <-  -theta_new[2]
    else if (test == 2)
      theta_new[2] <- 2 * 1.172 - theta_new[2]

    test <- findInterval(theta_new[3], c(0, 1.172))
    if (test == 1)
      results[3] <- TRUE
    else if (test == 0)
      theta_new[3] <-  -theta_new[3]
    else if (test == 2)
      theta_new[3] <- 2 * 1.172 - theta_new[3]
    
    test <- findInterval(theta_new[4], c(0, 50))
    if (test == 1)
      results[4] <- TRUE
    else if (test == 0)
      theta_new[4] <-  -theta_new[4]
    else if (test == 2)
      theta_new[4] <- 2 * 50 - theta_new[4]
    
    test <- findInterval(theta_new[5], c(0, 5))
    if (test == 1)
      results[5] <- TRUE
    else if (test == 0)
      theta_new[5] <-  -theta_new[5]
    else if (test == 2)
      theta_new[5] <- 2 * 5 - theta_new[5]
    
    test <- findInterval(theta_new[6], c(0, 1))
    if (test == 1)
      results[6] <- TRUE
    else if (test == 0)
      theta_new[6] <- -theta_new[6]
    else if (test == 2)
      theta_new[6] <- 2 - theta_new[6]
    
    test <- findInterval(theta_new[7], c(0, 1.884))
    if (test == 1)
      results[7] <- TRUE
    else if (test == 0)
      theta_new[7] <- -theta_new[7]
    else if (test == 2)
      theta_new[7] <- 2 * 1.884 - theta_new[7]
    
    test <- findInterval(theta_new[8], c(0, .786))
    if (test == 1)
      results[8] <- TRUE
    else if (test == 0)
      theta_new[8] <-  -theta_new[8]
    else if (test == 2)
      theta_new[8] <- 2 * .786 - theta_new[8]
    
    test <- findInterval(theta_new[9], c(0, 1.172))
    if (test == 1)
      results[9] <- TRUE
    else if (test == 0)
      theta_new[9] <-  -theta_new[9]
    else if (test == 2)
      theta_new[9] <- 2 * 1.172 - theta_new[9]
    
    test <- findInterval(theta_new[10], c(0, 1.172))
    if (test == 1)
      results[10] <- TRUE
    else if (test == 0)
      theta_new[10] <-  -theta_new[10]
    else if (test == 2)
      theta_new[10] <- 2 * 1.172 - theta_new[10]
    
    test <- findInterval(theta_new[11], c(0, 50))
    if (test == 1)
      results[11] <- TRUE
    else if (test == 0)
      theta_new[11] <-  -theta_new[11]
    else if (test == 2)
      theta_new[11] <- 2 * 50 - theta_new[11]
    
    test <- findInterval(theta_new[12], c(0, 5))
    if (test == 1)
      results[12] <- TRUE
    else if (test == 0)
      theta_new[12] <-  -theta_new[12]
    else if (test == 2)
      theta_new[12] <- 2 * 5 - theta_new[12]
    
    test <- findInterval(theta_new[13], c(0, 1))
    if (test == 1)
      results[13] <- TRUE
    else if (test == 0)
      theta_new[13] <- -theta_new[13]
    else if (test == 2)
      theta_new[13] <- 2 - theta_new[13]
    
    test <- findInterval(theta_new[14], c(0, 1.884))
    if (test == 1)
      results[14] <- TRUE
    else if (test == 0)
      theta_new[14] <- -theta_new[14]
    else if (test == 2)
      theta_new[14] <- 2 * 1.884 - theta_new[14]
    
    test <- findInterval(theta_new[15], c(0, .658))
    if (test == 1)
      results[15] <- TRUE
    else if (test == 0)
      theta_new[15] <- -theta_new[15]
    else if (test == 2)
      theta_new[15] <- 2 * .658 - theta_new[15]
    
    test <- findInterval(theta_new[16], c(0, .5))
    if (test == 1)
      results[16] <- TRUE
    else if (test == 0)
      theta_new[16] <- -theta_new[16]
    else if (test == 2)
      theta_new[16] <- 2 * .5 - theta_new[16]
    
    test <- findInterval(theta_new[17], c(0, 1.260))
    if (test == 1)
      results[17] <- TRUE
    else if (test == 0)
      theta_new[17] <- -theta_new[17]
    else if (test == 2)
      theta_new[17] <- 2 * 1.260 - theta_new[17]
  }
  return(theta_new)
}
