
check_support <- function(theta_new) {
  # Purpose: check proposed parameters are within proper sets
  # Input: numeric vector theta_new
  # Output: numeric vector
  
  drift <- theta_new[c(2, 12)] < theta_new[c(3, 13)]
    
  if (model == "normal") {
    rho <- theta_new[21:23]
    det_rho <- 1 - sum(rho ^ 2) + 2 * prod(rho)
    test <- if (det_rho > 0) TRUE else FALSE
    drift <- c(drift, test)
  } 
  all(drift)
}