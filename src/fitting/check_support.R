
check_beta <- function(lambda, gamma) {
  # Purpose: checks that mean and variance of a beta variable are proper
  # Input: numerical scalars lambda, gamma
  # Output: logical scalar
  
  restriction <- (gamma ^ 2 < lambda * (1 - lambda))
  return(restriction)
}

check_rho <- function(rho) {
  # Purpose: checks that correlations fall into positive definite set
  # Input: numerical vector rho
  # Ouput: logical scalar
  
  det_rho <- 1 - sum(rho ^ 2) + 2 * prod(rho)
  if (det_rho > 0)
    return(TRUE)
  else 
    return(FALSE)
}

check_support <- function(theta_new, model) {
  # Purpose: check proposed parameters are within proper sets
  # Input: numeric vector theta_new
  # Output: numeric vector
  
  if (!all(check_beta(theta_new[c(6, 13)], theta_new[16])))
    return(FALSE)
  if (model == "normal")
    if (!check_rho(theta_new[18:20]))
      return(FALSE)
  else
    return(TRUE)
}