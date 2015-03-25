
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
  # Purpose: checks that a DEMCMC proposal falls into the posterior support
  # Input: numerical vector theta_new, character scalar model
  # Output: logical scalar
  
  if (!all(check_beta(theta_new[c(6, 13)], theta_new[16])))
    return(FALSE)
  if (!(theta_new[1] >= 0 && theta_new[1] <= .786))
    return(FALSE)
  if (!(theta_new[2] >= 0 && theta_new[2] <= 1.172))
    return(FALSE)
  if (!(theta_new[3] >= 0 && theta_new[3] <= 1.172))
    return(FALSE)
  if (!(theta_new[4] >= 0 && theta_new[4] <= 50))
    return(FALSE)
  if (!(theta_new[5] >= 0 && theta_new[5] <= 5))
    return(FALSE)
  if (!(theta_new[6] >= 0 && theta_new[6] <= 1))
    return(FALSE)
  if (!(theta_new[7] >= 0 && theta_new[7] <= 1.884))
    return(FALSE)
  if (!(theta_new[8] >= 0 && theta_new[8] <= .786))
    return(FALSE)
  if (!(theta_new[9] >= 0 && theta_new[9] <= 1.172))
    return(FALSE)
  if (!(theta_new[10] >= 0 && theta_new[10] <= 1.172))
    return(FALSE)
  if (!(theta_new[11] >= 0 && theta_new[11] <= 50))
    return(FALSE)
  if (!(theta_new[12] >= 0 && theta_new[12] <= 5))
    return(FALSE)
  if (!(theta_new[13] >= 0 && theta_new[13] <= 1))
    return(FALSE)
  if (!(theta_new[14] >= 0 && theta_new[14] <= 1.884))
    return(FALSE)
  if (!(theta_new[15] >= 0 && theta_new[15] <= .658))
    return(FALSE)
  if (!(theta_new[16] >= 0 && theta_new[16] <= .5))
    return(FALSE)
  if (!(theta_new[17] >= 0 && theta_new[17] <= 1.260))
    return(FALSE)
  if (model == "normal") {
    if (!check_rho(theta_new[18:20]))
      return(FALSE)
  }
  return(TRUE)
}
