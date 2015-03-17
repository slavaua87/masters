
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
  
  if (model == "independent") {
    if (!theta_new[1:2] >= 0 && theta_new[1:2] <= .786)
      return(FALSE)
    if (!theta_new[3:6] >= 0 && theta_new[3:6] <= 1.172)
      return(FALSE)
    if (!theta_new[7:8] >= 0 && theta_new[7:8] <= 50)
      return(FALSE)
    if (!theta_new[9:10] >= 0 && theta_new[9:10] <= 5)
      return(FALSE)
    if (!theta_new[11:12] >= 0 && theta_new[11:12] <= 1)
      return(FALSE)
    if (!theta_new[13:14] >= 0 && theta_new[13:14] <= 1.884)
      return(FALSE)
    if (!all(theta_new[15:19] >= 0))
      return(FALSE)
    if (!theta_new[20] >= 0 && theta_new[20] <= .5)
      return(FALSE)
    if (!theta_new[21] >= 0)
      return(FALSE)
    if (!all(theta_new[22:51] > 0))
        return(FALSE)
    if (!theta_new[52:57] >= 0 && theta_new[52:57] <= 1)
      return(FALSE)
    if (!all(theta_new[58:63] > 0))
      return(FALSE)
    if (!theta_new[64:66] >= 0 && theta_new[64:66] <= .658)
      return(FALSE)
    if (!theta_new[67:69] >= 0 && theta_new[67:69] <= .860)
      return(FALSE)
    if (!theta_new[70:72] >= 0 && theta_new[70:72] <= 1.260)
      return(FALSE)
  }
  if (model == "normal") {
    if (!theta_new[1:2] >= 0 && theta_new[1:2] <= .786)
      return(FALSE)
    if (!theta_new[3:6] >= 0 && theta_new[3:6] <= 1.172)
      return(FALSE)
    if (!theta_new[7:8] >= 0 && theta_new[7:8] <= 50)
      return(FALSE)
    if (!theta_new[9:10] >= 0 && theta_new[9:10] <= 5)
      return(FALSE)
    if (!theta_new[11:12] >= 0 && theta_new[11:12] <= 1)
      return(FALSE)
    if (!theta_new[13:14] >= 0 && theta_new[13:14] <= 1.884)
      return(FALSE)
    if (!all(theta_new[15:19] >= 0))
      return(FALSE)
    if (!theta_new[20] >= 0 && theta_new[20] <= .5)
      return(FALSE)
    if (!all(theta_new[21:22] >= 0))
      return(FALSE)
    if (!all(theta_new[23:52] > 0))
      return(FALSE)
    if (!theta_new[53:58] >= 0 && theta_new[53:58] <= 1)
      return(FALSE)
    if (!all(theta_new[59:64] > 0))
      return(FALSE)
    if (!theta_new[65:67] >= 0 && theta_new[65:67] <= .658)
      return(FALSE)
    if (!theta_new[68:70] >= 0 && theta_new[68:70] <= .860)
      return(FALSE)
    if (!theta_new[71:73] >= 0 && theta_new[71:73] <= 1.260)
      return(FALSE)
    if (!check_rho(theta_new[74:76]))
      return(FALSE)
    if (!check_rho(theta_new[77:79]))
      return(FALSE)
    if (!check_rho(theta_new[80:82]))
      return(FALSE)
  }
  return(TRUE)
}
