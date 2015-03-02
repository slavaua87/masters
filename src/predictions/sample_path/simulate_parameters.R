

### Generates trial parameters from the joint density of drift and bias
library(copula)

simul_copula <- function(smpl_size, params, model) {
  # Draws a sample from a specified copula
  # Takes sample size, parameter values and model type to return a matrix
  # Params order: alpha, nu, eta, lambda, gamma, chi, phi, 3 rhos, omega
  
  if (model == 'independent') {
    cop_smpl <- matrix(runif(n = 2 * smpl_size), nrow = smpl_size, ncol = 2)
    return(cop_smpl)
  }
  
  cop_pdf <- ellipCopula(family = model, param = params["rho_db"], 
                         dim = 2, dispstr = 'un', df = params["omega"])
  cop_smpl <- rCopula(n = smpl_size, copula = cop_pdf)
  return(cop_smpl)
}

calc_param <- function(cop_smpl, params) {
  # Transforms copula draws to parameter draws using inverse cdfs
  # Takes a draws matrix, parameter vector, truncation points to output a matrix
  # Params order: alpha, nu, eta, lambda, gamma, chi, phi, 3 rhos, omega
  
  param_smpl <- matrix(0, nrow = dim(cop_smpl)[1], ncol = 3,
                       dimnames = list(row = seq_len(dim(cop_smpl)[1]),
                                       col = c("alpha", "delta", "beta")))
  
  param_smpl[, 1] <- as.numeric(params["alpha"])
  param_smpl[, 2] <- qnorm(cop_smpl[, 1], as.numeric(params["nu"]),
                           as.numeric(params["eta"]))
  shape1 <- as.numeric(((1 - params["lambda"]) / params["gamma"] ^ 2 -
              1 / params["lambda"]) * params["lambda"] ^ 2)
  shape2 <- as.numeric(shape1 * (1 / params["lambda"] - 1))
  param_smpl[, 3] <- qbeta(cop_smpl[, 2], shape1, shape2)
  return(param_smpl)
}

smpl_param <- function(params, smpl_size, model) {
  # Wraps simulation and transformation of copula draws to give parameter sample
  # Takes a scalar, vector, string, scalar to output a data frame
  
  cop_smpl <- simul_copula(smpl_size = smpl_size, params = params, model = model)
  param_smpl <- calc_param(cop_smpl = cop_smpl, params = params)
  return(as.data.frame(param_smpl))
}




