

simul_copula <- function(smpl_size, params, model) {
  # Purpose: draws a sample from a model-specified copula
  # Input: numeric scalar smpl_size, numeric vector params,
  # character scalar model
  # Output: numeric matrix cop_smpl
  # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
  # 4:lambda, 5:gamma, 6:chi, 7:phi, 8:rho_db, 9:rho_dt, 10:rho_bt
  
  if (model == 'independent') {
    cop_smpl <- matrix(runif(n = 3 * smpl_size), nrow = smpl_size)
    return(cop_smpl)
  }
  cop_pdf <- ellipCopula(family = model, 
                         param = as.numeric(
                           params[c(8, 9, 10)]),        
                         dim = 3, dispstr = 'un')
  cop_smpl <- rCopula(n = smpl_size, copula = cop_pdf)
  return(cop_smpl)
}

calc_param <- function(cop_smpl, params) {
  # Purpose: transforms copula draws to parameter draws using inverse cdfs
  # and adds the threshold parameter to completely define a random walk
  # Input: numeric matrix cop_smpl, numeric vector params, 
  # Output: numeric matrix param_smpl
  # Notes: params order- alpha, nu, eta, lambda, gamma, chi, phi, 3 rhos, omega
  
  param_smpl <- matrix(0, nrow = nrow(cop_smpl), ncol = 4,
                       dimnames = list(row = seq_len(nrow(cop_smpl)),
                                       col = c("alpha", "delta",
                                               "beta", "t_nd")))
  
  shape1 <- as.numeric(((1 - params[4]) / params[5] ^ 2 -
                         1 / params[4]) * params[4] ^ 2)
  shape2 <- as.numeric(shape1 * (1 / params[4] - 1))
  shape <- as.numeric((params[6] / params[7]) ^ 2)
  scale <- as.numeric(params[7] ^ 2 / params[6])
  param_smpl[, "alpha"] <- as.numeric(params[1])
  param_smpl[, "delta"] <- qnorm(p = cop_smpl[, 1], 
                                 mean = as.numeric(params[2]),
                                 sd = as.numeric(params[3]))
  param_smpl[, "beta"] <- qbeta(p = cop_smpl[, 2], 
                                shape1 = shape1, shape2 = shape2)
  param_smpl[, "t_nd"] <- qgamma(cop_smpl[, 3], shape = shape, scale = scale)
  return(param_smpl)
}

smpl_param <- function(params, smpl_size, model) {
  # Purpose: wraps simulation and transformation of copula draws
  # to give a parameter sample
  # Input: numeric vector params, numeric scalar smpl_size,
  # character scalar model
  # Ouput: numeric data.frame param_smpl
  cop_smpl <- simul_copula(smpl_size = smpl_size, 
                           params = params,
                           model = model)
  param_smpl <- as.data.frame(calc_param(cop_smpl = cop_smpl, 
                                         params = params))
  return(param_smpl)
}
