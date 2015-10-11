

simul_copula <- function(smpl_size, params, model) {
  # Purpose: draws a sample from a model-specified copula
  # Inputs: integer scalar smpl_size, double vector params,
  #        string model
  # Output: numeric matrix cop_smpl
  # Notes: params order- alpha, nu, eta, lambda, gamma,
  #        chi, phi, 3 rhos, omega
  
  if (model == 'independent') {
    cop_smpl <- matrix(runif(3 * smpl_size), nrow = smpl_size)
    return(cop_smpl)
  }
  cop_pdf <- ellipCopula(family = model, 
                         param = params[c(8, 9, 10)],        
                         dim = 3, dispstr = 'un')
  cop_smpl <- rCopula(n = smpl_size, copula = cop_pdf)
  cop_smpl
}

calc_param <- function(cop_smpl, params) {
  # Purpose: transforms copula draws to parameter draws using 
  #          inverse cdfs and adds the threshold parameter to 
  #          completely define a random walk
  # Inputs: double matrix cop_smpl, double vector params, 
  # Output: double matrix param_smpl
  # Notes: params order- alpha, nu, eta, lambda, gamma,
  #        chi, phi, 3 rhos, omega
  
  param_smpl <- matrix(0, nrow = nrow(cop_smpl), ncol = 4,
                       dimnames = list(row = seq_len(nrow(cop_smpl)),
                                       col = c("alpha", "delta",
                                               "beta", "t_nd")))
  
  shape1 <- as.numeric(((1 - params["lambda"]) / params["gamma"] ^ 2 -
                         1 / params["lambda"]) * params["lambda"] ^ 2)
  shape2 <- as.numeric(shape1 * (1 / params["lambda"] - 1))
  shape <- (as.numeric(params[1, "chi"]) / as.numeric(params[1, "phi"])) ^ 2
  scale <- as.numeric(params[1, "phi"]) ^ 2 / as.numeric(params[1, "chi"])
  
  param_smpl[, "alpha"] <- as.numeric(params["alpha"])
  param_smpl[, "delta"] <- qnorm(p = cop_smpl[, 1], 
                                 mean = as.numeric(params["nu"]),
                                 sd = as.numeric(params["eta"]))
  param_smpl[, "beta"] <- qbeta(p = cop_smpl[, 2], 
                                shape1 = shape1, shape2 = shape2)
  param_smpl[, "t_nd"] <- qgamma(p = cop_smpl[, 3], 
                                 shape = shape,
                                 scale = scale)
  param_smpl
}

smpl_param <- function(params, smpl_size, model) {
  # Purpose: wraps simulation and transformation of copula draws
  #          to give a parameter sample
  # Inputs: double vector params, integer scalar smpl_size,
  #         string scalar model
  # Ouput: double data_frame param_smpl

  cop_smpl <- simul_copula(smpl_size = smpl_size, 
                           params = params,
                           model = model)
  param_smpl <- as.data.frame(calc_param(cop_smpl = cop_smpl, 
                                         params = params))
  param_smpl
}
