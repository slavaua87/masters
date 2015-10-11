

### Generates trial parameters from the joint density of 
### drift, bias, nondecision
simul_copula <- function(smpl_size, params, model) {
  # Purpose: draws a sample from a model-specified copula
  # Input: numeric scalar smpl_size, numeric vector params,
  # character scalar model
  # Output: numeric matrix cop_smpl
  # Notes: params order- alpha, nu, eta, lambda, gamma, chi, phi, 3 rhos, omega
  
  if (model == 'independent') {
    cop_smpl <- matrix(runif(n = 3 * smpl_size), nrow = smpl_size)
    return(cop_smpl)
  }
  cop_pdf <- ellipCopula(family = model, 
                         param = 
                           as.numeric(params[c("rho_db", "rho_dt", "rho_bt")]), 
                         dim = 3, dispstr = 'un', 
                         df = as.numeric(params["omega"]))
  cop_smpl <- rCopula(n = smpl_size, copula = cop_pdf)
  return(cop_smpl)
}

#ptnorm_quant <- function(p, interval, mean, sigma, lower, upper, ...) {
  # Purpose: calculates the quantile function of a truncated normal variable
  # Input: numeric scalar p, numeric vector interval, numeric scalars mean,
  # sigma, lower, upper, optional arguments to uniroot ...
  # Output: numeric scalar quant

#  ptnorm_root <- function(q) {
#    lower_prob <- pnorm(q = lower, mean = mean,
#                        sd = sigma, lower.tail = TRUE)
#    upper_prob <- pnorm(q = upper, mean = mean,
#                        sd = sigma, lower.tail = TRUE)
#    root <- (pnorm(q = q, mean = mean, sd = sigma, lower.tail = TRUE) - 
#               lower_prob) / (upper_prob - lower_prob) - p
#    return(root)
#  }
#  quant <- uniroot(f = ptnorm_root, interval = interval, ...)$root
#  return(quant)
#}

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
  
#  param_smpl[, "t_nd"] <- sapply(X = cop_smpl[, 3],
#                                 FUN = ptnorm_quant,
#                                 interval = c(.Machine$double.xmin, 10), 
#                                 mean = as.numeric(params[1, "chi"]),
#                                 sigma = as.numeric(params[1, "phi"]),
#                                 lower = 0, upper = Inf,
#                                 extendInt = "upX", tol = 1e-4)
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
