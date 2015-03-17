
setwd("~/Dropbox/Slava/Masters/")

library("cubature")
library("RWiener")
source("src/predictions/behavior/calc_copula_dens.R")

calc_density_integrand <- function(dimensions, rt, choice,
                                   sigma, params, 
                                   shape1, shape2, trunc_const, sqrt_rho,
                                   lower, upper, 
                                   model) {
  # Purpose: Calculates integrand for a rt density for a given a 
  # choice and model
  # Input: numeric vector dimensions, numeric scalars rt, choice, sigma
  # numeric vector params, numeric scalars lower, upper, character scalars model
  # Output: numeric scalar integrand
  # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
  # 4:lambda, 5:gamma, 6:chi, 8:phi, 9:rho_db, 10:rho_dt, 11:rho_bt, 12:omega
  
  if (choice == 1)
    choice = "upper"
  if (choice == 2)
    choice = "lower"
  
  delta <- dimensions[1]
  beta <- dimensions[2]
  t_nd <- dimensions[3]
  
  delta_trans <- delta / (1 - delta ^ 2)
  
  behav_dens <- dwiener(q = rt,
                        alpha = params[1] / sigma, 
                        tau = t_nd, 
                        beta = beta, 
                        delta = delta_trans / sigma,
                        resp = choice)
  param_dens <- calc_density(delta = delta_trans, beta = beta,
                             t_nd = t_nd,
                             params = params,
                             shape1 = shape1, shape2 = shape2,
                             trunc_const = trunc_const, sqrt_rho = sqrt_rho,
                             lower = lower, upper = upper, model = model)
  jacobian <- (1 + delta ^ 2) / (1 - delta ^ 2) ^ 2 
  integrand <- behav_dens * param_dens * jacobian
  return(integrand)
}

integrate_density <- function(rt, choice, sigma, params, lower,
                              upper, model, 
                              tol = 1e-2, maxEval = 0) {
  # Integrates choice probability with respect to parameters' density
  # Input: numerical scalars rt, choice, sigma, numerical vector params,
  # numeric scalars lower, upper, character scalar model, 
  # logical scalar give_log, numerical scalar tol
  # Output: numerical scalar behav_mixture
  # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
  # 4:lambda, 5:gamma, 6:chi, 7:phi, 8:rho_db, 9:rho_dt, 10:rho_bt, 11:omega
  
  shape1 <- ((1 - params[4]) / params[5] ^ 2 -
               1 / params[4]) * params[4] ^ 2
  shape2 <- shape1 * (1 / params[4] - 1)
  
  trunc_const <- pnorm(q = c(lower, upper), 
                       mean = params[6],
                       sd = params[7],
                       lower.tail = TRUE)
  
  rho <- matrix(c(1, params[8], params[9],
                  params[8], 1, params[10],
                  params[9], params[10], 1), byrow = TRUE, ncol = 3)
  sqrt_rho <- cholesky(x = rho)
  
  behav_mixture <- adaptIntegrate(f = calc_density_integrand,  
                                  lowerLimit = c(-1, 0, 0), 
                                  upperLimit = c(1, 1, rt), 
                                  rt, choice, sigma, params,
                                  shape1, shape2, trunc_const, sqrt_rho,
                                  lower, upper,
                                  model, tol = tol, 
                                  fDim = 1, maxEval = maxEval)$integral
  return(behav_mixture)
}

integrate_density_vec <- function(rt, choice, sigma, params, lower,
                                  upper, model, tol, maxEval) {
  behav_mixture <- mapply(FUN = integrate_density_vec, rt, choice,
                          MoreArgs = list(sigma, params, lower,
                                          upper, model),
                          SIMPLIFY = FALSE,
                          USE.NAMES = FALSE)
 return(behav_mixture)
}

integrate_density_vec <- Vectorize(FUN = integrate_density, 
                                   vectorize.args = c("rt", "choice"))







