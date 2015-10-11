
setwd("~/Dropbox/Slava/Masters/")

library("cubature")
library("RWiener")
source("src/predictions/behavior/calc_copula_dens.R")

calc_prob_integrand <- function(dimensions, diffusion, params, model) {
  # Purpose: Calculates integrand for the probability of choice associated
  # with the lower bound
  # Input: numeric scalars drift, bias, nondec, diffusion, 
  # numeric vector params, character scalar model
  # Output: numeric scalar integrand
  
  drift <- dimensions[1]
  bias <- dimensions[2]
  nondec <- dimensions[3]
  
  drift_trans <- drift / (1 - drift ^ 2)
  nondec_trans <- nondec / (1 - nondec)
  
  choice_prob <- pwiener(q = .Machine$double.xmax, 
                         alpha = params[1, "alpha"] / diffusion, 
                         tau = nondec_trans, beta = bias, 
                         delta = drift_trans / diffusion,
                         resp = "lower")
  param_dens <- calc_density(drift = drift_trans, bias = bias,
                           nondec = nondec_trans,
                           params = params, model = model)
  jacobian <- (1 + drift ^ 2) / (1 - drift ^ 2) ^ 2 / (1 - nondec)
  integrand <- choice_prob * param_dens * jacobian
  return(integrand)
}

integrate_prob <- function(diffusion, params, model, tol, maxEval) {
  # Purpose: Integrates choice probability with respect to parameters' density
  # Input: numerical scalar diffusion, numerical vector params, 
  # character scalar model, numerical scalar tol
  # Output: numerical scalar integral_out
  
  integral_out <- adaptIntegrate(f = calc_prob_integrand, 
                                 lowerLimit = c(-1, 0, 0), 
                                 upperLimit = c(1, 1, 1), 
                                 diffusion, choice, params,
                                 model, tol = tol, fDim = 1,
                                 maxEval = maxEval)
  return(integral_out)
}


calc_quantile_integrand <- function(dimensions, prob, diffusion, 
                                    params, choice, model) {
  # Purpose: Calculates integrand for a rt quantile for a given a 
  # choice and model
  # Input: numeric vector dimensions, numeric scalars prob, diffusion, 
  # numeric vector params, character scalars choice, model
  # Output: numeric scalar integrand 
  
  drift <- dimensions[1]
  bias <- dimensions[2]
  nondec <- dimensions[3]
  
  drift_trans <- drift / (1 - drift ^ 2)
  nondec_trans <- nondec / (1 - nondec)
  
  norm_constant <- pwiener(q = .Machine$double.xmax,
                           alpha = params[1, "alpha"] / diffusion, 
                           tau = nondec_trans, beta = bias, 
                           delta = drift_trans / diffusion,
                           resp = choice)
  quant <- qwiener(p = prob * norm_constant, 
                   alpha = params[1, "alpha"] / diffusion, 
                   tau = nondec_trans, beta = bias, 
                   delta = drift_trans / diffusion,
                   resp = choice)
  if (is.nan(quant)) quant <- nondec_trans
  if (is.infinite(quant)) quant <- .Machine$double.xmax
  param_dens <- calc_density(drift = drift_trans, bias = bias,
                             nondec = nondec_trans,
                             params = params, model = model)
#  print(c(prob * norm_constant, params[1, "alpha"],
#          nondec_trans, bias, drift_trans / diffusion, quant, param_dens))
  jacobian <- (1 + drift ^ 2) / (1 - drift ^ 2) ^ 2 / (1 - nondec)
  integrand <- quant * param_dens * jacobian
  return(integrand)
}

integrate_quantile <- function(prob, diffusion, params, choice, 
                               model, tol, maxEval) {
  # Integrates choice probability with respect to parameters' density
  # Input: numerical scalar prob, diffusion, numerical vector params, 
  # character scalars choice, model, numerical scalar tol
  # Output: numerical scalar integral_out
  
  integral_out <- adaptIntegrate(f = calc_quantile_integrand,  
                                 lowerLimit = c(-1, 0, 0), 
                                 upperLimit = c(1, 1, 1), 
                                 prob, diffusion, params,
                                 choice, model, tol = tol, 
                                 fDim = 1, maxEval = maxEval)
  return(integral_out)
}

choice = "lower"

integrate_quantile(prob, diffusion, params, choice, 
                   model[2], tol = 1e-4, maxEval = 0) * 1000


root <- function(q, prob, alpha, tau, beta, delta, resp, eps = 1-3) {
  fn <- lower_CDF(q, a = alpha, v = delta, s = .1,
                  z = beta * alpha, t_nd = tau, eps) - prob
  return(fn)
}

root <- function(q, prob, alpha, tau, beta, delta, resp) {
  fn <- pwiener(q, alpha, tau, beta, delta, resp) - prob
  return(fn)
}

q <- 0.5140067
prob <- .9
alpha <- .223
tau <- .5
beta <- .5
delta <- -100
resp <- "upper"

uniroot(f = root, interval = c(0, 3), 
        prob, alpha, tau, beta, delta, resp, 
        extendInt = "upX", check.conv = TRUE, tol = 1e-4, maxiter = 10000)

qwiener(prob, alpha / diffusion, tau, beta, delta / diffusion, resp)

microbenchmark(uniroot(f = root, interval = c(0, 3), 
                       prob, alpha, tau, beta, delta, resp, 
                       extendInt = "upX", check.conv = TRUE, 
                       tol = 1e-4, maxiter = 1000),
               qwiener(prob, alpha, tau, beta, delta, resp),
               times = 100, unit = "eps")

pwiener(q, alpha / diffusion, tau, beta, delta / diffusion, resp)
lower_CDF(q, a = alpha, v = delta, s = .1,
          z = beta * alpha, t_nd = tau, eps)












