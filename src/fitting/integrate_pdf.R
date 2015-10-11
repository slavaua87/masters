
integrator <- function(f, lowerLimit, upperLimit, ..., tol = 1e-05, fDim = 1, 
                       maxEval = 0, absError = 0, doChecking = FALSE) {
  fn <- function(x) f(x, ...)
  result <- .Call("doCubature", as.integer(fDim), body(fn), 
                  as.double(lowerLimit), as.double(upperLimit),
                  as.integer(maxEval), as.double(absError),
                  as.double(tol), environment(), PACKAGE = "cubature")
  names(result) <- c("integral", "error", "functionEvaluations", 
                     "returnCode")
  return(result)
}

calc_density_integrand <- function(dimensions, rt, choice,
                                   sigma, alpha, nu, eta, 
                                   shape1, shape2, 
                                   shape, scale,
                                   sqrt_rho, model) {
  # Purpose: Calculates integrand for a rt density for a given a 
  # choice and model
  # Input: numeric vector dimensions, numeric scalars rt, choice, sigma
  # numeric vector params, character scalar model
  # Output: numeric scalar integrand
  # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
  # 4:lambda, 5:gamma, 6:chi, 8:phi, 9:rho_db, 10:rho_dt, 11:rho_bt, 12:omega
  
  delta_trans <- dimensions[1]
  beta <- dimensions[2]
  t_nd <- dimensions[3]
  
  delta <- delta_trans / (1 - delta_trans ^ 2)
  
  behav_dens <- dwiener(q = rt,
                        alpha = alpha / sigma, 
                        tau = t_nd, 
                        beta = beta, 
                        delta = delta / sigma,
                        resp = choice)
  
  if (behav_dens <= 0 || is.infinite(behav_dens)) 
    behav_dens <- .Machine$double.xmin
  
  param_logdens <- calc_logdensity(delta = delta, beta = beta, t_nd = t_nd,
                                   nu = nu, eta = eta,
                                   shape1 = shape1, shape2 = shape2,
                                   shape = shape, scale = scale,
                                   sqrt_rho = sqrt_rho,
                                   model = model)
  jacobian <- (1 + delta_trans ^ 2) / (1 - delta_trans ^ 2) ^ 2 
  integrand <- exp(sum(log(behav_dens), param_logdens, log(jacobian)))
  if (is.infinite(integrand))
    integrand <- .Machine$double.xmax
  return(integrand)
}

integrate_density <- function(rt, choice, sigma, alpha, nu, eta,
                              shape1, shape2, shape, scale, rho_db,
                              rho_dt, rho_bt, model, 
                              tol = 1e-2, maxEval = 0) {
  # Integrates choice probability with respect to parameters' density
  # Input: numerical scalars rt, choice, sigma, alpha, nu, eta,
  # shape1, shape2, shape, scale, rho_db, rho_dt, rho_bt, character scalar model, 
  # numerical scalar tol, integer scalar maxEval
  # Output: numerical scalar behav_mixture
  
#   if (choice == 1)
#     choice = "upper"
#   if (choice == 2)
#     choice = "lower"
  
  if (model == "normal") {
    rho <- matrix(c(1, rho_db, rho_dt,
                    rho_db, 1, rho_bt,
                    rho_dt, rho_bt, 1),
                  byrow = TRUE, ncol = 3)
    sqrt_rho <- cholesky(x = rho)
  }
  eps_neg <- .Machine$double.neg.eps
  eps_pos <- .Machine$double.eps
  eps_min <- .Machine$double.xmin
  lower_lim <- c(-1 + eps_pos, eps_min, eps_min)
  upper_lim <- c(1 - eps_neg, 1 - eps_neg, rt)
  
  behav_mixture <- integrator(f = calc_density_integrand_c,  
                              lowerLimit = lower_lim, 
                              upperLimit = upper_lim, 
                              rt, choice, sigma, 
                              alpha, nu, eta,
                              shape1, shape2, shape, scale, 
                              sqrt_rho, model, tol = tol, 
                              fDim = 1, maxEval = maxEval)
  error_test <- behav_mixture$error / behav_mixture$integral
  if (error_test > tol || is.nan(error_test) ||
      behav_mixture$functionEvaluations < 300)
    return(eps_min)
  return(behav_mixture$integral)  
}

integrate_density_vec <- Vectorize(FUN = integrate_density, 
                                   vectorize.args = c("rt", "choice", "alpha", 
                                                      "nu", "eta", "shape1",
                                                      "shape2", "shape", "scale",
                                                      "rho_db", "rho_dt", "rho_bt"))



# calc_density_integrand <- function(dimensions, rt, choice,
#                                    sigma, alpha, nu, eta, 
#                                    shape1, shape2, 
#                                    shape, scale,
#                                    sqrt_rho, model) {
#   # Purpose: Calculates integrand for a rt density for a given a 
#   # choice and model
#   # Input: numeric vector dimensions, numeric scalars rt, choice, sigma
#   # numeric vector params, character scalar model
#   # Output: numeric scalar integrand
#   # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
#   # 4:lambda, 5:gamma, 6:chi, 8:phi, 9:rho_db, 10:rho_dt, 11:rho_bt, 12:omega
#   
#   delta_trans <- dimensions[1]
#   beta <- dimensions[2]
#   t_nd <- dimensions[3]
#   delta <- delta_trans / (1 - delta_trans ^ 2)
#   
#   behav_dens <- dwiener(q = rt,
#                         alpha = alpha / sigma, 
#                         tau = t_nd, 
#                         beta = beta, 
#                         delta = delta / sigma,
#                         resp = choice)
# #   if (choice == "lower") rt <- -rt
# #   wiener <- dwiener_cpp(q = rt,
# #                         alpha = alpha / sigma, 
# #                         tau = t_nd, 
# #                         beta = beta, 
# #                         delta = delta / sigma,
# #                         resp = choice)
# #   if (round(behav_dens, 6) != round(wiener, 6)) { 
# #     print(paste("wiener", round(behav_dens, 6), round(wiener, 6)))
# #     print(c(rt, alpha, sigma, delta, beta, t_nd,
# #     nu, eta, shape1, shape2, shape, scale))
# #   }
#   if (behav_dens <= .Machine$double.xmin || is.infinite(behav_dens)) {
#     behav_dens <- .Machine$double.xmin
#   #  wiener <-.Machine$double.xmin
#   }
#   
#   param_logdens <- calc_logdensity(delta = delta, beta = beta, t_nd = t_nd,
#                                    nu = nu, eta = eta,
#                                    shape1 = shape1, shape2 = shape2,
#                                    shape = shape, scale = scale,
#                                    sqrt_rho = sqrt_rho,
#                                    model = model)
# #   dens <- calc_logdensity_c(delta = delta, beta = beta, t_nd = t_nd,
# #                             nu = nu, eta = eta,
# #                             shape1 = shape1, shape2 = shape2,
# #                             shape = shape, scale = scale,
# #                             sqrt_rho = sqrt_rho,
# #                             model = model)
# #   if (round(param_logdens, 6) != round(dens, 6)) {
# #     print(paste("dens", round(param_logdens, 6), round(dens, 6)))
# #     print(c(rt, alpha, sigma, delta, beta, t_nd,
# #             nu, eta, shape1, shape2, shape, scale))
# #   }
#   jacobian <- (1 + delta_trans ^ 2) / (1 - delta_trans ^ 2) ^ 2 
#   integrand <- exp(sum(log(behav_dens), param_logdens, log(jacobian)))
# #  intg <- exp(sum(log(wiener), dens, log(jacobian)))
#   if (is.infinite(integrand)) {
#     integrand <- .Machine$double.xmin
# #    intg <- .Machine$double.xmin
#   }
# #   if (round(integrand, 6) != round(intg, 6)) {
# #     print(paste("integrand", integrand, intg))
# #   }
#   return(integrand)
# }
# 
# integrate_density <- function(rt, choice, sigma, alpha, nu, eta,
#                               shape1, shape2, shape, scale, sqrt_rho, model, 
#                               tol, maxEval) {
#   # Integrates choice probability with respect to parameters' density
#   # Input: numerical scalars rt, choice, sigma, alpha, nu, eta,
#   # shape1, shape2, shape, scale, rho_db, rho_dt, rho_bt, character scalar model, 
#   # numerical scalar tol, integer scalar maxEval
#   # Output: numerical scalar behav_mixture
#   
#   eps_neg <- .Machine$double.neg.eps
#   eps_pos <- .Machine$double.eps
#   eps_min <- .Machine$double.xmin
#   lower_lim <- c(-1 + eps_pos, eps_min, eps_min)
#   upper_lim <- c(1 - eps_neg, 1 - eps_neg, rt)
#   
#   behav_mixture <- integrator(f = calc_density_integrand_c,  
#                               lowerLimit = lower_lim, 
#                               upperLimit = upper_lim, 
#                               rt, choice, sigma, 
#                               alpha, nu, eta, 
#                               shape1, shape2, shape, scale,  
#                               sqrt_rho, model,  
#                               tol = tol, fDim = 1, maxEval = maxEval)
# #   error_test <- behav_mixture$error / behav_mixture$integral
# #   if (error_test > tol || is.nan(error_test) || 
# #       behav_mixture$functionEvaluations < 1.5e3) {
# #     tol <- tol * 1e-1
# #     maxEval <- maxEval * 4
# #     behav_mixture <- integrator(f = calc_density_integrand_c,  
# #                                 lowerLimit = lower_lim, 
# #                                 upperLimit = upper_lim, 
# #                                 rt, choice, sigma, 
# #                                 alpha, nu, eta, 
# #                                 shape1, shape2, shape, scale,  
# #                                 sqrt_rho, model,
# #                                 tol = tol, fDim = 1, maxEval = maxEval)
# #   }
#   return(behav_mixture$integral)  
# }
# 
# integrate_density_vec <- Vectorize(FUN = integrate_density, 
#                                    vectorize.args = c("rt", "choice", "alpha", 
#                                                       "nu", "eta", "shape1",
#                                                       "shape2", "shape", "scale"))


