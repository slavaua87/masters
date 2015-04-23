

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
  if (behav_dens <= 0) 
    behav_dens <- .Machine$double.xmin
  if (is.infinite(behav_dens)) 
    behav_dens <- .Machine$double.xmax
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
                              lambda, gamma, chi, phi, rho_db,
                              rho_dt, rho_bt, model, 
                              tol = 1e-2, maxEval = 0) {
  # Integrates choice probability with respect to parameters' density
  # Input: numerical scalars rt, choice, sigma, alpha, nu, eta,
  # lambda, gamma, chi, phi, rho_db, rho_dt, rho_bt, character scalar model, 
  # numerical scalar tol, integer scalar maxEval
  # Output: numerical scalar behav_mixture
    
  if (choice == 1)
    choice = "upper"
  if (choice == 2)
    choice = "lower"
  
  shape1 <- ((1 - lambda) / gamma ^ 2 -
               1 / lambda) * lambda ^ 2
  shape2 <- shape1 * (1 / lambda - 1)

  shape <- (chi / phi) ^ 2
  scale <- phi ^ 2 / chi
  
  if (model == "normal") {
    rho <- matrix(c(1, rho_db, rho_dt,
                    rho_db, 1, rho_bt,
                    rho_dt, rho_bt, 1),
                  byrow = TRUE, ncol = 3)
    sqrt_rho <- cholesky(x = rho)
  }
  eps_neg <- .Machine$double.neg.eps
  eps_pos <- .Machine$double.eps
  lower_lim <- c(-1 + eps_pos, 0 + eps_pos, .Machine$double.xmin)
  upper_lim <- c(1 - eps_neg, 1 - eps_neg, rt)
  
  timer <- proc.time()
  behav_mixture <- adaptIntegrate(f = calc_density_integrand,  
                                  lowerLimit = lower_lim, 
                                  upperLimit = upper_lim, 
                                  rt, choice, sigma, 
                                  alpha, nu, eta,
                                  shape1, shape2, shape, scale, 
                                  sqrt_rho, model, tol = tol, 
                                  fDim = 1, maxEval = maxEval)
  timer <- proc.time() - timer
  return(behav_mixture$integral)
}

integrate_density_vec <- Vectorize(FUN = integrate_density, 
                                   vectorize.args = c("rt", "choice", "alpha", 
                                                      "nu", "eta", "lambda", 
                                                      "gamma", "chi", "phi",
                                                      "rho_db", "rho_dt", "rho_bt"))

