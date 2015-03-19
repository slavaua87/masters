
library("copula")
library("tmvtnorm")

backsub <- function(r, x, upper.tri = TRUE, transpose = FALSE) {
  y <- .Internal(call = backsolve(r, x, 3, upper.tri, transpose))
  return(y)
}

cholesky <- function(x, pivot = FALSE, LINPACK = FALSE, tol = -1) {
  .Internal(call = La_chol(x, pivot, tol))
}

t_logcopula <- function(cop_value, sqrt_rho, df) {
  x <- qt(p = cop_value, df = df)
  sqrt_rho_x <- backsub(r = sqrt_rho, x = x, transpose = TRUE)
  rho_x <- sum(sqrt_rho_x ^ 2)
  logdens <- (lgamma((3 + df) / 2) - 
                (lgamma(df / 2) + sum(log(as.numeric(sqrt_rho)[c(1, 5, 9)])) + 
                   1.5 * log(pi * df)) - 0.5 * (df + 3) * 
                log1p(rho_x / df)) - sum(log(dt(x, df)))  
  return(logdens)
}

normal_logcopula <- function(cop_value, sqrt_rho) {
  x <- qnorm(p = cop_value)
  sqrt_rho_x <- backsub(r = sqrt_rho, x = x, transpose = TRUE)
  rho_x <- sum(sqrt_rho_x ^ 2)
  logdens <- (-1.5 * log(2 * pi) - 
              sum(log(as.numeric(sqrt_rho)[c(1, 5, 9)])) - 
              0.5 * rho_x) - sum(log(dnorm(x)))  
  return(logdens)
}


calc_logdensity <- function(delta, beta, t_nd, params,
                            shape1, shape2, shape, scale,
                            sqrt_rho, model) {
  # Calculates probability density for independent, normal and 
  # t copula-based models
  # Input is a numeric scalars delta, beta, t_nd, 
  # numeric data.frame params, character scalar model
  # Output is a numeric scalar dens
  # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
  # 4:lambda, 5:gamma, 6:chi, 7:phi, 8:rho_db, 9:rho_dt, 10:rho_bt, 11:omega
    
  marginal_dens <- c(
    dnorm(x = delta, mean = params[2], sd = params[3]), 
    dbeta(beta, shape1 = shape1, shape2 = shape2),
    dgamma(x = t_nd, shape = shape, scale = scale))
  if (any(marginal_dens == 0)) {
    marginal_dens[marginal_dens == 0] <- marginal_dens[
      marginal_dens == 0] + .Machine$double.eps
  }
  marginal_logdens <- sum(log(marginal_dens))
  if (model == "independent") {
    return(marginal_logdens)
  }
  
  cop_value <- c(pnorm(q = delta, mean = params[2], 
                       sd = params[3]),
                 pbeta(q = beta, shape1 = shape1, shape2 = shape2),
                 pgamma(q = t_nd, shape = shape, 
                        scale = scale, lower.tail = TRUE))
  
  if (any(cop_value <= sqrt(.Machine$double.xmin))) {
    cop_value[cop_value <= sqrt(.Machine$double.xmin)] <- 
      cop_value[cop_value <= sqrt(.Machine$double.xmin)] + .Machine$double.eps
  }
  if (any(cop_value == 1)) {
    cop_value[cop_value == 1] <- cop_value[cop_value == 1] -
      .Machine$double.neg.eps
  }
  
  if (model == "normal") {
    joint_logdens <- normal_logcopula(cop_value = cop_value,
                                      sqrt_rho)
  }
  if (model == "t") {
    joint_logdens <- t_logcopula(cop_value = cop_value,
                                 sqrt_rho,
                                 df = params[11])
  }
  logdens <- sum(joint_logdens, marginal_logdens)
  return(logdens)
}

