
library("copula")
library("tmvtnorm")

backsub <- function(r, x, upper.tri = TRUE, transpose = FALSE) {
  y <- .Internal(backsolve(r, x, 3, upper.tri, transpose))
  return(y)
}

cholesky <- function(x, pivot = FALSE, LINPACK = FALSE, tol = -1) {
  .Internal(La_chol(x, pivot, tol))
}

t_copula <- function(cop_value, sqrt_rho, df) {
  x <- qt(p = cop_value, df = df)
  sqrt_rho_x <- backsub(r = sqrt_rho, x = x, transpose = TRUE)
  rho_x <- sum(sqrt_rho_x ^ 2)
  dens <- exp(lgamma((3 + df) / 2) - 
                (lgamma(df / 2) + sum(log(as.numeric(sqrt_rho)[c(1, 5, 9)])) + 
                   1.5 * log(pi * df)) - 0.5 * (df + 3) * 
                log1p(rho_x / df)) / prod(dt(x, df))  
  return(dens)
}

normal_copula <- function(cop_value, sqrt_rho) {
  x <- qnorm(p = cop_value)
  sqrt_rho_x <- backsub(r = sqrt_rho, x = x, transpose = TRUE)
  rho_x <- sum(sqrt_rho_x ^ 2)
  dens <- exp(-1.5 * log(2 * pi) - 
              sum(log(as.numeric(sqrt_rho)[c(1, 5, 9)])) - 
              0.5 * rho_x) / prod(dnorm(x))  
  return(dens)
}


calc_density <- function(delta, beta, t_nd, params,
                         shape1, shape2, trunc_const, sqrt_rho,
                         lower, upper, model) {
  # Calculates probability density for independent, normal and 
  # t copula-based models
  # Input is a numeric scalars delta, beta, t_nd, lower, upper
  # numeric data.frame params, character scalar model
  # Output is a numeric scalar dens
  # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
  # 4:lambda, 5:gamma, 6:chi, 7:phi, 8:rho_db, 9:rho_dt, 10:rho_bt, 11:omega
  
  marginal_dens <- dnorm(x = delta, mean = params[2], 
                         sd = params[3]) *
    dbeta(beta, shape1 = shape1, shape2 = shape2) * 
    dnorm(x = t_nd, mean = params[6], sd = params[7]) /
    (trunc_const[2] - trunc_const[1])
  
  if (model == "independent") {
    return(marginal_dens)
  }
  
  cop_value <- c(pnorm(q = delta, mean = params[2], 
                       sd = params[3]),
                 pbeta(q = beta, shape1 = shape1, shape2 = shape2),
                 (pnorm(q = t_nd, mean = params[6],
                        sd = params[7], lower.tail = TRUE) - 
                    trunc_const[1]) / (trunc_const[2] - trunc_const[1]))
  
  if (any(cop_value == 0)) {
    cop_value[cop_value == 0] <- cop_value[cop_value == 0] + 
      .Machine$double.eps
  }
  if (any(cop_value == 1)) {
    cop_value[cop_value == 1] <- cop_value[cop_value == 1] -
      .Machine$double.neg.eps
  }
  
  if (model == "normal") {
    joint_dens <- normal_copula(cop_value = cop_value,
                                sqrt_rho)
  }
  if (model == "t") {
    joint_dens <- t_copula(cop_value = cop_value,
                           sqrt_rho,
                           df = params[11])
  }
  dens <- joint_dens * marginal_dens
  return(dens)
}

