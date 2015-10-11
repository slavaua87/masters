
backsub <- function(r, x, upper.tri = TRUE, transpose = FALSE) {
  y <- .Internal(call = backsolve(r, x, 3, upper.tri, transpose))
  return(y)
}

cholesky <- function(x, pivot = FALSE, LINPACK = FALSE, tol = -1) {
  .Internal(call = La_chol(x, pivot, tol))
}

normal_logcopula <- function(cop_value, sqrt_rho) {
  x <- qnorm(p = cop_value)
  marginals <- dnorm(x)
  zero_test <- marginals == 0
  if (any(zero_test)) {
    #print(paste("marginals", marginals))
    marginals[zero_test] <- .Machine$double.xmin
  }
  sqrt_rho_x <- backsub(r = sqrt_rho, x = x, transpose = TRUE)
  rho_x <- sum(sqrt_rho_x ^ 2)
  rho_diag <- as.numeric(sqrt_rho)[c(1, 5, 9)]
  zero_test <- rho_diag == 0
  if (any(zero_test)) {
   # print(paste("rho_diag", rho_diag))
    rho_diag[zero_test] <- .Machine$double.xmin
  }
  logdens <- (-1.5 * log(2 * pi) - 
              sum(log(rho_diag)) - 
              0.5 * rho_x) - sum(log(marginals))  
  return(logdens)
}

calc_logdensity <- function(delta, beta, t_nd, nu, eta,
                            shape1, shape2, shape, scale,
                            sqrt_rho, model) {
  # Calculates probability density for independent and normal
  # Input is a numeric scalars delta, beta, t_nd, 
  # numeric data.frame params, character scalar model
  # Output is a numeric scalar dens
  # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
  # 4:lambda, 5:gamma, 6:chi, 7:phi, 8:rho_db, 9:rho_dt, 10:rho_bt, 11:omega
    
  marginal_dens <- c(
    dnorm(x = delta, mean = nu, sd = eta), 
    dbeta(beta, shape1 = shape1, shape2 = shape2),
    dgamma(x = t_nd, shape = shape, scale = scale))
  zero_test <- marginal_dens == 0
  if (any(zero_test)) {
    #print(paste("cop marginals", marginal_dens))
    marginal_dens[zero_test] <- .Machine$double.xmin
  }
  inf_test <- is.infinite(marginal_dens)
  if (any(inf_test)) {
   # print(paste("cop marginals", marginal_dens))
    marginal_dens[inf_test] <- .Machine$double.xmax
  }
  marginal_logdens <- sum(log(marginal_dens))
  if (model == "independent") {
    return(marginal_logdens)
  }
  
  cop_value <- c(pnorm(q = delta, mean = nu, sd = eta),
                 pbeta(q = beta, shape1 = shape1, shape2 = shape2),
                 pgamma(q = t_nd, shape = shape, scale = scale))
  zero_test <- cop_value == 0
  if (any(zero_test)) {
    #print(paste("copula", cop_value))
    cop_value[zero_test] <- .Machine$double.xmin
  }
  one_test <- cop_value == 1
  if (any(one_test)) {
    #print(paste("copula", cop_value))
    cop_value[one_test] <- 1 - .Machine$double.neg.eps
  }
  
  joint_logdens <- normal_logcopula(cop_value = cop_value, sqrt_rho)
  #print(c(marginal_logdens, joint_logdens))
  logdens <- sum(joint_logdens, marginal_logdens)
  return(logdens)
}

