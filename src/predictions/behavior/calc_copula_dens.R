
library("copula")
library("tmvtnorm")

calc_density <- function(drift, bias, nondec, params, model) {
  # Calculates probability density for independent, normal and 
  # t copula-based models
  # Input is a numeric scalar drift, numeric scalar bias, numeric scalr nondec,
  # numeric data.frame params, character scalar model
  # Output is a numeric scalar dens
  
  shape1 <- ((1 - params[1, "lambda"]) / params[1, "gamma"] ^ 2 -
               1 / params[1, "lambda"]) * params[1, "lambda"] ^ 2
  shape2 <- shape1 * (1 / params[1, "lambda"] - 1)
  
  marginal_dens <- dnorm(x = drift, mean = params[1, "nu"], 
                         sd = params[1, "eta"]) *
    dbeta(bias, shape1 = shape1, shape2 = shape2) * 
    dtmvnorm.marginal(xn = nondec, n = 1, mean = params[1, "chi"],
                      sigma = params[1, "phi"], lower = 0, upper = Inf)
  
  if (model == "independent") {
    return(marginal_dens)
  }
  
  else {
    cop <- ellipCopula(family = model, 
                       param = as.numeric(params[1, 
                                          c("rho_db", "rho_dt", "rho_bt")]), 
                       dim = 3, dispstr = 'un', df = params[1, "omega"])
    cop_value <- c(pnorm(q = drift, mean = params[1, "nu"], 
                         sd = params[1, "eta"]),
                   pbeta(q = bias, shape1 = shape1, shape2 = shape2),
                   ptmvnorm.marginal(xn = nondec, n = 1, mean = params[1, "chi"],
                                     sigma = params[1, "phi"], 
                                     lower = 0, upper = Inf))
    if (any(cop_value == 0)) {
      cop_value[cop_value == 0] <- cop_value[cop_value == 0] + 
        .Machine$double.eps
    }
    if (any(cop_value == 1)) {
      cop_value[cop_value == 1] <- cop_value[cop_value == 1] -
        .Machine$double.neg.eps
    }
    joint_dens <- dCopula(u = cop_value, copula = cop)
    dens <- joint_dens * marginal_dens
    return(dens)
  }
}

