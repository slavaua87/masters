
# Algorithm to evaluate cumulative FPT distribution of the Wiener process
# conditioned on upper or lower absorbption boundary

prob_absorb <- function(a, z, v) {
  # Calculates probality of the path being absorbed at the lower bound
  # Takes numeric scalars a, z, v parameterizing threshold, starting point
  # drift coefficient
  
  w <- z / a
  
  if (abs(v) < .Machine$double.xmin) {
    prob <- 1 - w
    return(prob)
  }
  
  else {
    prob <- min(1, (1 - exp(-2 * v * a * (1 - w))) / 
                    (exp(2 * v * a * w) - 
                     exp(-2 * v * a * (1 - w))),
                na.rm = TRUE)
    return(prob)
  }
}

iter_large <- function(rt, a, v, z, eps) {
  # Calculates number of iterations for large-time representation
  # Takes numeric scalars rt, a, v, z, eps parameterizing time, threshold,
  # drift coefficient, starting point and approximation precision
  
  w <- z / a
  
  cond1 <- sqrt(a ^ 2 / (rt * pi ^ 2)) 
  cond2 <- sqrt(pmax(1, -2 * a ^ 2 / (rt * pi ^ 2) * 
                    (log(pmax(.Machine$double.xmin, eps * pi * rt *
                    (v ^ 2 + pi ^ 2 / a ^ 2) / 2)) + 
                     v * a * w + v ^ 2 * rt / 2)))
  n_iter <- ceiling(pmax(1, cond1, cond2))
  return(n_iter) 
}

iter_small <- function(rt, a, v, z, eps) {
  # Calculates number of iterations for small-time representation
  # Takes numeric scalars rt, a, v, z, eps parameterizing time, threshold,
  # drift coefficient, starting point and approximation precision
  
  w = z / a
  
  if (v > 0) {
    n_iter <- iter_small(rt, a, -v, z, eps * exp(-2 * v * a * w))
    return(n_iter)
  }

  if (v < 0) {  
    cond1 <- w - 1 + log(pmax(.Machine$double.xmin, eps * (1 - exp(2 * v * a)) /
                              2)) / (2 * v * a)
    cond2 = (.535 * sqrt(2 * rt) + v * rt + a * w) / (2 * a)
    cond3 = w / 2 - sqrt(rt) * 
      qnorm(pmax(.Machine$double.xmin, 
                 pmin(1 - .Machine$double.neg.eps,
                 eps * a * exp(v ^ 2 * rt / 2 + v * a * w ) /
                 (.3 * sqrt(2 * pi * rt))))) /
      (2 * a)
    n_iter <- ceiling(pmax(cond1, cond2, cond3, 1)) 
    return(n_iter)
  }
  
  if (abs(v) < .Machine$double.xmin) {  
    cond <- w / 2 - sqrt(rt) / (2 * a) * 
      qnorm(pmin(1 - .Machine$double.neg.eps,
                 pmax(.Machine$double.xmin, eps / (2 - 2 * w))))
    n_iter <- ceiling(pmax(cond, 1))
    return(n_iter)
  }
}

largecdf_lower <- function(rt, k, a, v, z) {
  # Cumulative distribution approximation of a large-time 
  # representation for lower bound
  # Takes numeric scalars rt, k, a, v, z parameterizing time, 
  # number of iterations, threshold, drift coefficient, starting point
  
  w <- z / a
  norm_const <- prob_absorb(a, z, v)
  n_iter <- seq_len(k) 
  cdf_lower <- pmax(0, norm_const - 2 * pi / a ^ 2 * 
                       exp(-v * a * w - v ^ 2 * rt / 2) *
                       sum(n_iter * sin(pi * n_iter * w) * 
                       exp(-.5 * rt * n_iter ^ 2 * pi ^ 2 / a ^ 2) / 
                       (v ^ 2 + (n_iter * pi / a) ^ 2)))
  return(cdf_lower)
}
  

smallcdf_lower <- function(rt, k, a, v, z) {
  # Cumulative distribution approximation of a small-time 
  # representation for lower bound  
  # Takes numeric scalars rt, k, a, v, z parameterizing time, 
  # number of iterations, threshold, drift coefficient, starting point
  
  w <- z / a
  norm_const <- prob_absorb(a, z, v)
  
  if (abs(v) < .Machine$double.xmin) {  
    n_iter <- seq_len(k)
    cdf_upper <- 2 * sum(pnorm((-2 * a * n_iter - a * w) / sqrt(rt)) -
                         pnorm((-2 * a * n_iter + a * w - 2 * a) / sqrt(rt)))
    return(cdf_upper)   
  }
    
  n_iter <- -floor((k - 1) / 2):ceiling((k - 1) / 2)
  cdf_upper <-  pmax(0, norm_const - sign(v) * 
                     sum(exp(-2 * v * a * w - 2 * v * a * n_iter) * 
                     pnorm(sign(v) * (2 * a * n_iter + a * w - v * rt)
                           / sqrt(rt)) -
                     exp(2 * v * a * n_iter) * pnorm(sign(v) * 
                     (-2 * a * n_iter - a * w - v * rt) / sqrt(rt)))) 
  return(cdf_upper)
}

lower_CDF <- function(rt, a, v, s, z, t_nd, eps) { 
  # Evaluate lower bound CDF
  # Takes numeric vector rt, numeric scalars a, v, s, z, t_nd, eps 
  # parameterizing time, threshold, drift coefficient, 
  # diffusion coefficient, starting point and approximation precision
  
  if (min(rt) < t_nd) return(0)
    #stop("rt must be larger than t_nd")
  
  probs <- numeric(length(rt))
  a <- a / s
  v <- v / s
  z <- z / s
  norm_constant <- prob_absorb(a, z, v)
  decision <- rt - t_nd
  
  n_large <- iter_large(decision, a, v, z, eps)
  n_small <- 10 * iter_small(decision, a, v, z, eps)
  ind_large <- n_large <= n_small
  if (any(ind_large)) {
    probs[ind_large] <- largecdf_lower(decision[ind_large], 
                                       max(n_large[ind_large]), a, v, z)
  }
  if (any(!ind_large)) {
  probs[!ind_large] <- smallcdf_lower(decision[!ind_large],
                                      max(n_small[!ind_large]), a, v, z)
  }
  cdf <- probs / norm_constant
  cdf[is.nan(cdf)] <- 0
  return(cdf)
}

upper_CDF <- function(rt, a, v, s, z, t_nd, eps) {
  # Evaluate upper bound CDF
  # Takes numeric scalars rt, a, v, s, z, eps parameterizing time, threshold,
  # drift coefficient, diffusion coefficient, starting point and 
  # approximation precision
  cdf <- lower_CDF(rt, a, -v, s,  a - z, t_nd, eps)
  return(cdf)
}









  