
# Calculates weibull function on a constrained domain
weibull <- function(bright, lower, upper, scale, shape) {
  # bright - covariate
  # lower, upper, scale, shape - parameters
  drifts <- lower + (upper - lower) * (1 - 
            exp((-bright ^ shape) / (scale ^ shape)))
  return(drifts)
}

