
# Calculates weibull function on a constrained domain
weibull <- function(bright, lower, upper, shape, scale) {
  # bright - covariate
  # lower, upper, scale, shape - parameters
  drifts <- -lower + (upper + lower) * (1 - 
            exp(-((bright / scale) ^ shape)))
  return(drifts)
}

