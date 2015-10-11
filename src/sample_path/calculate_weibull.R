
weibull <- function(bright, lower, upper, scale, shape) {
  # Purpose: calculates weibull function that connects
  #          brightness level to drift rates
  # Inputs: numeric vector bright, numeric scalars lower, upper, scale, shape
  # Outputs: numeric vector drifts

  drifts <- lower + 
    (upper - lower) * 
    (1 - exp((-bright ^ shape) / (scale ^ shape)))
  drifts
}

