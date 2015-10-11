
weibull <- function(bright, lower, upper, shape, scale) {
  # Purpose: calculates weibull function that connects
  #          brightness level to drift rates
  # Inputs: numeric vector bright, numeric scalars lower, 
  #          upper, scale, shape
  # Output: numeric vector drifts
  
  drifts <- lower + 
    (upper - lower) * 
    (1 - exp((-bright ^ shape) / (scale ^ shape)))
  drifts
}
