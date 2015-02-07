
library(compiler)

rndwalk <- function(smpl_size = 1, delta = .01, sigma = 1, beta = .75, 
                    alpha = 1.3, low_bound = 0, time_unit = 1e-3) {
  # Random walk simulation of a bounded Wiener process
  # Takes scalar process parameters, number of paths to give a list of sample paths
  path_smpl <- vector(mode = "list", length = smpl_size)
  state_unit <- sigma * sqrt(time_unit)
  prob.up <- .5 * (1 + delta * sqrt(time_unit) / sigma)
  if (prob.up > 1) prob.up <- 1
  if (prob.up < 0) prob.up <- 0
  for (draw in seq_len(smpl_size)) {
    path <- vector(mode = "double", length = 1e4)
    path[1] <- position <- beta * alpha
    timer <- 0
    counter <- 2
    while (all(position < alpha, position > low_bound)) {    
      if (rbinom(1, 1, prob.up) == 1) {
        position <- position + state_unit
        path[counter] <- position
      }
      else {
        position <- position - state_unit
        path[counter] <- position
      } 
      timer <- timer + time_unit
      counter <- counter + 1
    }    
    path_smpl[[draw]] <- path[path != 0]
  }
  return(path_smpl)
}

# Compiles just-in-time
rndwalk_cmp <- cmpfun(f = rndwalk)
# Vectorizes with respect to variable parameters
rndwalk_vec <- Vectorize(FUN = rndwalk_cmp, 
                         vectorize.args = list("delta", "beta", "alpha"))





