

source("src/fitting/initiate_sampler.R")
source("src/fitting/sample_demcmc.R")

sample_posterior <- function(param_n, chain_n, draw_n, 
                             cores = 1, seed = 1800968452) {
  # Purpose: samples from the posterior using parallelized demcmc transitios
  # Input: integer scalars chain_n, draw_n, cores, seed
  # Output: numeric array posterior
  
  posterior <- array(data = 0, dim = c(chain_n, param_n, draw_n))
  
  set.seed(seed = seed)
  posterior[, , 1] <- initialize_chains(model, chain_n)
  
  registerDoParallel(cores = cores)
  for (draw in 2:draw_n)) {
    posterior[, , draw] <- sample_demcmc(posterior[, , draw - 1])
  }
  return(posterior)
}


