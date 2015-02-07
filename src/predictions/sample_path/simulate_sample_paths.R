
simul_paths <- function(model, smpl_size, seed = 2132326000,
                        sigma = 1, time_unit = 1e-3, cores = 1) {
  # Simulates sample paths in parallel for different copula models 
  # Takes a model string, simulation parameter numerical scalars and 
  # returns a list of lists of numeric vectors
  library(package = "doParallel")
  library(package = "doRNG")
  source("src/predictions/sample_path/wiener_parameters.R")
  source("src/predictions/sample_path/combine_parameters.R")
  source("src/predictions/sample_path/simulate_parameters.R")
  source("src/predictions/sample_path/simulate_rndwalk.R")
  
  # Maintains progress log during the simulation
  name <- paste0("results/sample_path/progress-log-", Sys.time(), ".txt")
  writeLines(text = "", con = name)
  
  ind_param <- combine_param(nu = nu, wiener = wiener,
                             rho = rho, omega = omega)
  set.seed(seed)
  registerDoParallel(cores = cores)
  results <- foreach(params = iter(obj = ind_param,
                                   by = 'row')) %dorng% {
             cat(paste(1, "\n"), file = name, append = TRUE)
             trial_param <- smpl_param(params = params,
                                       smpl_size = smpl_size,
                                       model = model)
             rndwalk_vec(smpl_size = 1, 
                         delta = trial_param$delta, 
                         sigma = sigma, trial_param$beta, 
                         trial_param$alpha,
                         low_bound = 0, time_unit = time_unit)
  }
  file.remove(name)
  return(results)
}

