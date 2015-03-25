
library(package = "dplyr")
library(package = "magrittr")
library(package = "doParallel")
library(package = "doRNG")

source(file = "src/predictions/behavior/wiener_parameters.R")
source(file = "src/predictions/behavior/combine_parameters.R")
source(file = "src/predictions/behavior/simulate_wiener_parameters.R")
source(file = "src/predictions/behavior/simulate_rndwalk_rts.R")
source(file = "src/predictions/behavior/calculate_summary.R")

simul_behavior <- function(model, smpl_size, seed = 1771363045,
                           sigma = .1, cores = 1) {
  # Purpose: simulates reaction times and choices in parallel for
  # different copula models 
  # Input: character scalar model, numerical scalars smpl_size,
  # seed, sigma, cores
  # Ouput: numeric data.frame results
  
  name <- paste0("results/behavior/progress-log-", Sys.time(), ".txt")
  writeLines(text = "", con = name)
  
  ind_param <- combine_param(nu = nu, wiener = wiener,
                             rho = rho, omega = omega)
  
  set.seed(seed = seed)
  registerDoParallel(cores = cores)
  results <- foreach(params = iter(obj = ind_param,
                                   by = 'row')) %dorng% {
             cat(paste(smpl_size, "\n"), file = name, append = TRUE)
             trial_param <- smpl_param(params = params,
                                       smpl_size = smpl_size,
                                       model = model)
             behav_smpl <- smpl_rts(n = 1,
                           alpha = trial_param$alpha,
                           tau = trial_param$t_nd,
                           beta = trial_param$beta,
                           delta = trial_param$delta, 
                           sigma = sigma)
             calc_quant_prob(behav_smpl = behav_smpl)
  }
  
  results <- bind_rows(results) %>%
    mutate(instruction = c("accuracy", "speed")[rep(c(1, 2), each = 60)] %>%
             rep(times = 9), 
           condition = nrow(ind_param) %>% seq_len %>% 
             rep(each = 10), model = model)
  file.remove(name)
  return(results)
}




