source(file = "src/fitting/calculate_weibull.R")
source(file = "src/fitting/simulate_wiener_parameters.R")
source(file = "src/fitting/simulate_rndwalk_rts.R")

sample_behavior <- function(params, smpl_size, model, sigma = .1) {
  trial_param <- smpl_param(params = params,
                            smpl_size = smpl_size,
                            model = model)

  behav_smpl <- smpl_rts(alpha = trial_param$alpha,
                         tau = trial_param$t_nd,
                         beta = trial_param$beta,
                         delta = trial_param$delta, 
                         sigma = sigma)
  return(behav_smpl)
}

sample_experiment <- function(theta, model, prop, smpl_size) {
  
  prop_n <- length(prop)
  instr <- c(0, 1)
  
  exper_theta <- 
    bind_rows(data.frame(alpha = theta[1], 
                         nu = weibull(bright = prop, lower = theta[2],             
                                      upper = theta[3], shape = 9,
                                      scale = .5),
                         eta = theta[15],
                         lambda = theta[6],
                         gamma = theta[16],
                         chi = theta[7],
                         phi = theta[17]),
              data.frame(alpha = theta[8], 
                         nu = weibull(bright = prop, lower = theta[9],
                                      upper = theta[10], shape = 14,
                                      scale = .6),
                         eta = theta[15],
                         lambda = theta[13],
                         gamma = theta[16],
                         chi = theta[14],
                         phi = theta[17]))
  if (model == "normal") {
    exper_theta <- 
      bind_cols(exper_theta, 
                data.frame(rho_db = rep(x = theta[18],                  
                                        times = 2 * prop_n),
                           rho_dt = rep(x = theta[19], 
                                        times = 2 * prop_n),
                           rho_bt = rep(x = theta[20], 
                                        times = 2 * prop_n)))
  }
  
  behav_data <- apply(X = exper_theta, MARGIN = 1, FUN = sample_behavior,
                      smpl_size = smpl_size, model = model) %>% bind_rows
  colnames(behav_data) <- c("rt", "resp")
  prop_data <- data.frame(prop = rep(x = prop, each = smpl_size) %>% rep(times = 2))
  instr_data <- data.frame(instr = rep(x = c(0, 1), each = prop_n * smpl_size))
  data_mat <- bind_cols(behav_data, prop_data, instr_data)
  return(data_mat)
}


