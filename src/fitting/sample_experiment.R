

sample_behavior <- function(params, smpl_size, model, sigma = .1) {
  # Purpose: samples response times and responses
  # Inputs: double vector params, integer scalar smpl_size,
  #         character scalar model, double scalar sigma
  # Output: double matrix behav_smpl
  
  trial_param <- smpl_param(params = params,
                            smpl_size = smpl_size,
                            model = model)

  behav_smpl <- smpl_rts(alpha = trial_param$alpha,
                         tau = trial_param$t_nd,
                         beta = trial_param$beta,
                         delta = trial_param$delta, 
                         sigma = sigma)
  behav_smpl
}

sample_experiment <- function(theta, model, prop, smpl_size) {
  # Purpose: samples data from a brightness discrimination task
  # Inputs: double vector theta, character scalar model,
  #         double vector prop, integer scalar smpl_size
  # Output: data_matrix data_mat
  
  prop_n <- length(prop)
  instr <- c(1, 0)
  theta[1:20] <- exp(theta[1:20])
  
  exper_theta <- 
    bind_rows(data_frame(alpha = theta[1], 
                         nu = weibull(bright = prop, lower = theta[2],             
                                      upper = theta[3], shape = theta[4],
                                      scale = theta[5]),
                         eta = theta[6],
                         shape1 = theta[7],
                         shape2 = theta[8],
                         shape = theta[9],
                         scale = theta[10]),
              data_frame(alpha = theta[11], 
                         nu = weibull(bright = prop, lower = theta[12],
                                      upper = theta[13], shape = theta[14],
                                      scale = theta[15]),
                         eta = theta[16],
                         shape1 = theta[17],
                         shape2 = theta[18],
                         shape = theta[19],
                         scale = theta[20]))
  if (model == "normal") {
    exper_theta <- 
      bind_cols(exper_theta, 
                data_frame(rho_db = rep(x = theta[21],                  
                                        times = 2 * prop_n),
                           rho_dt = rep(x = theta[22], 
                                        times = 2 * prop_n),
                           rho_bt = rep(x = theta[23], 
                                        times = 2 * prop_n)))
  }
  
  behav_data <- apply(X = exper_theta, MARGIN = 1, FUN = sample_behavior,
                      smpl_size = smpl_size, model = model) %>% bind_rows
  colnames(behav_data) <- c("rt", "resp")
  prop_data <- data_frame(prop = rep(prop, each = smpl_size) %>%
                            rep(times = 2))
  # 1 for accuracy, 0 speed
  instr_data <- data_frame(instr = rep(instr, each = prop_n * smpl_size))
  data_mat <- bind_cols(behav_data, prop_data, instr_data)
  data_mat
}


