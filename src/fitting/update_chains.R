

combine_data <- function(theta) {
  # Purpose: forms a data matrix with observations and parameters
  # Input: data_frame train_data, double vector theta
  # Output: data_frame dat_mat
  
  theta[c(1, 4:11, 14:20)] <- exp(theta[c(1, 4:11, 14:20)])
  
  acc <- data_frame(rt = filter(train_data, instr == 1) %>%
                    dplyr::select(rt) %>%
                    unlist(use.names = FALSE),
                  choice = filter(train_data, instr == 1) %>% 
                    dplyr::select(resp) %>% 
                    unlist(use.names = FALSE),
                  alpha = theta[1], 
                  nu = weibull(filter(train_data, instr == 1) %>% 
                                 dplyr::select(prop) %>% 
                                 unlist(use.names = FALSE),
                               theta[2], theta[3], theta[4], theta[5]),
                  eta = theta[6], shape1 = theta[7], shape2 = theta[8],
                  shape = theta[9], scale = theta[10], rho_db = theta[21],
                  rho_dt = theta[22], rho_bt = theta[23])
  spd <- data_frame(rt = filter(train_data, instr == 0) %>%
                      dplyr::select(rt) %>%
                      unlist(use.names = FALSE),
                    choice = filter(train_data, instr == 0) %>%
                      dplyr::select(resp) %>%
                      unlist(use.names = FALSE),
                    alpha = theta[11], 
                  nu = weibull(filter(train_data, instr == 0) %>%
                                 dplyr::select(prop) %>% 
                                 unlist(use.names = FALSE),
                               theta[12], theta[13], theta[14], theta[15]),
                  eta = theta[16], shape1 = theta[17], shape2 = theta[18],
                  shape = theta[19], scale = theta[20], rho_db = theta[21],
                  rho_dt = theta[22], rho_bt = theta[23])
  as.matrix(rbind.data.frame(acc, spd))
}

update_chains <- function() {
  # Purpose: explores the posterior with MCMC
  # Inputs: double matrix posterior_cur, integer scalars chain_n, theta_n
  # Output: double matrix posterior_updated
  e1 <- parent.frame(1)
  
  theta_proposal <- sample_am()
  joint_old <- e1$posterior_all[e1$chain_idx, theta_n + 1]
  
  if (theta_proposal[theta_n + 1] == 0) 
    return(e1$posterior_all[e1$chain_idx, ])
  else 
    data_mat <- combine_data(theta_proposal[e1$theta_idx])
  
  joint_new <- integral$joint_logdensity_cpp(data_mat, 
                                             theta_proposal[e1$theta_idx],
                                             model, thread_n, chunk_n, 
                                             tol, maxEvals)

  alpha_log <- min(joint_new - joint_old, 0)
  tune_updated <- e1$posterior_all[e1$chain_idx, theta_n + 2]
  if (e1$chain_idx > 200)
    tune_updated <- exp(log(tune_updated) + 
      1 / e1$chain_idx ^ (1 / 2) * 
      (exp(alpha_log) - alpha))  
  
  if (log(runif(1)) < alpha_log)
    c(theta_proposal[e1$theta_idx], joint_new, tune_updated)
  else
    c(e1$posterior_all[e1$chain_idx, e1$theta_idx], joint_old, tune_updated)
}
