
sample_demcmc <- function(chain_idx, theta_idx, posterior_all,
                          history_n, tune_constant) {
  # Purpose: proposal mechanism according to DE-MCMCz
  # Input: integer vectors chain_idx, theta_idx, numeric matrix posterior_all,
  # integer scalar history_n, numeric scalar tune_constant
  # Output: numeric vector theta_new
  
  theta_old <- posterior_all[chain_idx, theta_idx]
  while (TRUE) {
    diff_idx <- sample(seq_len(history_n), 2)
    diff_vectors <- posterior_all[diff_idx, theta_idx]
    theta_new <- theta_old + 
      tune_constant * 
      (diff_vectors[1, ] - diff_vectors[2, ]) +
      runif(length(theta_idx), -1e-3, 1e-3)
    if (check_support(theta_new, model))
      break
  }
  return(theta_new)
}

combine_data <- function(behav_data, theta) {
  # Purpose: forms a data matrix with observations and parameters
  # Input: numeric data.frame behav_data, numeric vector theta
  # Output: numeric data.frame dat_mat
  
  data_mat <- bind_rows(
    transmute(filter(behav_data, instr == 0), rt = rt, choice = resp,
              alpha = theta[1], 
              nu = weibull(prop, theta[2], theta[3], theta[4], theta[5]),
              eta = theta[15], lambda = theta[6], gamma = theta[16],
              chi = theta[7], phi = theta[17], rho_db = theta[18],
              rho_dt = theta[19], rho_bt = theta[20]),
    transmute(filter(behav_data, instr == 1), rt = rt, choice = resp,
              alpha = theta[8], 
              nu = weibull(prop, theta[9], theta[10], theta[11], theta[12]),
              eta = theta[15], lambda = theta[13], gamma = theta[16],
              chi = theta[14], phi = theta[17], rho_db = theta[18],
              rho_dt = theta[19], rho_bt = theta[20]))
  return(data_mat)
}

update_chains <- function(train_data, posterior_all, chain_n,
                          theta_n, history_n, cores, chunk_n) {
  # Purpose: parallelizes simulation of the DEMCMC transition kernel for n chains
  # Input: numeric matrix posterior, integer scalars chain_n, theta_n
  # Output: numeric matrix posterior_updated
  
  theta_idx <- -(theta_n + 1)
  chain_idx <- seq_len(history_n) %>% rev %>% "["(seq_len(chain_n)) %>% rev
  tune_constant <- 2.38 / sqrt(2 * theta_n)
  theta_proposal <- apply(matrix(chain_idx, nrow = chain_n), 1,
                          sample_demcmc, theta_idx, posterior_all,
                          history_n, tune_constant)
  cat(paste(c("proposed chains", theta_proposal)), 
      file = "results/fitting/progress-log-ind-test.txt",
      sep = "\n",
      append = TRUE)
  
  data_mat <- apply(theta_proposal, 1, combine_data, behav_data = train_data)
  joint_new <- joint_logdensity(data_mat, theta_proposal, model,
                                chain_n, cores, chunk_n)
  joint_old <- posterior_all[chain_idx, -theta_idx]
  mh_rule <- log(runif(chain_n)) < joint_new - joint_old
  
  theta_updated <- posterior_all[chain_idx, theta_idx]
  theta_updated[mh_rule, ] <- theta_proposal[mh_rule]
  joint_updated <- joint_old
  joint_updated[mh_rule] <- joint_new[mh_rule]
  posterior_updated <- cbind(theta_updated, joint_updated)
  return(posterior_updated)
}
