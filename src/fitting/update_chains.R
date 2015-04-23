

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

update_chains <- function(train_data, posterior_cur, chain_n,
                          theta_n, cores, chunk) {
  # Purpose: parallelizes simulation of the DEMCMC transition kernel for n chains
  # Input: numeric matrix posterior_cur, integer scalars chain_n, theta_n
  # Output: numeric matrix posterior_updated
  
  theta_idx <- -(theta_n + 1)
  chain_idx <- nrow(posterior_cur)
  tune_constant <- .1 ^ 2 / theta_n
  theta_proposal <- apply(posterior_cur, 3, sample_am, 
                          chain_idx, theta_idx, theta_n,
                          tune_constant) %>% t
  new_idx <- theta_proposal[, -theta_idx] == 1
  if (sum(new_idx) == 0)
    return(posterior_cur[chain_idx, , ])
  if (sum(new_idx) == 1) {
    data_mat <- combine_data(train_data, theta_proposal[new_idx, theta_idx])
  } 
  else {
    data_mat <- theta_proposal[new_idx, theta_idx] %>%
      apply(1, combine_data, behav_data = train_data) %>% 
      rbind_all
  }  

  joint_old <- posterior_cur[chain_idx, -theta_idx, ]
  joint_new <- numeric(chain_n)
  joint_new[new_idx] <- joint_logdensity(data_mat, 
                                         theta_proposal[new_idx, theta_idx],
                                         model, sum(new_idx), cores, chunk)
  joint_new[!new_idx] <- joint_old[!new_idx]
  mh_rule <- runif(chain_n) %>% log < joint_new - joint_old
  cat("*** new ***",
      joint_new,
      "*** old ***",
      joint_old,
      file = "results/fitting/progress-likelihoods.txt",
      sep = "\n",
      append = TRUE)

  theta_updated <- posterior_cur[chain_idx, theta_idx, ] %>% t
  theta_updated[mh_rule, ] <- theta_proposal[mh_rule, theta_idx]
  joint_updated <- joint_old
  joint_updated[mh_rule] <- joint_new[mh_rule]
  posterior_updated <- cbind(theta_updated, joint_updated) %>% t
  #write.table(rbind(posterior_cur[chain_idx, theta_idx, 3],
  #            posterior_updated[3, ]),
  #            "results/fitting/thetas_last.txt", append = TRUE)
  return(posterior_updated)
}
