
demcmc_z <- function(chain_idx, theta_idx, posterior_all,
                          history_n, tune_constant) {
  # Purpose: proposal mechanism according to DE-MCMCz
  # Input: integer vectors chain_idx, theta_idx, numeric matrix posterior_all,
  # integer scalar history_n, numeric scalar tune_constant
  # Output: numeric vector theta_new
  
  theta_old <- posterior_all[chain_idx, theta_idx]
  while (TRUE) {
    diff_idx <- sample(seq_len(history_n)[-chain_idx], 2)
    diff_vectors <- posterior_all[diff_idx, theta_idx]
    theta_new <- theta_old + 
      tune_constant * 
      (diff_vectors[1, ] - diff_vectors[2, ]) +
      runif(length(theta_idx), -1e-3, 1e-3)
    check <- check_support(theta_new, model)
    if (as.logical(check[1]))
      break
  }
  return(check[-1])
}