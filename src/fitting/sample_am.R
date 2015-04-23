

sample_am <- function(posterior_cur, chain_idx, theta_idx, 
                      theta_n, tune_constant) {
  # Purpose: proposal mechanism according to adaptive metropolis
  # Input: integer vectors chain_idx, theta_idx, numeric matrix posterior_all,
  # numeric scalar tune_constant
  # Output: numeric vector theta_new
  
  if (chain_idx == 1)
    cov_estimate <- tune_constant * diag(theta_n)
  else
    cov_estimate <- cov(posterior_cur[, theta_idx]) + tune_constant * diag(theta_n)
  theta_old <- posterior_cur[chain_idx, theta_idx]
  theta_new <- rmvnorm(1, theta_old, cov_estimate)
  theta_reflected <- reflect(theta_new, model)
  check <- check_support(theta_reflected, model)
#  write.table(rbind(theta_old, theta_new, theta_reflected),
#              "results/fitting/thetas_proposed.txt", append = TRUE)
  if (check)
    return(c(theta_reflected, 1))
  else
    return(c(theta_old, 0))
}