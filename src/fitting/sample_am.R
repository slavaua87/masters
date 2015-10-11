

sample_am <- function() {
  # Purpose: proposal mechanism according to adaptive metropolis
  # Inputs: integer vectors chain_idx, theta_idx, 
  #         double matrix posterior_all,
  #         double scalar tune_constant (all global)
  # Output: double vector c(theta_reflected, 1) or 
  #         double vector c(theta_old, 1)
  e1 <- parent.frame(1)
  e2 <- parent.frame(2)
  
  tune_constant_2 <- e2$posterior_all[e2$chain_idx, theta_n + 2]
  if (e2$chain_idx < 200)
    cov_estimate <- e2$tune_constant_1 * diag(theta_n)
  else
    cov_estimate <- tune_constant_2 * 
                    cov(e2$posterior_all[e2$iter_idx, e2$theta_idx])
  theta_old <- e2$posterior_all[e2$chain_idx, e2$theta_idx]
  theta_new <- rmvnorm(1, theta_old, cov_estimate)
  theta_reflected <- reflect(theta_new)
  check <- check_support(theta_reflected)
  if (check)
    c(theta_reflected, 1)
  else
    c(theta_old, 0)
}