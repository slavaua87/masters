
lkj_kernel <- function(rho, omega) {
  kernel <- (1 - sum(rho ^ 2) + 2 * prod(rho)) ^ (omega - 1)
  return(kernel)
}

calc_prior <- function(theta, model) {
  prior_density <- dnorm(theta[1:20], 0, 100)
  if (model == "normal") {
    prior_density <- c(prior_density, 
                       lkj_kernel(rho = theta[21:23], omega = 1))
  }
  return(prior_density)
}

calc_likelihood <- function(data_mat, model, thread_n, chunk_n, tol, maxEval) {
  
  
#   rho_db <- data_mat[1, 10]
#   rho_dt <- data_mat[1, 11]
#   rho_bt <- data_mat[1, 12]
#   
#   if (model == "normal") {
#     rho <- matrix(c(1, rho_db, rho_dt,
#                     rho_db, 1, rho_bt,
#                     rho_dt, rho_bt, 1),
#                   byrow = TRUE, ncol = 3)
#     sqrt_rho <- cholesky(x = rho)
#   }
  
  behav_density <- foreach(data_chunk = iter(as.matrix(data_mat),
                                             by = "row",
                                             chunksize = chunk_n), 
                           .combine = "c", .multicombine = TRUE,
                           .options.multicore = list(cores = thread_n),
                           .maxcombine = 1e3, .inorder = FALSE) %dopar% {
    res <- integrate_density_vec(rt = data_chunk[, 1], 
                             choice = data_chunk[, 2],
                             sigma = .1, 
                             alpha = data_chunk[, 3], nu = data_chunk[, 4],
                             eta = data_chunk[, 5], shape1 = data_chunk[, 6],
                             shape2 = data_chunk[, 7], shape = data_chunk[, 8],
                             scale = data_chunk[, 9], rho_db = data_chunk[, 10],
                             rho_dt = data_chunk[, 11], rho_bt = data_chunk[, 12],
                             model = model, tol = tol, maxEval = maxEval)
    return(res)
  }
#  attr(behav_density, "rng") <- NULL
  zero_test <- behav_density <= 0
  if (any(zero_test))
    behav_density[zero_test] <- .Machine$double.xmin
  return(behav_density)
}

calc_joint <- function(prior_density, behav_density, chain_n) {
  
  chain_idx <- rep(seq_len(chain_n), each = length(behav_density))
  behav_density <- data_frame(dens = behav_density, idx = chain_idx)
  prior_logdensity <- log(prior_density) %>% 
    as.matrix %>%
    colSums %>% 
    as.matrix
  behav_logdensity <- group_by(behav_density, idx) %>% 
    summarise(ll = sum(log(dens))) %>%
    transmute(ll) %>% 
    as.matrix 
  bayes_logdensity <- cbind(prior_logdensity, behav_logdensity) %>% rowSums
  return(bayes_logdensity)
}

joint_logdensity <- function(data_mat, theta_prop, model,
                             chain_n, thread_n, chunk_n, tol, maxEval) {
  
  behav_density <- integral$calc_likelihood_cpp(as.matrix(data_mat), model, thread_n,
                                                chunk_n, tol, maxEval)
  if (chain_n == 1) {
    prior_density <- calc_prior(theta_prop, model)
    bayes_logdensity <- sum(log(c(prior_density, behav_density)))
    return(bayes_logdensity)
  }
  else {
    prior_density <- apply(theta_prop, 1, calc_prior, model = model)
    bayes_logdensity <- calc_joint(prior_density, behav_density, chain_n)
    return(bayes_logdensity)
  }
}





