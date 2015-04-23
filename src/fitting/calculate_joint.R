
lkj_kernel <- function(rho, omega) {
  kernel <- (1 - sum(rho ^ 2) + 2 * prod(rho)) ^ (omega - 1)
  return(kernel)
}

calc_prior <- function(theta, model) {
  prior_density <- c(dunif(x = theta[1], min = 0.001, max = 0.786),
                     dunif(x = theta[2], min = 0, max = 1.172),
                     dunif(x = theta[3], min = 0, max = 1.172),
                     dunif(x = theta[4], min = 0, max = 50),
                     dunif(x = theta[5], min = 0, max = 5),
                     dunif(x = theta[6], min = 0, max = 1),
                     dunif(x = theta[7], min = 0, max = 1.884),
                     dunif(x = theta[8], min = 0.001, max = 0.786),
                     dunif(x = theta[9], min = 0, max = 1.172),
                     dunif(x = theta[10], min = 0, max = 1.172),
                     dunif(x = theta[11], min = 0, max = 50),
                     dunif(x = theta[12], min = 0, max = 5),
                     dunif(x = theta[13], min = 0, max = 1),
                     dunif(x = theta[14], min = 0, max = 1.884),
                     dunif(x = theta[15], min = 0, max = 0.658),
                     dunif(x = theta[16], min = 0, max = 0.5),
                     dunif(x = theta[17], min = 0, max = 1.260))
  if (model == "normal") {
    prior_density <- c(prior_density, 
                       lkj_kernel(rho = theta[18:20], omega = 1))
  }
  return(prior_density)
}


calc_likelihood <- function(data_mat, model, cores, chunk_n) {
  
  behav_density <- foreach(data_chunk = iter(as.matrix(data_mat),
                                             by = "row",
                                             chunksize = chunk_n),
                           .combine = "c",
                           .options.multicore = list(cores = cores)) %dorng% {
    res <- integrate_density_vec(rt = data_chunk[, 1], 
                                choice = data_chunk[, 2],
                                sigma = .1, 
                                alpha = data_chunk[, 3], nu = data_chunk[, 4],
                                eta = data_chunk[, 5], lambda = data_chunk[, 6],
                                gamma = data_chunk[, 7], chi = data_chunk[, 8],
                                phi = data_chunk[, 9], rho_db = data_chunk[, 10],
                                rho_dt = data_chunk[, 11], rho_bt = data_chunk[, 12],
                                model = model, tol = 1e-2, maxEval = 1e4)
    res
  } 
  zero_test <- behav_density <= 0
  if (any(zero_test))
    behav_density[zero_test] <- .Machine$double.xmin
  attr(behav_density, "rng") <- NULL
  return(behav_density)
}

calc_joint <- function(prior_density, behav_density, chain_n) {
  
  chain_idx <- rep(seq_len(chain_n), each = length(behav_density) / chain_n)
  behav_density <- data.frame(dens = behav_density, idx = chain_idx)
  prior_logdensity <- log(prior_density) %>% 
    as.matrix %>%
    colSums %>% 
    as.matrix
  behav_logdensity <- group_by(behav_density, idx) %>% 
    summarise(ll = sum(log(dens))) %>%
    transmute(ll) %>% 
    as.matrix 
  bayes_logdensity <- cbind(prior_logdensity,
                            behav_logdensity) %>% rowSums
  return(bayes_logdensity)
}

joint_logdensity <- function(data_mat, theta_prop, model,
                             chain_n, cores, chunk) {
  if (chain_n == 1)
    prior_density <- calc_prior(theta_prop, model)
  else
    prior_density <- apply(theta_prop, 1, calc_prior, model = model)
  behav_density <- calc_likelihood(data_mat, model, cores, chunk)
  bayes_logdensity <- calc_joint(prior_density, behav_density, chain_n)
  return(bayes_logdensity)
}

  






