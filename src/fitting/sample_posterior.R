
source("src/fitting/initiate_sampler.R")
source("src/fitting/sample_demcmc.R")
source("src/fitting/calculate_joint_parallel.R")

sample_posterior <- function(train_data, theta_n, chain_n, draw_n, 
                             model, cores = 1, seed = 1800968452) {
  # Purpose: samples from the posterior using parallelized demcmc transitios
  # Input: integer scalars theta_n, chain_n, draw_n, cores, seed
  # Output: numeric array posterior
  
  progress_record <- paste0("results/fitting/progress-log-test", ".txt")
  writeLines(text = "first draw calculation", con = progress_record)
  
  posterior <- array(data = 0, dim = c(chain_n, theta_n + 1, draw_n))
  
  set.seed(seed = seed)
  registerDoParallel()
  posterior[, seq_len(theta_n), 1] <- initialize_chains(model, chain_n)
  posterior[, theta_n + 1, 1] <- joint_parallel(train_data,
                                                posterior[, seq_len(theta_n), 1],
                                                chain_n,
                                                model,
                                                cores)
  cat("continue sampling", file = progress_record,
      sep = "\n", append = TRUE)
  progress <- progress_estimated(n = draw_n - 1)
  
  for (draw in 2:draw_n) {
    posterior[, , draw] <- update_chains(train_data, posterior[, , draw - 1],
                                         chain_n, theta_n, cores)
    cat(capture.output(progress$tick()),
        file = progress_record,
        sep = "\n",
        append = TRUE)
  }
  cat("finished sampling", file = progress_record,
      sep = "\n", append = TRUE)
  #file.remove(progress_record)
  return(posterior)
}

