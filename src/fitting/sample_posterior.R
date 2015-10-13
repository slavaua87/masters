

sample_posterior <- function() {
  # Purpose: samples from the posterior using parallelized demcmc transitios
  # Inputs: global variables defined in a script
  # Output: returns double array posterior_all, 
  #          saves double matrix partial_res to disc
  # Note: theta_n + 1 contains likelihood calculations,
  #       theta_n + 2 contains acceptance indicator
  
  progress_record <- paste0("results/fitting/jf-progress-log-norm-fit", ".txt")
  writeLines("initialize chains", progress_record)
  
  if (continue) {
    set.seed(seed)
    theta_idx <- -(theta_n + 1:2)
    partial_n <- nrow(partial_res[!partial_res[, 1] == 0, ])
    posterior_all <- partial_res
    
    save_points <- seq(5, draw_n + 5, 5)
    for (draw in seq(partial_n + 1, draw_n)) {
      chain_idx <- draw - 1
      iter_idx <- seq_len(chain_idx)
      
      timer <- proc.time()
      posterior_all[draw, ] <- update_chains()
      timer <- proc.time() - timer
      cat(paste("draw", draw, "time ", timer["elapsed"]),               
          file = progress_record,
          sep = "\n",
          append = TRUE)
      
      if (draw %in% save_points) {        
        partial_res[seq_len(draw), ] <- posterior_all[seq_len(draw), ]
        save(partial_res, 
             file = "results/fitting/jf-posterior-chains-norm-fit.RData")
      }
    }
    return(posterior_all)
  } else {
  
    set.seed(seed)
  
    tune_constant_1 <- .001 ^ 2 / theta_n
    tune_constant_2 <- 2.38 ^ 2 / theta_n
    theta_idx <- -(theta_n + 1:2)
    
    posterior_all <- matrix(0, draw_n, theta_n + 2)
    partial_res <- matrix(0, draw_n, theta_n + 2)
  
    posterior_all[1, theta_n + 2] <- tune_constant_2
    posterior_all[1, theta_idx] <- sample_prior(model)
    data_mat <- combine_data(posterior_all[1, theta_idx])
      
    timer <- proc.time()
    posterior_all[1, theta_n + 1] <- 
      integral$joint_logdensity_cpp(data_mat,
                                    posterior_all[1, theta_idx],
                                    model, thread_n, chunk_n,
                                    tol, maxEvals)
    timer <- proc.time() - timer
    
    cat(paste("draw 1 time ", timer["elapsed"]),
        file = progress_record,
        sep = "\n", append = TRUE)
    
    save_points <- seq(5, draw_n + 5, 5)
    for (draw in seq(2, draw_n)) {
      chain_idx <- draw - 1
      iter_idx <- seq_len(chain_idx)
      
      timer <- proc.time()
      posterior_all[draw, ] <- update_chains()
      timer <- proc.time() - timer
      cat(paste("draw", draw, "time ", timer["elapsed"]),               
          file = progress_record,
          sep = "\n",
          append = TRUE)
      
      if (draw %in% save_points) {        
        partial_res[seq_len(draw), ] <- posterior_all[seq_len(draw), ]
        save(partial_res, 
             file = "results/fitting/jf-posterior-chains-norm-fit.RData")
      }
    }
  }
  posterior_all
}
