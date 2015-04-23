

sample_posterior <- function(train_data, theta_n, chain_n,  
                             draw_n, model, cores, 
                             chunk, seed,
                             continue, chain_last) {
  # Purpose: samples from the posterior using parallelized demcmc transitios
  # Input: integer scalars theta_n, init_n, chain_n, draw_n, cores, seed
  # Output: numeric array posterior
  # Note: theta_n + 1 contains likelihood calculations,
  # theta_n + 2 contains acceptance indicator
  
  progress_record <- paste0("results/fitting/progress-log-norm-fit", ".txt")
  writeLines(text = "initialize chains", con = progress_record)
  
  set.seed(seed = seed)
  registerDoParallel()
  
  if (continue == TRUE) {
    history_n <- nrow(history_last)
    init_idx <- seq_len(init_n)
    posterior_all <- rbind(history_last,
                           matrix(data = 0,
                                  nrow = chain_n * draw_n,
                                  ncol = theta_n + 1))
    draw_last <- (history_n - init_n) / chain_n
    save_points <- seq(draw_last + 5, draw_last + draw_n, 5)
    for (draw in 2:draw) {
      timer <- proc.time()
      posterior_all[draw_idx, ] <- update_chains(train_data, 
                                                 posterior_all[seq_len(history_n), ],
                                                 chain_n, theta_n,
                                                 history_n, cores,
                                                 chunk)
      timer <- proc.time() - timer
      cat(c(draw_idx, draw, "draws has been finished", timer["elapsed"]),
          file = progress_record,
          sep = "\n",
          append = TRUE)
      if (any(draw == save_points)) {
        partial_res <- posterior_all[-init_idx, ]
        save(partial_res, 
             file = "results/fitting/posterior-chains-norm-fit.RData")
      }
    }
  }
  if (continue == FALSE) {
    theta_idx <- -(theta_n + 1)
    posterior_all <- array(data = 0, dim = c(draw_n, -theta_idx, chain_n))
    posterior_all[1, theta_idx, ] <- t(initialize_chains(model, chain_n))
    data_init <- apply(posterior_all[1, theta_idx, ], 2,
                       combine_data, behav_data = train_data) %>% rbind_all
    timer <- proc.time()
    posterior_all[1, -theta_idx, ] <- 
      joint_logdensity(data_init, t(posterior_all[1, theta_idx, ]),
                       model, chain_n, cores, chunk)
    timer <- proc.time() - timer
    cat(c("init has been finished", timer["elapsed"]),
        file = progress_record,
        sep = "\n", append = TRUE)
    
    save_points <- seq(5, draw_n, 5)
    for (draw in 2:draw_n) {
      if (draw == 2) {
        posterior_cur <- posterior_all[seq_len(draw - 1), , ] 
        dim(posterior_cur) <- c(1, -theta_idx, chain_n)
      } 
      else  
        posterior_cur <- posterior_all[seq_len(draw - 1), , ]
      
      timer <- proc.time()
      posterior_all[draw, , ] <- update_chains(train_data, 
                                               posterior_cur,
                                               chain_n, theta_n,
                                               cores, chunk)
      timer <- proc.time() - timer
      cat(c(draw, "draws has been finished", timer["elapsed"]),
          file = progress_record,
          sep = "\n",
          append = TRUE)
      if (any(draw == save_points)) {
        partial_res <- posterior_all[seq_len(draw), , ]
        save(partial_res, 
             file = "results/fitting/posterior-chains-norm-fit.RData")
      }
    }
  }
  #file.remove(progress_record)
  return(posterior_all)
}
