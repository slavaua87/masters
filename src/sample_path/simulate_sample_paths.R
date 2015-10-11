
par_sim <- function(group_n, ind_param, tasks, task_n,
                    sim_size, smpl_size, model, time_unit) {
  # Purpose: run a parallel simulation to obtain sample paths
  # Inputs: integer scalar group_n, double matrix ind_param, 
  #         integer vector tasks, scalar integer task_n, 
  #         integer scalars sim_size, smpl_size, 
  #         string model, double scalar time_unit
  # Ouputs: list of lists of double vectors
     
  for (i in seq_len(group_n)) {
    timer <- proc.time()
    params <- nextElem(ind_param)$value
    results <- foreach(tasks, 
                       seeds = rng[(params[[11]] - 1) * task_n + tasks],
                       .combine = "c", .multicombine = TRUE, .inorder = FALSE,
                       .maxcombine = task_n) %dopar% {
      rngtools::setRNG(seeds)
      trial_param <- smpl_param(params = params,
                                smpl_size = sim_size,
                                model = model)
      path_sample <- rndwalk_vec(smpl_size = 1, 
                                 delta = trial_param$delta, 
                                 sigma = sigma, trial_param$beta, 
                                 trial_param$alpha,
                                 low_bound = 0, time_unit = time_unit)
      remove(trial_param)
      return(path_sample)
    } %>% 
      calc_model_paths_stats(smpl_size) %>% 
      list %>% 
      append(results, .)
    timer <- proc.time() - timer
    cat(paste("simulated", timer["elapsed"], "\n"), 
        file = "results/sample_path/progress.txt",
        append = TRUE)
  }
  results
}

simul_paths <- function(model, smpl_size, sim_size, seed = 2132326000,
                        sigma = .1, time_unit = 1e-3, cores = 1, group_n = 1) {
  # Purpose: simulates sample paths in parallel for different copula models 
  # Inputs: string model, integer scalars smpl_size, seed, double scalars
  #         sigma, time_unit, integer scalars cores, group_n
  # Outputs: list of lists of numeric vectors results
  
  name <- paste0("results/sample_path/progress.txt")
  writeLines(text = "", con = name)
  
  cond <- seq_len(72)
  ind_param <- combine_param(nu = nu, wiener = wiener,
                             rho = rho, omega = omega) %>% slice(cond)
  alpha <- c(0, rev(unique(ind_param$alpha)))
  
  task_n <- smpl_size / sim_size
  tasks <- seq_len(task_n)
  outer <- nrow(ind_param)
  inner <- length(tasks)
  cond_n <- outer / group_n
  rng <- RNGseq(outer * inner, seed)
  results <- list()
  group_tags <- seq_len(group_n) %>% rep(each = cond_n)
  ind_param <- cbind(ind_param, cond) %>% isplit(group_tags)
  
 cat(paste("initiating", "\n"), 
     file = "results/sample_path/progress.txt",
     append = TRUE)
  registerDoParallel(cores = cores)
  results <- par_sim(group_n, ind_param, tasks, task_n, params,
                     sim_size, smpl_size, model, time_unit)
  results
}

