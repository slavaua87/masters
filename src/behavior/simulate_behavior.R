


sample_par <- function(ind_param, tasks, group_n, cond_n, sim_size, task_n, 
                       outer, sigma, model) {
  # Purpose: samples behavioral data in parallel
  # Inputs: double matrix ind_param, integer vector tasks,
  #         integer scalars group_n, cond_n, sim_size,
  #         task_n, outer, double scalar sigma, string model
  # Output: double matrix results
  
  results <- matrix(0, ncol = 3, nrow = outer * 10)
  for (i in seq_len(group_n)) {
    cat(paste("group simulating", "\n"), 
        file = "results/behavior/progress.txt",
        append = TRUE)
    timer <- proc.time()
    results[((i - 1) * cond_n * 10 + 1):(i * cond_n * 10), ] <- foreach(
      params = iter(obj = nextElem(ind_param)$value, by = "row"), 
      .multicombine = TRUE, .combine = "combine_quant",
      .maxcombine = cond_n) %:% 
      foreach(tasks, seeds = rng[(params[[11]] - 1) * task_n + tasks],
              .combine = "rbind", .multicombine = TRUE,
              .maxcombine = sim_size, .inorder = FALSE) %dopar% {
        rngtools::setRNG(seeds)
        cat(paste("simulating", "\n"), 
            file = "results/behavior/progress.txt",
            append = TRUE)
        trial_param <- smpl_param(params = params,
                                  smpl_size = sim_size,
                                  model = model)
        behav_smpl <- smpl_rts(n = 1,
                               alpha = trial_param$alpha,
                               tau = trial_param$t_nd,
                               beta = trial_param$beta,
                               delta = trial_param$delta, 
                               sigma = sigma)
        remove(trial_param)
        behav_smpl
    }
    results
}

simul_behavior <- function(model, smpl_size, sim_size, seed = 1771363045,
                           sigma = .1, cores = 1, group_n = 1) {
  # Purpose: simulates reaction times and choices in parallel for
  #          different copula models 
  # Input: string model, integer scalars smpl_size, seed, 
  #        double scalar sigma, integer scalars cores, group_n 
  # Ouput: numeric data_frame
  
  name <- paste0("results/behavior/progress.txt")
  writeLines(text = "", con = name)
  
  cond <- seq_len(108)
  task_n <- smpl_size / sim_size
  tasks <- seq_len(task_n)
  outer <- length(cond)
  inner <- task_n
  cond_n <- outer / group_n
  rng <- RNGseq(outer * inner, seed)
  group_tags <- seq_len(group_n) %>% rep(each = cond_n)
  ind_param <- cbind(combine_param(nu = nu, wiener = wiener,
                                   rho = rho, omega = omega),
                     cond) %>% isplit(group_tags)

  cat(paste("initiating", "\n"), 
      file = "results/behavior/progress.txt",
      append = TRUE)
  set.seed(seed = seed)
  registerDoParallel(cores = cores)
  results <- sample_par(ind_param, tasks, group_n, cond_n, sim_size,
                        task_n, outer, sigma, model)
  timer <- proc.time() - timer
  cat(paste("simulated", timer["elapsed"], "\n"), 
      file = "results/behavior/progress.txt",
      append = TRUE)
  }
  results <- as.data.frame(results)
  colnames(results) <- c("quant", "prob", "quant_rank")
  results %<>%
    mutate(instruction = c("accuracy", "speed")[rep(c(1, 2), each = 60)] %>%
             rep(times = 9), 
           condition = rep(cond, each = 10),
           model = model)
  results
}




