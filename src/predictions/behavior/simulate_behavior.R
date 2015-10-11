

simul_behavior <- function(model, smpl_size, sim_size, seed = 1771363045,
                           sigma = .1, cores = 1, group_n = 1) {
  # Purpose: simulates reaction times and choices in parallel for
  # different copula models 
  # Input: character scalar model, numerical scalars smpl_size,
  # seed, sigma, cores
  # Ouput: numeric data.frame
  
  name <- paste0("results/behavior/progress.txt")
  writeLines(text = "", con = name)
  
  cond <- seq_len(108)
  task_n <- smpl_size / sim_size
  tasks <- seq_len(task_n)
  outer <- length(cond)
  inner <- task_n
  cond_n <- outer / group_n
  rng <- RNGseq(outer * inner, seed)
  results <- matrix(0, ncol = 3, nrow = outer * 10)
  group_tags <- seq_len(group_n) %>% rep(each = cond_n)
  ind_param <- cbind(combine_param(nu = nu, wiener = wiener,
                                   rho = rho, omega = omega),
                     cond) %>% isplit(group_tags)

  
  set.seed(seed = seed)
  registerDoParallel(cores = cores)
  for (i in seq_len(group_n)) {
  results[((i - 1) * cond_n * 10 + 1):(i * cond_n * 10), ] <- foreach(
    params = iter(obj = nextElem(ind_param)$value, by = "row"), 
    .multicombine = TRUE, .combine = "combine_quant",
    .maxcombine = cond_n) %:% 
    foreach(tasks, seeds = rng[(params[[11]] - 1) * task_n + tasks],
            .combine = "rbind", .multicombine = TRUE,
            .maxcombine = sim_size, .inorder = FALSE) %dopar% {
      rngtools::setRNG(seeds)
      print("bam")
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
    return(behav_smpl)
    }
  }
  results <- as.data.frame(results)
  colnames(results) <- c("quant", "prob", "quant_rank")
  results %<>%
    mutate(instruction = c("accuracy", "speed")[rep(c(1, 2), each = 60)] %>%
             rep(times = 9), 
           condition = rep(cond, each = 10),
           model = model)
# file.remove(name)
  return(results)
}




