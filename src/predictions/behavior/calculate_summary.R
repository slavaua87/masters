

calc_quant_prob <- function(behav_smpl) {
  # Purpose: calculates point estimates of reaction time quantiles
  # and choice probabilities from a sample of draws from the joint density
  # Input: numerical data.frame behav_smpl
  # Output: numerical data.frame behav_sum
  
  if (ncol(behav_smpl) == 3)
    return(behav_smpl)

  probs <- filter(behav_smpl, choice == 1) %>% 
    summarise(., prob = n() / nrow(behav_smpl)) %>% unlist
  probs <- c(probs, 1 - probs)

  behav_sum <- group_by(behav_smpl, choice) %>% do({
                 
                 behavior <- quantile(x = .$rt, 
                                      probs = seq(.1, .9, .2), 
                                      type = 8) %>% 
                 as.data.frame})
  
  choice_n <- attr(behav_sum, "group_sizes")
  quant <- behav_sum[, 2] %>% unlist
  quant_rank <- rep(x = seq_len(5), times = 2)
  prob <- rep(probs, times = choice_n)
  behav_sum <- data.frame(quant, prob, quant_rank)
  return(behav_sum)
}

combine_quant <- function(...) {
  old_res <- list(...)
  results <- mclapply(old_res, calc_quant_prob,
                      mc.cores = length(old_res))
  behav_sum <- bind_rows(results) %>% unlist
  cat(paste("combined", "\n"), 
      file = "results/behavior/progress.txt",
      append = TRUE)
  return(behav_sum)
}

