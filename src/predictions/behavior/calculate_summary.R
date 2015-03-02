

calc_quant_prob <- function(behav_smpl) {
  # Purpose: calculates point estimates of reaction time quantiles
  # and choice probabilities from a sample of draws from the joint density
  # Input: numerical data.frame behav_smpl
  # Output: numerical data.frame behav_sum
  
  probs <- filter(behav_smpl, choice == 1) %>% 
    summarise(., prob = n() / nrow(behav_smpl)) %>% unlist
  probs <- c(probs, 1 - probs)

  behav_sum <- group_by(behav_smpl, choice) %>% do({
                 
                 behavior <- quantile(x = .$rt, 
                                      probs = seq(.1, .9, .2), 
                                      type = 8) %>% 
                 as.data.frame})
  
  choice_n <- attr(behav_sum, "group_sizes")
  behav_sum %<>% ungroup %>% 
    mutate(prob = rep(probs, times = choice_n), 
           quant_rank = rep(x = seq_len(5), times = 2)) %>%
    rename(., quant = .) %>% select(-choice)
  return(behav_sum)
}