
# Contaminant removal functions for use in the preprocessing script

remove_slow <- function(rr_processed, center_meas, n_sd) {
  # Purpose: removes slow contamints using a number of standard deviations
  # from a measure of central tendency as a cut off
  # Input: numeric data.frame rr_data, function center_meas, 
  # scalar integer n_sd
  # Output: numeric data.frame cleaned
  
  cleaned <- group_by(rr_processed, subj, instr, prop) %>%
    do({
      upper_bound <- center_meas(.$rt) + n_sd * sd(.$rt)
      filter(., rt <= upper_bound)
    })
  cat("removed: ", nrow(cleaned))
  return(cleaned)
}

bin_bound <- function(rt, acc, window_size = 50, p = 0.53) {
  # Purpose: finds lower bound for removing fast reaction times
  # Input: numeric vector rt, integer vector acc, 
  # numeric scalars window_size, p
  # Ouput: numeric scalar bin_breaks[bin + 1]
  
  bin_breaks <- seq(from = 0, to = 1000, by = window_size)
  bin_tags <- findInterval(x = rt, vec = bin_breaks)
  bin_n <- ceiling(1000 / window_size)
  for (bin in seq_len(bin_n)) {
    bin_rts <- which(bin_tags == bin)
    if (length(bin_rts) == 0) 
      next
    prob <- mean(acc[bin_rts])
    if (prob <= p)
      next
    else
      return(bin_breaks[bin])
  }
  return(0)
}

remove_fast <- function(rr_processed, window_size, p) {
  # Purpose: removes fast outliers using accuracy over bins
  # Input: numeric scalars window_size, p
  # Output: numeric data.frame cleaned
  
  cleaned <- do(rr_processed, {
    lower_bound <- bin_bound(.$rt, .$acc, window_size, p)
    filter(., rt >= lower_bound)
  })
  cat("removed: ", nrow(cleaned))
  return(cleaned)
}
