
rwiener_num <- function(alpha, tau, beta, delta) {
  # Purpose: converts data.frame output of rwiener into a matrix to enable
  #          vectorization
  # Input: integer scalar n, numeric scalars alpha, tau, beta, delta
  # Output: numeric matrix rt_r
  
  repeat {
    rt_r <- rwiener(1, alpha, tau, beta, delta) %>% 
      unlist %>%
      matrix(nrow = 1)
    if (rt_r[1, 1] < 5)
      break
  }
  if (identical(rt_r[1, 2], 1))
    rt_r[1, 2] <- 2
  rt_r
}

rwiener_vec <- Vectorize(FUN = rwiener_num, 
                         vectorize.args = c("alpha", "tau", "beta", "delta"))

smpl_rts <- function(alpha, tau, beta, delta, sigma) {
  # Purpose: rescales parameters and cleans up the reaction time/choice sample
  # Input: integer scalar n, double scalars alpha, tau,
  #        beta, delta, sigma
  # Output: numeric data_frame rt_r with row and column names
  
  rts <- rwiener_vec(alpha / sigma, tau, beta, delta / sigma) %>% 
    matrix(nrow = length(alpha), byrow = TRUE) %>%
    as.data.frame
  colnames(rts) <- c("rt", "choice")
  rownames(rts) <- seq_along(alpha)
  rts
}


