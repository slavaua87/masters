
### Vectorizes rwiener sampler
library(RWiener)
library(magrittr)

rwiener_num <- function(n, alpha, tau, beta, delta) {
  # Purpose: converts data.frame output of rwiener into a matrix to enable
  # vectorization
  # Input: integer scalar n, numeric scalars alpha, tau, beta, delta
  # Output: numeric matrix rt_r
  rt_r <- rwiener(n, alpha, tau, beta, delta) %>% unlist %>%
    matrix(nrow = n)
  return(rt_r)
}

rwiener_vec <- Vectorize(FUN = rwiener_num, 
                         vectorize.args = c("alpha", "tau", "beta", "delta"))

smpl_rts <- function(n, alpha, tau, beta, delta, sigma) {
  # Purpose: rescales parameters and cleans up the reaction time/choice sample
  # Input: integer scalar n, numeric scalars alpha, tau, beta, delta, sigma
  # Output: numeric data.frame rt_r with row and column names
  rts <- rwiener_vec(n, alpha / sigma, tau, beta, delta / sigma) %>% 
    matrix(nrow = length(alpha), byrow = TRUE) %>% as.data.frame
  colnames(rts) <- c("rt", "choice")
  rownames(rts) <- length(alpha) %>% seq_len
  return(rts)
}





