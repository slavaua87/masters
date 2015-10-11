
calc_response_mean_path <- function(path_sample, alpha) {
  # Calculates mean sample path for each decision
  # Takes a list of numeric vectors path_sample, numeric vector alpha,
  # and returns a numeric vector mean_path
  
  if (length(path_sample) == 0) return(0)
  
  mean_path <- path_sample %>% sapply(length) %>% max %>%
               as.data.frame %>% 
               do({sapply(X = path_sample, 
                          FUN = function(path) {
                                bound <- findInterval(x = rev(path)[1],
                                                      vec = c(-.1, .01, .1, .3))
                                c(path, rep(x = NA,
                                            times = as.numeric(.) - length(path)))             
               },
               simplify = TRUE) %>% as.data.frame}) %>%
               rowMeans(na.rm = TRUE) 
  return(mean_path)
}

calc_mean_paths <- function(path_sample, alpha) {
  resp_index <- sapply(X = path_sample, 
                             FUN = function(x) {
                               min(x) %>% findInterval(vec = c(-.1, 0), 
                                              rightmost.closed = T) == 1
                      })
  mean_paths <- list(upper = calc_response_mean_path(path_sample[!resp_index],
                                                     alpha),
                     lower = calc_response_mean_path(path_sample[resp_index],
                                                     alpha))
  return(mean_paths)
}


calc_model_paths <- function(paths, alpha, cores = 1) {
  # Calculates mean sample paths for all simulated paths
  # Takes a list of lists of numeric vectors paths, numeric vector alpha,
  # integer scalar cores, and returns a list of numeric vectors mean_path_model
  registerDoParallel(cores = cores)
  mean_path_model <- foreach(paths_model = iter(obj = paths, 
                                                by = "cell")) %dopar% {
                     lapply(X = paths_model, FUN = calc_mean_paths, alpha = alpha)
  }
  attr(mean_path_model, "names") <- c("independent", "normal", "t")
  return(mean_path_model)
}
