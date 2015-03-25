
calc_response_path_stats <- function(path_sample, alpha) {
  # Calculates mean sample path for each decision
  # Takes a list of numeric vectors path_sample, numeric vector alpha,
  # and returns a pair of numeric vectors response_stats
  
  if (length(path_sample) == 0) {
    response_stats <- list(mean_path = 0, path_sd = 0)
    return(response_stats)
  }
  
  path_sample %<>% sapply(length) %>% max %>%
               as.data.frame %>% 
               do({sapply(X = path_sample, 
                          FUN = function(path) {
                                bound <- findInterval(x = rev(path)[1],
                                                      vec = c(-.1, .01, .1, .3))
                                c(path, rep(x = alpha[bound],
                                            times = as.numeric(.) - length(path)))             
               },
               simplify = TRUE) %>% as.data.frame}) 

  mean_path <- path_sample %>% rowMeans
  if (length(path_sample) == 1) {
    response_stats <- list(mean_path = mean_path, path_sd = 0)
    return(response_stats)
  }
  
  path_sd <- path_sample %>% apply(MARGIN = 1, FUN = sd)
  response_stats <- list(mean_path = mean_path, path_sd = path_sd)
  return(response_stats)
}

calc_model_paths_stats <- function(path_sample, alpha) {
  # Splits paths into lower and upper bound hitting paths and calculates their
  # mean and sd functions
  # Takes a list of n numeric vectors path_sample, numeric scalar alpha and 
  # returns a list of 2 numeric vectors for each response paths_stats
  
  resp_index <- sapply(X = path_sample, 
                             FUN = function(x) {
                               min(x) %>% findInterval(vec = c(-.1, 0), 
                                              rightmost.closed = T) == 1
                      })
  paths_stats <- list(upper = calc_response_path_stats(path_sample[!resp_index],
                                                       alpha),
                      lower = calc_response_path_stats(path_sample[resp_index],
                                                       alpha))
  return(paths_stats)
}


calc_paths_stats <- function(paths, alpha, cores = 1) {
  # Calculates mean sample paths for all simulated paths
  # Takes a list of lists of numeric vectors paths, numeric vector alpha,
  # integer scalar cores, and returns a list of numeric vectors model_path_stats
  registerDoParallel(cores = cores)
  all_path_stats <- foreach(paths_model = iter(obj = paths, 
                                                by = "cell")) %dopar% {
                      lapply(X = paths_model, FUN = calc_model_paths_stats, 
                             alpha = alpha)
  }
  attr(all_path_stats, "names") <- c("independent", "normal", "t")
  return(all_path_stats)
}
