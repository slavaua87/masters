
calc_response_path_stats <- function(path_sample, alpha) {
  # Calculates mean sample path for each decision
  # Takes a list of numeric vectors path_sample, numeric vector alpha,
  # and returns a pair of numeric vectors response_stats
  
  if (length(path_sample) == 0) {
    response_stats <- list(fast_mean = 0, middle_mean = 0, slow_mean = 0,
                           fast_sd = 0, middle_sd = 0, slow_sd = 0)
    return(response_stats)
  }
  
  
  bounds <- sapply(path_sample, length) %>% quantile(c(.33, .66))
  path_sample %<>% sapply(length) %>% max %>%
               as.data.frame %>% 
               do({sapply(X = path_sample, 
                          FUN = function(path) {
                                vec <- sort(c(0, bounds[1],
                                              bounds[2]))
                                speed <- findInterval(x = length(path),
                                                      vec = vec)
                                alpha <- sort(alpha)
                                last_state <- rev(path)[1]
                                threshold <- findInterval(x = last_state,
                                                          vec = c(-0.3, .001,
                                                                  alpha[2] + .1, 
                                                                  alpha[3] + .1))
                                c(path, 
                                  rep(x = alpha[threshold],
                                      times = as.numeric(.) - length(path)),
                                  speed)             
               },
               simplify = TRUE) %>% as.data.frame}) 

  idx_row <- nrow(path_sample)
  if (length(path_sample) == 1) {
    response_stats <- list(fast_mean = path_sample[-idx_row, 
                                                   path_sample[idx_row, ] == 1] %>% 
                             as.matrix %>% rowMeans(na.rm = TRUE), 
                           middle_mean = path_sample[-idx_row, 
                                                   path_sample[idx_row, ] == 2] %>% 
                             as.matrix %>% rowMeans(na.rm = TRUE) ,
                           slow_mean = path_sample[-idx_row, 
                                                   path_sample[idx_row, ] == 3] %>% 
                             as.matrix %>% rowMeans(na.rm = TRUE),
                           fast_sd = 0,
                           middle_sd = 0,
                           slow_sd = 0)
  response_stats <- lapply(response_stats, function(x) x[!is.nan(x)])
    return(response_stats)
  }
  
  response_stats <- list(fast_mean = path_sample[-idx_row, 
                                                 path_sample[idx_row, ] == 1] %>% 
                           as.matrix %>% rowMeans(na.rm = TRUE) , 
                         middle_mean = path_sample[-idx_row, 
                                                   path_sample[idx_row, ] == 2] %>% 
                           as.matrix %>% rowMeans(na.rm = TRUE),
                         slow_mean = path_sample[-idx_row, 
                                                 path_sample[idx_row, ] == 3] %>% 
                           as.matrix %>% rowMeans(na.rm = TRUE),
                         fast_sd = path_sample[-idx_row, 
                                               path_sample[idx_row, ] == 1] %>% 
                           as.matrix %>% apply(MARGIN = 1, FUN = sd, na.rm = TRUE), 
                         middle_sd = path_sample[-idx_row, 
                                                 path_sample[idx_row, ] == 2] %>% 
                           as.matrix %>% apply(MARGIN = 1, FUN = sd, na.rm = TRUE),
                         slow_sd = path_sample[-idx_row, 
                                               path_sample[idx_row, ] == 3] %>% 
                           as.matrix %>% apply(MARGIN = 1, FUN = sd, na.rm = TRUE))
  response_stats <- lapply(response_stats, function(x) x[!is.na(x)])
  return(response_stats)
}

calc_model_paths_stats <- function(path_sample, alpha) {
  # Splits paths into lower and upper bound hitting paths and calculates their
  # mean and sd functions
  # Takes a list of n numeric vectors path_sample, numeric scalar alpha and 
  # returns a list of 2 numeric vectors for each response paths_stats
  
  resp_index <- sapply(X = path_sample, 
                       FUN = function(x) {
                         min(x, na.rm = TRUE) %>% findInterval(vec = c(-.1, 0), 
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
                     path_stats <- lapply(X = paths_model, 
                                          FUN = calc_model_paths_stats,
                                          alpha = alpha)
  }
  attr(all_path_stats, "names") <- c("independent", "normal", "t")
  return(all_path_stats)
}
