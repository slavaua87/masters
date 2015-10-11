

extend_paths <- function(path_sample) {
  # Purpose: extends all paths with boundary at which they absorbed
  # Inputs: list of double vectors path_sample
  # Outputs: list of double vectors path_sample
  
  bounds <- sapply(path_sample, length) %>% quantile(c(.35, .7))
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
  path_sample
}


calc_response_path_stats <- function(path_sample) {
  # Purpose: calculates sample path means and standard deviations
  #          for each speed and both decisions
  # Inputs: list of double vectors path_sample
  # Outputs: list of double vectors response_stats
  
  if (length(path_sample) %in% c(0, 1)) {
    stop("Need more than one path in a sample")
  }
  
  alpha <- c(0, .050, .221)
  path_sample <- extend_paths(path_sample)
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
  response_stats
}

calc_model_paths_stats <- function(path_sample, smpl_size) {
  # Purpose: splits paths into lower and upper bound hitting paths and
  #          calculates their means and standard deviations
  # Inputs: list of n double vectors path_sample, integer scalar smpl_size 
  # Outputs: list of lists of 2 double vectors paths_stats
  
  if (length(path_sample) != smpl_size)
    return(path_sample)
  timer <- proc.time()
  resp_index <- sapply(X = path_sample, 
                       FUN = function(x) {
                         min(x, na.rm = TRUE) %>% 
                           findInterval(vec = c(-.1, 0), 
                                        rightmost.closed = T) == 1
                       })
  timer <- proc.time() - timer
  cat(paste("calculated", timer["elapsed"], "\n"), 
      file = "results/sample_path/progress.txt",
      append = TRUE)
  paths_stats <- list(upper = calc_response_path_stats(path_sample[!resp_index]),
                      lower = calc_response_path_stats(path_sample[resp_index]))
  paths_stats
  
  
}

combine_stats <- function(...) { 
  # Purpose: combines path means and standard deviations across conditions
  # Input: lists of double vectors
  # Output: list of lists of numeric vectors stats
  
  paths <- list(...)
  cat(paste("combining", "\n"), 
      file = "results/sample_path/progress.txt",
      append = TRUE)
  smpl_size <- length(paths[[2]])
  timer <- proc.time()
  stats <- mclapply(paths, calc_model_paths_stats,
                    smpl_size = smpl_size, 
                    mc.cores = 1) 
  timer <- proc.time() - timer
  cat(paste("combined", timer["elapsed"], "\n"), 
      file = "results/sample_path/progress.txt",
      append = TRUE)
  if (length(stats[[1]][[1]]) == 6) 
    return(stats)
  else 
    stats <- append(stats[[1]], stats[-1])
  stats
}
