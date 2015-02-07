
library(package = "magrittr")
library(package = "doParallel")

calc_mean_path <- function(path_sample, alpha) {
  # Calculates mean sample path centered on
  # Takes a list of numeric vectors path_sample, numeric vector alpha,
  # and returns a numeric vector mean_path
  mean_path <- path_sample %>% sapply(length) %>% max %>%
    as.data.frame %>% 
    do({sapply(X = path_sample, 
               FUN = function(path) {
                 bound <- findInterval(x = path[.],
                                       vec = c(-.1, .01, .1, .3))
                 c(path, rep(x = alpha[bound], times = . - length(path)))
               },
               simplify = TRUE) %>% as.data.frame}) %>% 
    rowMeans
  return(mean_path)
}

calc_model_paths <- function(paths, alpha, cores = 1) {
  # Calculates mean sample paths for all simulated paths
  # Takes a list of lists of numeric vectors paths, numeric vector alpha,
  # integer scalar cores, and returns a list of numeric vectors mean_path_model
  registerDoParallel(cores = cores)
  mean_path_model <- foreach(paths_model = iter(obj = paths, 
                                                by = "cell")) %dopar% {
                     lapply(X = paths_model, FUN = calc_mean_path, alpha = alpha)
  }
  attr(mean_path_model, "names") <- c("independent", "normal", "t")
  return(mean_path_model)
}
