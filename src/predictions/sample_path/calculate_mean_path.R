
#load(file = "results/sample_path/paths.RData")
library(package = "magrittr")
library(package = "doParallel")

calc_mean_path <- function(path_sample) {
  # Calculates mean sample path
  # Takes a list of numeric vectors and returns a numeric vector
  min_length <- path_sample %>% 
                sapply(X = ., FUN = length, simplify = TRUE) %>% 
                min
  mean_path <- path_sample %>% lapply(X = ., FUN = rev) %>% 
               sapply(X = ., FUN = . %>% extract(1:min_length),
                      simplify = TRUE) %>%
               rowMeans %>% rev
  return(mean_path)
}

calc_model_paths <- function(paths, cores = 1) {
  # Calculates mean sample paths for all simulated paths
  # Takes a list of lists of numeric vectors, scalar integer cores,
  # and returns a list of numeric vectors
  registerDoParallel(cores = cores)
  mean_path_model <- foreach(paths_model = iter(obj = paths, 
                                                by = "cell")) %dopar% {
                     lapply(X = paths_model, FUN = calc_mean_path)
  }
  attr(mean_path_model, "names") <- c("independent", "normal", "t")
  return(mean_path_model)
}