# # root directory - all files are sourced based on relative paths
# setwd("~/Masters")
# 
# # load dependencies
# source("src/fitting/load_dependencies.R")
# 
# # load data
# train_data <- read.table("data/train_data.txt") %>% 
#   filter(subj == "jf") %>%
#   dplyr::select(rt, resp, prop, instr) %>%
#   mutate(rt = rt / 1000, prop = prop / 32, resp = resp)
# 
# set.seed(245166)
# theta <- initialize_chains(model = model, chain_n = 1) %>% as.numeric
# theta_n <- length(theta)


# # parameter sample from the prior
# data_mat <- combine_data(train_data, theta)
# save(data_mat, file = "src/fitting/data_mat.rdata")

setwd("~/Masters")

library("integral")
load("src/fitting/data_mat.rdata")
#source("src/fitting/load_dependencies.R")

model <- "normal"

data_mat <- as.matrix(data_mat)
tol <- 1e-3
maxEvals <- 5e5
thread_n <- 60
chunk_n <- 1 #; ceiling(nrow(data_mat) / thread_n)

# registerDoParallel()
# timer_r <- proc.time()
# x_r <- calc_likelihood(data_mat[1:1000, ], model, 
#                        thread_n, chunk_n, tol, maxEvals)
# timer_r <- proc.time() - timer_r
# cat(paste("R code: ", timer_r["elapsed"]))

# 
#system(export OMP_SCHEDULE=, paste0(schedule, ", ", chunk_n)))
#cat(c(schedule, thread_n, chunk_n), sep = "\n")

# 
# n <- 7000
# theta <- partial_res[n, 1:23]
# data_mat <- combine_data(train_data, theta)


timer_cpp <- proc.time()
x_cpp <- integral$calc_likelihood_cpp(as.matrix(data_mat), model, thread_n,
                                      chunk_n, tol, maxEvals)
timer_cpp <- proc.time() - timer_cpp
cat(paste("compiled code: ", timer_cpp["elapsed"]))
# sum(log(c(x_cpp, calc_prior(theta, model))))
# partial_res[n, 24]


