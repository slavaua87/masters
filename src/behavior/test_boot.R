

library("RWiener")
library("dplyr")
library("doParallel")
library("magrittr")
library("boot")
setwd("~/mount/giverny/Masters/")
source("src/fitting/wiener_parameters.R")
source("src/predcombine_parameters.R")
source("src/predictions/behavior/load_dependencies.R")

n <- 1e3
m <- 1e3
k <- 3
model <- "normal"
ind_param <- combine_param(nu, wiener, rho, omega)
params <- ind_param[1, ]

dat_sim <- function(params, n, m, k, model) {
  registerDoParallel(cores = k)
  dat <- foreach(seq_len(m), .combine = "rbind",
                 .multicombine = TRUE, .inorder = FALSE) %dopar% {
    trial_param <- smpl_param(params = params,
                              smpl_size = n,
                              model = model)
    behav_smpl <- smpl_rts(n = 1,
                           alpha = trial_param$alpha,
                           tau = trial_param$t_nd,
                           beta = trial_param$beta,
                           delta = trial_param$delta, 
                           sigma = .1)
  return(behav_smpl)
  }
} 

timer <- proc.time()
dat <- dat_sim(params, n, m, k, model)
timer <- proc.time() - timer
timer

names(dat) <- c("q", "resp")

stat_boot <- function(data, indices) {
  data_boot <- slice(data, indices)
  rt <- dplyr::select(data_boot, q) * 1000
  p <- dplyr::select(data_boot, resp) %>% unlist %>% as.numeric %>% subtract(1)
  stat <- c(quantile(rt$q, probs = .1, type = 8), 1 - mean(p))
  return(stat)
}

res <- boot(dat, stat_boot, 200)
boot.ci(res, type = c("basic", "norm", "perc"), index = 1)
boot.ci(res, type = c("basic", "norm", "perc"), index = 2)


## Accuracy
# rts
# 1e4
(636 - 578) / 4 = 14.5 #1 2.5
(899 - 838) / 4 = 15.3 #1 2.5
(1208 - 1127) / 4 = 20.3 #1.25 3
(1660 - 1465) / 4 = 48.8 #3 7.5
(2678 - 2394) / 4 = 71 #4 10 

# 1e5
(599.5 - 593.6) / 4 = 1.475 #1
(869.6 - 861.6) / 4 = 2 #1.3
(1199 - 1186) / 4 = 3.25 #1.6
(1687 - 1671) / 4 = 4 #3
(2738 - 2704) / 4 = 8.5 #3.5 

# 1e6
(596.2 - 594.6) / 4 = .4 /  0.325 #1
(863.7 - 861.4) / 4 = .575 / 0.475 #1.3
(1188 - 1185) / 4 = .75 / .625 #1.8
(1677 - 1672) / 4 = 1.25 / 1 #3.1
(2722 - 2712) / 4 = 2.5 / 2.25 #6.2

# probability
# 1e6
(0.6244 - 0.6224) / 4 = 5e-4

## Speed
# rts
# 1e6
(506.1 - 505) / 4 = .275
(692.4 - 690.7) / 4 = .425
(918.5 - 915.7) / 4 = .7
(1278 - 1274) / 4 = 1
(2130 - 2120) / 4 = 2.5

# probability
# 1e6
(.5577 - .5557) / 4 = 5e-04







