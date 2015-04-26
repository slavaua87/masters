

library("RWiener")
library("dplyr")
library("doParallel")
library("magrittr")
library("boot")
source("src/fitting/wiener_parameters.R")
source("src/fitting/combine_parameters.R")

n <- 1e5
m <- 1e3
k <- 50
ind_param <- combine_param(nu, wiener, rho, omega)
pars <- ind_param[1, ]

dat_sim <- function(pars, n, m, k) {
  registerDoParallel(cores = k)
  dat <- foreach(seq_len(m), .combine = "rbind") %dopar% 
    with(pars, rwiener(n, alpha / .1, chi, lambda, nu / .1))
  return(dat)
}

stat_boot <- function(data, indices) {
  data_boot <- slice(data, indices)
  rt <- dplyr::select(data_boot, q) * 1000
  p <- dplyr::select(data_boot, resp) %>% unlist %>% as.numeric %>% subtract(1)
  stat <- c(quantile(rt$q, probs = .1, type = 8), 1 - mean(p))
  return(stat)
}
timer <- proc.time()
dat <- dat_sim(pars, n, m, k)
timer <- proc.time() - timer

#res <- boot(dat, stat_boot, 100, parallel = "multicore", ncpus = 10)
#boot.ci(res, type = c("basic", "norm", "perc"), index = 1)
#boot.ci(res, type = c("basic", "norm", "perc"), index = 2)
