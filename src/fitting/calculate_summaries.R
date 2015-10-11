

# loads samples
setwd("~/mount/giverny/Masters/")
library("coda")
load("results/fitting/posterior-chains-norm-fit.RData")

setwd("~/Dropbox/Slava/Masters/")
source("src/fitting/update_chains.R")
library("coda")
library("integral")
load("results/fitting/nh-normal-samples.RData")


# cleans 0 and burn-in
partial_res <- partial_res[!partial_res[, 1] == 0, ]
burn <- if (F) seq_len(nrow(partial_res)) else -(1:8e3)
results <- partial_res[burn, ]

# load data
train_data <- read.table("data/train_data.txt") %>% 
  filter(subj == "nh") %>%  # nh may be better
  dplyr::select(rt, resp, prop, instr) %>%
  mutate(rt = rt / 1000, prop = prop / 32) 
# levels(train_data$resp) <- c("upper", "lower")
# train_data %<>% as_data_frame

# model comparison function
calc_ic <- function(results, loglikelihoods) {
  D_bar <- -2 * mean(loglikelihoods)
  theta_bar <- summary(as.mcmc(results))$statistics[, "Mean"]
  dat <- as.matrix(combine_data(as.data.frame(train_data), theta_bar))
  dens <- integral$calc_likelihood_cpp(dat, model, thread_n,
                                       chunk_n, tol, maxEvals)
  D_hat <- -2 * sum(log(dens))
  pD <- D_bar - D_hat
  pV <- var(-2 * loglikelihoods) / 2
  list(dens = dens, DIC1 = pD + D_bar, DIC2 = pV + D_bar,
       BPIC1 = 2 * pD + D_bar, BPIC2 = 2 * pV + D_bar,
       pD = pD, pV = pV, Dbar = D_bar, Dhat = D_hat)
}

calc_waic <- function(results) {
  iter_n <- nrow(results)
  obs_n <- nrow(dat)
  dens <- matrix(0, nrow = obs_n, ncol = iter_n)
  for (i in seq_len(iter_n) {
    dens[, i] <- integral$calc_likelihood_cpp(dat, model, thread_n,
                                              chunk_n, tol, maxEvals)
  }
  lpd <- sum(log(rowMeans(ic_norm$dens)))
  pwaic <- sum(apply(ic_norm$dens, 1, var))
  waic <- -2 * (lpd - pwaic)
  list(waic = waic, lpd = lpd, pwaic = pwaic, dens = dens)
}

# model probabilities function
calc_probs <- function(cors, comp1, comp2, comp3) {
  fn1 <- match.fun(comp1)
  fn2 <- match.fun(comp2)
  fn3 <- match.fun(comp3)
  
  test <- cbind(fn1(cors[, 1], 0),
                fn2(cors[, 2], 0),
                fn3(cors[, 3], 0))
  probs <- colMeans(test)
  probs
}

model <- c("independent", "normal")[1]
thread_n <- 3
chunk_n <- 1
tol <- 1e-4
maxEvals <- 6e5

ic_norm <- calc_ic(results[, 1:23], results[, ncol(results) - 1])
probs_norm <- calc_probs(results[, 21:23], ">", "<", ">")
quant_norm <- list(ic_norm, probs_norm)
load(file = "~/Dropbox/Slava/Masters/results/fitting/quant_norm.RData")

ic_ind <- calc_ic(results[, 1:20], results[, ncol(results) - 1])


