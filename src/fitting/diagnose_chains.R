
setwd("~/Dropbox/Slava/Masters/")
library("coda")
load("results/fitting/posterior-chains-norm-fit.RData")

partial_res <- res[-(1:6170), ]
theta_idx <- 21
chains_all <- mcmc.list(mcmc(partial_res[seq(1, 2556 * 3, 3), -theta_idx]),
                        mcmc(partial_res[seq(2, 2556 * 3, 3), -theta_idx]),
                        mcmc(partial_res[seq(3, 2556 * 3, 3), -theta_idx]))
1 - rejectionRate(chains_all)
gelman.diag(chains_all)

chains_all <- chains_all[-c(1:500), ] 

chains_all <- NULL
for (chain in seq_len(nrow(res))) {
  chains_all <- rbind(chains_all, mcmc(t(res[chain, ,])))
}
chains_all <- mcmc(chains_all)
1 - rejectionRate(chains_all)

par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
traceplot(chains_all[, 1:12])
traceplot(chains_all[, -(1:12)])

cumuplot(chains_all[, 1:12])
cumuplot(chains_all[, 12:20])

acfplot(chains_all, lag.max = 100)
autocorr.plot(chains_all, lag.max = 100)
autocorr.diag(chains_all, c(5, 10, 25, 35))
crosscorr.plot(chains_all)

heidel.diag(chains_all) # length control
geweke.diag(chains_all) # after burn-in
raftery.diag(chains_all) # pilot chain for burn-in/length control
effectiveSize(chains_all) # information quality / length control / mixing

par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
densplot(chains_all[, 1:12])
densplot(chains_all[, -c(1:12)])

HPDinterval(chains_all)




