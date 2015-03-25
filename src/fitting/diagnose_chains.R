
setwd("~/Dropbox/Slava/Masters/")
library("coda")
load("results/fitting/posterior-chains-test-2015-03-24.RData")

chains_all <- apply(res[, , ], 2, c)
chains_mcmc <- mcmc(chains_all[, -18])

par(mfrow = c(4, 4), mar = c(1, 1, 1, 1))
traceplot(chains_mcmc[, 1:16])
traceplot(chains_mcmc[, 17])

cumuplot(chains_mcmc[, 1:16])
cumuplot(chains_mcmc[, 17])

acf(chains_mcmc[, 1:16])
acf(chains_mcmc[, 17])

heidel.diag(chains_mcmc)
geweke.diag(chains_mcmc)
raftery.diag(chains_mcmc)
effectiveSize(chains_mcmc)

densplot(chains_mcmc[, 1:16])
densplot(chains_mcmc[, 17])

HPDinterval(chains_mcmc[, 1:16])
HPDinterval(chains_mcmc[, 17])




