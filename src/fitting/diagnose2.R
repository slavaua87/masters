

setwd("~/Dropbox/Slava/Masters/")
library("coda")
load("results/fitting/posterior-chains-norm-fit.RData")

dim(partial_res)
theta_idx <- 21
burn <- -seq_len(1)
chains_all <- mcmc.list(mcmc(partial_res[burn, , 1]),
                        mcmc(partial_res[burn, , 2]),
                        mcmc(partial_res[burn, , 3]))

1 - rejectionRate(chains_all)
gelman.diag(chains_all[, 13:20])
gelman.diag(chains_all[, 1:12])

chains_all <- NULL
for (chain in seq_len(nrow(res))) {
  chains_all <- rbind(chains_all, mcmc(t(res[chain, ,])))
}
chains_all <- mcmc(chains_all)
1 - rejectionRate(chains_all)

lattice::xyplot(chains_all[, 1:6])
lattice::xyplot(chains_all[, 7:12])
lattice::xyplot(chains_all[, 13:17])
lattice::xyplot(chains_all[, 18:20])
lattice::xyplot(chains_all[, 21], ylim = c(-100, 0))



acfplot(chains_all, lag.max = 100)
lattice::levelplot(chains_all[[1]])
lattice::levelplot(chains_all[[2]])
lattice::levelplot(chains_all[[3]])

lattice::densityplot(chains_all[, 1:8])
lattice::densityplot(chains_all[, 9:17])
lattice::densityplot(chains_all[, 18:21])

lattice::qqmath(chains_all)

autocorr.plot(chains_all, lag.max = 100)
autocorr.diag(chains_all, c(5, 10, 25, 35))
crosscorr.plot(chains_all)

cumuplot(chains_all[, 1:12])
cumuplot(chains_all[, 12:20])

heidel.diag(chains_all) # length control
geweke.diag(chains_all) # after burn-in
raftery.diag(chains_all) # pilot chain for burn-in/length control
effectiveSize(chains_all) # information quality / length control / mixing

HPDinterval(chains_all)
summary(chains_all)


library(ggmcmc)
chains_ggs <- ggs(chains_all[, 18:20]) 
ggs_traceplot(chains_ggs, original_burnin = 1000) + theme_solarized_2(light = T)
ggs_compare_partial(chains_ggs) + theme_solarized_2(light = T)
ggmcmc(chains_ggs), plot = c("density", "running", "caterpillar"))
ci(chains_ggs)
ggs_running(chains_ggs)











