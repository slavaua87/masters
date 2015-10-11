

res <- function(partial_res) {
  ind <- c(1, 4:11, 14:20)
  partial_res[, ind] <- exp(partial_res[, ind])
  partial_trans <- partial_res
  partial_res[, 7] <- partial_trans[, 7] / (partial_trans[, 7] + partial_trans[, 8])
  partial_res[, 17] <- partial_trans[, 17] / (partial_trans[, 17] + partial_trans[, 18])
  partial_res[, 8] <- sqrt(partial_trans[, 7] * partial_trans[, 8] / 
                      (partial_trans[, 7] + partial_trans[, 8]) ^ 2 / 
                      (partial_trans[, 7] + partial_trans[, 8] + 1))
  partial_res[, 18] <- sqrt(partial_trans[, 17] * partial_trans[, 18] / 
                       (partial_trans[, 17] + partial_trans[, 18]) ^ 2 / 
                       (partial_trans[, 17] + partial_trans[, 18] + 1))
  partial_res[, 9] <- partial_trans[, 9] * partial_trans[, 10]
  partial_res[, 10] <- sqrt(partial_trans[, 9]) * partial_trans[, 10]
  partial_res[, 19] <- partial_trans[, 19] * partial_trans[, 20]
  partial_res[, 20] <- sqrt(partial_trans[, 19]) * partial_trans[, 20]
  partial_res  
}

# fin <- data_frame(cota = 13.50,
#                   genfee = 184,
#                   health = 1277,
#                   instr = 5780,
#                   recr = 123,
#                   act = 37.50,
#                   union = 74.40)
# 
# sum(fin)

setwd("~/mount/giverny/Masters/")
library("coda")
load("results/fitting/jf-posterior-chains-norm-fit.RData")
#partial_res[1:18060, ] <- y
#partial_res[-c(1:18060), ] <- 0
partial_res <- partial_res[!partial_res[, 1] == 0, ]
partial_res <- res(partial_res)

#save(partial_res, file = "results/fitting/kr-posterior-chains-norm-fit.RData")

dim(partial_res)
theta_idx <- 23
burn <- if (F) seq_len(nrow(partial_res)) else -(1:100000)
chains_all <- mcmc(partial_res[burn, ])

1 - rejectionRate(chains_all)

# Can I obtain correlations from the correlations between two consequitive
# trials? No knowledge of the starting point, only stimulus for each trial.

lattice::xyplot(chains_all[, 1:6])
lattice::xyplot(chains_all[, 7:12])
lattice::xyplot(chains_all[, 13:18])
lattice::xyplot(chains_all[, 19:23])

lattice::xyplot(chains_all[, 24]) 
lattice::xyplot(chains_all[, 25])

acfplot(chains_all, lag.max = 200)
lattice::levelplot(chains_all)
lattice::levelplot(chains_all[[2]])
lattice::levelplot(chains_all[[3]])

lattice::densityplot(chains_all[, 1:6], plot.points = FALSE)
lattice::densityplot(chains_all[, 7:12], plot.points = FALSE)
lattice::densityplot(chains_all[, 13:18], plot.points = FALSE)
lattice::densityplot(chains_all[, 19:23], plot.points = FALSE)

lattice::densityplot(chains_all[, c(7, 8, 17, 18)], plot.points = FALSE)
lattice::densityplot(chains_all[, c(9, 10, 19, 20)], plot.points = FALSE)
lattice::densityplot(chains_all[, c(6, 16)], plot.points = FALSE)
lattice::densityplot(chains_all[, c(1, 11)], plot.points = FALSE)
lattice::densityplot(chains_all[, c(2, 3, 12, 13)], plot.points = FALSE)
lattice::densityplot(chains_all[, c(4, 5, 14, 15)], plot.points = FALSE)
lattice::densityplot(chains_all[, c(21, 22, 23)], plot.points = FALSE)


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


library("ggmcmc")
chains_ggs <- ggs(chains_all[, 1:23]) 
ggs_compare_partial(chains_ggs) + theme_solarized_2(light = T)
ggmcmc(chains_ggs, plot = c("density", "running", "caterpillar"))
ci(chains_ggs)
ggs_running(chains_ggs) + theme_solarized_2()
ggs_caterpillar(chains_ggs) + theme_solarized_2()
ggs_density(chains_ggs) + theme_solarized_2()
ggs_geweke(chains_ggs)



pdf("results/fitting/traceplots.pdf")
ggs_traceplot(chains_ggs) + theme_solarized_2()
dev.off()

pdf("results/fitting/geweke.pdf")
ggs_geweke(chains_ggs) + theme_solarized_2()
dev.off()

pdf("results/fitting/crosscor.pdf")
ggs_crosscorrelation(chains_ggs) + theme_solarized_2()
dev.off()

pdf("results/fitting/autocor.pdf")
ggs_autocorrelation(chains_ggs, nLags = 100) + theme_solarized_2()
dev.off()

pdf("results/fitting/density.pdf")
ggs_density(chains_ggs) + theme_solarized_2()
dev.off()











