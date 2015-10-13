res <- function(partial_res) {
  # Purpose: map sampled parameter from sampling to model space
  # Inputs: double matrix partial_res
  # Output: double matrix partial_res
  partial_res <- cbind(exp(partial_res[, 1:20]), partial_res[, 21:25])
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

# load, transform, and remove burn-in from sampled values
setwd("~/Masters/")
load("results/fitting/posterior-chains-norm-fit.RData")
partial_res <- partial_res[!partial_res[, 1] == 0, ]
partial_res <- res(partial_res)
burn <- if (F) seq_len(nrow(partial_res)) else -(1:24e3)

# make diagnostic plots
library("coda")
library("ggmcmc")
chains_all <- mcmc(partial_res[burn, ])
chains_ggs <- ggs(chains_all[, 21:23])

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
