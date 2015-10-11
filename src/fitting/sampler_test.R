
library("mvtnorm")

sample_am <- function(x, chain_idx) {
  if (chain_idx < 200)
    cov_estimate <- .1 ^ 2 / ncol(x) * diag(ncol(x))
  else
    cov_estimate <- 2.38 ^ 2 / ncol(x) * cov(x)
  rmvnorm(1, x[chain_idx, -21], cov_estimate)
}

set.seed(394398)
mu <- rnorm(20)
sig <- matrix(0, 20, 20)
covs <- runif(190)
sig[upper.tri(sig, FALSE)] <- covs
sig[lower.tri(sig, FALSE)] <- t(sig)[lower.tri(sig, FALSE)]
diag(sig) <- 1:20
isSymmetric(sig)

x <- if (F) matrix(0, nrow = 5e4, ncol = 21) else rbind(x, matrix(0, nrow = 1e5, ncol = 21))

x[1, 21] <- dmvnorm(x[1, 1:20], mu, sig)

for (i in c(1.5e5 + 1):nrow(x)) {
  chain_idx <- i - 1
  old <- x[chain_idx, -21]
  new <- sample_am(x[seq_len(chain_idx), -21, drop = F], chain_idx)
  old_dens <- x[chain_idx, 21]
  new_dens <- dmvnorm(new, mu, sig)
  ratio <- min(new_dens / old_dens, 1)
  if (runif(1) < ratio) 
    x[i, ] <- c(new, new_dens)
  else 
    x[i, ] <- c(old, old_dens)
  
}

save(x, file = "~/Dropbox/Slava/Masters/src/fitting/sampler_test.RData")
load("~/Dropbox/Slava/Masters/src/fitting/sampler_test.RData")

plot.ts(x[, 1:10])
plot.ts(x[, 11:20])
plot.ts(log(x[, 21]))

burn <- -c(1:1.5e5)
round(colMeans(x[burn, ]), 3)
round(cov(x[burn, ]), 3)[1, ]
round(diag(cov(x[burn, ])), 3)





