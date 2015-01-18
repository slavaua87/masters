

# load package
library(copula)
# directory for the plot image
setwd('/home/sn/Dropbox/Slava/Masters/Presentations/')

# Independent copula
unif.n <- matrix(runif(2000, 0, 1), 1000, 2)

# 1x2 plot 
#x11()
jpeg('diff_stand.jpeg', width = 580, height = 250, res = 100)
par(mfrow = c(1, 2), mar = c(4, 4.5, 1.5, 1.5), mgp = c(1.5, .5, 0))
plot(unif.n[, 1], unif.n[, 2], xlab = expression(u[1]),
     ylab = expression(u[2]), main = 'Independent Copula',
     cex.lab = 1, cex.axis = .5, xlim = c(0, 1), ylim = c(0, 1))
plot(qnorm(unif.n[, 1], 0.3, 1), qunif(unif.n[, 2], .2, .8),
     xlab = expression(delta), ylab = expression(xi), 
     main = 'Standard Model', cex.lab = 1, cex.axis = .5,
     xlim = c(-2.5, 2.5), ylim = c(0, 1))
graphics.off()

# Explicit copula
frank.cop <- archmCopula('frank', 15, 2)
frank.n <- rCopula(1000, frank.cop)

# 1x2 plot 
#x11()
jpeg('diff_copulas.jpeg', width = 580, height = 250, res = 100)
par(mfrow = c(1, 2), mar = c(4, 4.5, 1.5, 1.5), mgp = c(1.5, .5, 0))
plot(frank.n[, 1], frank.n[, 2], xlab = expression(u[1]),
     ylab = expression(u[2]), main = 'Frank Copula',
     cex.lab = 1, cex.axis = .5, xlim = c(0, 1), ylim = c(0, 1))
plot(qnorm(frank.n[, 1], 0.3, 1), qunif(frank.n[, 2], .2, .8),
     xlab = expression(delta), ylab = expression(xi), 
     main = 'Frank-based Model', cex.lab = 1, cex.axis = .5,
     xlim = c(-2.5, 2.5), ylim = c(0, 1))
graphics.off()

### Frechet-Hoeffding bounds
#u <- rnorm(1000)
#sig <- 100
#x <- exp(u)
#y <- exp(sig * u)
#cor(x,y)
#cor(x, y, method = 'kendall')#
#cor(x, y, method = 'spearman')








