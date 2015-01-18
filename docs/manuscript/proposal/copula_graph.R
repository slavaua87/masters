

# load package
library(copula)
# directory for the plot image
setwd('/home/sn/Dropbox/Slava/Masters/latex_version/')
# Explicit copula
frank.cop <- archmCopula('frank', -15, 2)
frank.n <- rCopula(1000, frank.cop)
# Implicit copula
t.cop <- ellipCopula('t', -.95, 2, 'un')
t.n <- rCopula(1000, t.cop)
# check for kendall's tau
cor(frank.n, method = 'kendall')
cor(t.n, method = 'kendall')


# 2x2 plot 
#x11()
jpeg('sdt_copulas.jpeg')
par(mfrow = c(2, 2))
plot(frank.n[, 1], frank.n[, 2], xlab = expression(u[1]),
     ylab = expression(u[2]), main = 'Frank Copula', cex.lab = 1.3)
plot(t.n[, 1], t.n[, 2], xlab = expression(u[1]),
     ylab = expression(u[2]), main = 't Copula', cex.lab = 1.3)
plot(qt(frank.n[, 1], 3, 1.5), qgamma(frank.n[, 2], 3, 6),
     xlab = expression(mu[s]), ylab = expression(sigma[s]), 
     main = 'Frank-based Model', cex.lab = 1.3)
plot(qt(t.n[, 1], 3, 1.5), qgamma(t.n[, 2], 3, 6), 
     xlab = expression(mu[s]), ylab = expression(sigma[s]), 
     main = 't-based Model', cex.lab = 1.3)
graphics.off()

### Frechet-Hoeffding bounds
#u <- rnorm(1000)
#sig <- 100
#x <- exp(u)
#y <- exp(sig * u)
#cor(x,y)
#cor(x, y, method = 'kendall')#
#cor(x, y, method = 'spearman')








