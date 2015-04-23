

# load package
library("copula")
library("gridExtra")
# directory for the plot image
setwd('/home/sn/Dropbox/Slava/Masters/')
# Explicit copula
frank.cop <- archmCopula('frank', -15, 2)
frank.n <- rCopula(1000, frank.cop)
# Implicit copula
t.cop <- ellipCopula('t', -.95, 2, 'un')
t.n <- rCopula(1000, t.cop)
# check for kendall's tau
cor(frank.n, method = 'kendall')
cor(t.n, method = 'kendall')

dat <- data.frame(u1_f = frank.n[, 1], u2_f = frank.n[, 2],
                  u1_t = t.n[, 1], u2_t = t.n[, 2],
                  mu_f = qt(frank.n[, 1], 3, 1.5), 
                  sig_f = qgamma(frank.n[, 2], 3, 6),
                  mu_t = qt(t.n[, 1], 3, 1.5),
                  sig_t = qgamma(t.n[, 2], 3, 6))

u_f <- ggplot(dat) + geom_point(aes(x = u1_f, y = u2_f), size = 1) +
  xlab(expression(u[1])) + ylab(expression(u[2])) +
  ggtitle('Frank Copula') + theme_solarized_2()
u_t <- ggplot(dat) + geom_point(aes(x = u1_t, y = u2_t), size = 1) +
  xlab(expression(u[1])) + ylab(expression(u[2])) +
  ggtitle('t Copula') + theme_solarized_2()
p_f <- ggplot(dat) + geom_point(aes(x = mu_f, y = sig_f), size = 1) +
  xlab(expression(mu[s])) + ylab(expression(sigma[s])) +
  ggtitle('Frank-based Probability Density') + theme_solarized_2()
p_t <- ggplot(dat) + geom_point(aes(x = mu_t, y = sig_t), size = 1) +
  xlab(expression(mu[s])) + ylab(expression(sigma[s])) +
  ggtitle('t-based Probability Density') + theme_solarized_2()
pdf("results/exploratory/sdt_copula.pdf")
grid.arrange(u_f, u_t, p_f, p_t)
dev.off()








