
####
# multivariate integration package
library("R2Cuba")
library("copula")
library("mvtnorm")
library("cubature")
# mixing density
mixing.dens <- function(x, r, m, s, a, b) {
  cop.code <- ellipCopula('normal', param = r,
  dim = 2, dispstr = 'un') 
  dens <- dCopula(c(pnorm(x[1], m, s), pbeta(x[2], a, b)), cop.code) *
  dnorm(x[1], m, s) * dbeta(x[2], a, b)
  return(dens)
}

test.dens <- function(x, r, m, s, a, b, dim1, dim2) {
  cop.code <- ellipCopula('normal', param = r,
                          dim = 2, dispstr = 'un')
  y1 <- dim1[1] + (dim2[1] - dim1[1]) * x[1]
  y2 <- dim1[2] + (dim2[2] - dim1[2]) * x[2]
  dens <- dCopula(c(pnorm(y1, m, s), pnorm(y2, a, b)), cop.code) *
    dnorm(x[1], m, s) * dnorm(x[2], a, b) * prod(dim2 - dim1)
  return(dens)
}

test.norm <- function(y, mu, sig) {
  # change variable
  x <- y / (1 - y ^ 2)
  # calculate jacobian
  jacob <- (y ^ 2 + 1) / (1 - y ^ 2) ^ 2
  #calculate density
  dens <- dnorm(x, mu, sig) * jacob
  return(dens)
}

dnorm(x[2], mu, sig) * 
  prod(jacob)


mu <- 0
sig <- 1
cubature::adaptIntegrate(test.norm, c(-.5), .5, mu = mu, sig = sig)

cubature::adaptIntegrate(test.norm, c(-1, -1), c(1, 0), mu, sig)
                        
cuhre(ndim = 2, ncomp = 1, integrand = test.norm, mu, sig,
      lower = c(-11, -100), upper = c(0, 100), key = 13)
# parameters
r = 0
m <- .5
s <- .5
a <- .5
b <- .5
dim1 <- c(-1, -1)
dim2 <- c(1, 1)
# numerical integration
cuhre(ndim = 2, ncomp = 1, integrand = test.dens, r, m, s, a, b,
      dim1, dim2, lower = c(0, 0), upper = c(1, 1), key = 13)

pmvnorm(lower = dim1, upper = dim2, mean = c(m, a), 
        sigma = matrix(c(s, r, r, b), 2, 2, byrow = TRUE))



fn(c(1, .5), m, s, a, b)

cop.code <- ellipCopula('normal', param = r,
dim = 2, dispstr = 'un', df = 5)


x <- rnorm(100, 5, 5)
y <- 1 * x + 1 + rnorm(100, 0, 1)
cor(x, y, method = c('pearson'))

windows()
plot(x, y)


# BiNormal - Normal model
cop.code <- ellipCopula('normal', param = -.5,
dim = 2, dispstr = 'un', df = 1000)
cop.sam <- rCopula(100000, cop.code)
con.sam <- cbind(qnorm(cop.sam[, 1], 3, 1), qexp(cop.sam[, 2], 3))
data.sam <- rnorm(100000, mean = con.sam[, 1], sd = con.sam[, 2])

# BiNormal - BiNormal model
cop.code <- ellipCopula('normal', param = -.10,
dim = 2, dispstr = 'un', df = 1000)
cop.sam <- rCopula(10000, cop.code)
con.sam <- cbind(qnorm(cop.sam[, 1], 3, 1), qnorm(cop.sam[, 2], 3, 1))
data.sam <- rnorm(10000, mean = con.sam[, 1], sd = 1)
data.sam2 <- rnorm(10000, mean = con.sam[, 2], sd = 1)

# BiNormal - Normal-Bern model
cop.code <- ellipCopula('normal', param = 1,
dim = 2, dispstr = 'un', df = 1000)
cop.sam <- rCopula(10000, cop.code)
con.sam <- cbind(qnorm(cop.sam[, 1], 3, 1), qbeta(cop.sam[, 2], 5, 1))
data.sam <- rnorm(10000, mean = con.sam[, 1], sd = 1)
data.sam2 <- rbinom(n = 10000, size = 1, p = con.sam[, 2])


df <- 5
rho <- .99
# Kendalls tau
2*asin(rho)/pi
cor(cop.sam, method = 'kendall')
cor(con.sam, method = 'kendall')

# Tail dependence
2 * (pt(-sqrt((df + 1) * (1 - rho) / (1 + rho)), df + 1))

scatter.smooth(cop.sam)

# Bi-uni
plot(density(data.sam))
lines(seq(-3, 6, .1), dnorm(seq(-3, 6, .1), mean(data.sam), sd(data.sam))
acf(data.sam)
mean(data.sam)
var(data.sam)
skewness(data.sam)
kurtosis(data.sam) - 3

# Bi-Bi
acf(data.sam)
acf(data.sam2)
ccf(data.sam, data.sam2)

cor(data.sam, data.sam2, method = 'pearson')
scatter.smooth(data.sam, data.sam2)

### microbenchmark of eigenvalue - 27k / sec
library("microbenchmark")
x <- matrix(c(3, 2, 1,
              2, 1, .5,
              1, .5, 1.5), 3, 3, byrow = TRUE)

times <- microbenchmark(eigen(x, symmetric = TRUE, only.values = TRUE), 
                        times = 1000, unit = "eps")


