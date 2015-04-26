
library("boot")

params <- data.frame(lambda = .464, gamma = .065)
shape1 <- as.numeric(((1 - params[1, "lambda"]) / params[1, "gamma"] ^ 2 -
                        1 / params[1, "lambda"]) * params[1, "lambda"] ^ 2)
shape2 <- as.numeric(shape1 * (1 / params[1, "lambda"] - 1))

dat <- rbeta(6e3, shape1, shape2)

stat_boot <- function(data, indices) {
  perm <- data[indices]
  stat <- mean(perm)
  return(stat)
}

res <- boot(dat, stat_boot, 1000)
boot.ci(res, type = c("basic", "norm", "perc"))
