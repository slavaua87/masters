
source("src/predictions/behavior/integrate_cdf.R")
x <- lineprof(code = integrate_density(.1, 2, .1, 
                                       params = as.numeric(params), 
                                       lower = 0, upper = Inf,
                                       model = "normal", 
                                       tol = 1e-2, maxEval = 0))
shine(x)


integrate_density(.1, 2, .1, params = params, lower = .1, upper = .2,
                  model = "normal", tol = 1e-3, maxEval = 0)


x <- seq(.2, 4.5, .1)


den_ind <- sapply(x, integrate_density, 
                  choice = choice, sigma = .1, params = params,
                  lower = .08, upper = .38,
                  model = model, tol = 1e-3, maxEval = 0)["integral", ] %>%
  unlist
den_nor <- sapply(seq(.05, 5, .1), integrate_density, 
                  choice = choice, sigma = .1, params = params,  
                  model = "normal", tol = 1e-3, maxEval = 0)["integral", ] %>%
  unlist

den_t <- sapply(seq(.05, 5, .1), integrate_density, 
                choice = choice, sigma = .1, params = params,  
                model = "t", tol = 1e-3, maxEval = 0) %>% unlist

den_simp <- dwiener(x, alpha = params[1, "alpha"] / .1, 
                    tau = params[1, "chi"], 
                    beta = params[1, "lambda"], 
                    delta = params[1, "nu"] / .1,
                    resp = rep("upper", length(x)),
                    give_log = F)


ggplot(mapping = aes(x = x, y = den_ind)) + geom_line() +
  geom_line(mapping = aes(x = x, y = den_simp, colour = "blue")) +
  geom_line(mapping = aes(x = x, y = den_nor, colour = "red")) +
  geom_line(mapping = aes(x = x, y = den_t, colour = "green"))


source(file = "src/predictions/behavior/wiener_parameters.R")
source(file = "src/predictions/behavior/combine_parameters.R")
source(file = "src/predictions/behavior/simulate_wiener_parameters.R")
source(file = "src/predictions/behavior/simulate_rndwalk_rts.R")

params2 <- combine_param(nu = nu, wiener = wiener,
                        rho = rho, omega = omega)[5, ]

rts <- function(params, smpl_size, model, sigma = .1) {
  trial_param <- smpl_param(params = params,
                            smpl_size = smpl_size,
                            model = model)
  behav_smpl <- smpl_rts(n = 1,
                         alpha = trial_param$alpha,
                         tau = trial_param$t_nd,
                         beta = trial_param$beta,
                         delta = trial_param$delta, 
                         sigma = sigma)
  return(behav_smpl)
}


model <- "t"

timer <- proc.time()
sim <- rts(params, 100000, model)
proc.time() - timer

x <- seq(0.1, 5, .05)

timer <- proc.time()
den_ind <- sapply(x, integrate_density, 
                  choice = choice, sigma = .1, params = as.numeric(params),
                  lower = 0, upper = Inf,
                  model = model, tol = 1e-2, maxEval = 0)
proc.time() - timer

den_ind <- den_ind["integral", ] %>% unlist

den_sim <- density(x = filter(sim, choice == 1) %>% 
                       select(rt) %>%
                       unlist, adjust = .8, from = 0.1, to = 5)
prob_sim <- filter(sim, choice == 1) %>% summarise(mu = n() / 100000) %>% unlist

ggplot() +
  geom_line(mapping = aes(x = x, y = den_ind, colour = "green"),size = 1) +
  geom_line(mapping = aes(x = den_sim$x, y = den_sim$y * prob_sim)) +
  geom_vline(x = .1) + xlim(c(.1, 5))



den_tnd <- dtmvnorm.marginal(x, n = 1, mean = params[1, "chi"],
                             sigma = params[1, "phi"] ^ 2,
                             lower = 0, upper = Inf)
ggplot() + geom_line(mapping = aes(x = x, y = den_tnd))

  geom_line(mapping = aes(x = x, y = den_simp, colour = "blue")) +
  geom_density(mapping = aes(x = filter(sim, choice == 2) %>% select(rt)))


















obs <- data.frame(x1 = seq(0, 1, .1), x2 = 2 * seq(0, 1, .1), y = rnorm(11))

root <- function(q, prob, alpha, tau, beta, delta, resp, eps = 1-3) {
  fn <- lower_CDF(q, a = alpha, v = delta, s = .1,
                  z = beta * alpha, t_nd = tau, eps) - prob
  return(fn)
}

root <- function(q, prob, alpha, tau, beta, delta, resp) {
  fn <- pwiener(q, alpha, tau, beta, delta, resp) - prob
  return(fn)
}

q <- 0.5140067
prob <- .9
alpha <- .223
tau <- .5
beta <- .5
delta <- -100
resp <- "upper"

uniroot(f = root, interval = c(0, 3), 
        prob, alpha, tau, beta, delta, resp, 
        extendInt = "upX", check.conv = TRUE, tol = 1e-4, maxiter = 10000)

qwiener(prob, alpha / diffusion, tau, beta, delta / diffusion, resp)

microbenchmark(uniroot(f = root, interval = c(0, 3), 
                       prob, alpha, tau, beta, delta, resp, 
                       extendInt = "upX", check.conv = TRUE, 
                       tol = 1e-4, maxiter = 1000),
               qwiener(prob, alpha, tau, beta, delta, resp),
               times = 100, unit = "eps")

pwiener(q, alpha / diffusion, tau, beta, delta / diffusion, resp)
lower_CDF(q, a = alpha, v = delta, s = .1,
          z = beta * alpha, t_nd = tau, eps)




x1 = matrix(c(1, .3, .5, 
              1, .6, .1,
              1, 1.5, .3), byrow = T, 3)


