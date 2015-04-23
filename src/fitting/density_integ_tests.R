
source("src/predictions/behavior/integrate_cdf.R")
library("lineprof")
library("microbenchmark")
x <- lineprof(code = replicate(1000, calc_logdensity(delta = delta_trans, 
                                                     beta = beta,
                                     t_nd = t_nd,
                                     params = as.numeric(params),
                                     shape1 = shape1, shape2 = shape2,
                                     shape = shape, scale = scale,
                                     sqrt_rho = sqrt_rho,
                                     model = model)))
shine(x)

x <- lineprof(code = R2Cuba::cuhre(ndim = 3, ncomp = 1, 
                                   integrand = calc_density_integrand,
                                   rt, choice, sigma, params,
                                   shape1, shape2, shape, scale, 
                                   sqrt_rho, model,
                                   lower = c(-1, 0, 0), 
                                   upper = c(1, 1, rt),
                                   rel.tol = 1e-2,
                                   flags = list(verbose = 0)))
shine(x)

microbenchmark(calc_logdensity(delta = delta_trans, beta = beta,
                               t_nd = t_nd,
                               params = as.numeric(params),
                               shape1 = shape1, shape2 = shape2,
                               shape = shape, scale = scale,
                               sqrt_rho = sqrt_rho,
                               model = model),
               times = 1000,
               unit = "ms")

microbenchmark(R2Cuba::cuhre(ndim = 3, ncomp = 1, 
                             integrand = calc_density_integrand,
                             rt, choice, sigma, params,
                             shape1, shape2, shape, scale, 
                             sqrt_rho, model,
                             lower = c(-1, 0, 0), 
                             upper = c(1, 1, rt),
                             rel.tol = 1e-1,
                             flags = list(verbose = 0)),
               adaptIntegrate(f = calc_density_integrand,                   
                              lowerLimit = c(-1, 0, 0), 
                                               upperLimit = c(1, 1, rt), 
                                               rt, choice, sigma, params,
                                               shape1, shape2, shape, scale, 
                                               sqrt_rho, model, tol = 1e-1, 
                                               fDim = 1, maxEval = maxEval),
               times = 10,
               unit = "eps")

system.time(expr = integrate_density_vec(seq(.2, 4.5, .001), 2, .1, 
                                     params = as.numeric(params), 
                                     model = "normal", 
                                     tol = 1e-2, maxEval = 0))

microbenchmark(integrate_density(1, 2, .1, 
                                 params = as.numeric(params), 
                                 model = "independent", 
                                 tol = 1e-2, maxEval = 0),
               integrate_density(1, 2, .1, 
                                 params = as.numeric(params), 
                                 model = "normal", 
                                 tol = 1e-2, maxEval = 0),
               integrate_density(1, 2, .1, 
                                 params = as.numeric(params), 
                                 model = "t", 
                                 tol = 1e-2, maxEval = 0),
               times = 10,
               unit = "eps")



source(file = "src/predictions/behavior/wiener_parameters.R")
source(file = "src/predictions/behavior/combine_parameters.R")
source(file = "src/predictions/behavior/simulate_wiener_parameters.R")
source(file = "src/predictions/behavior/simulate_rndwalk_rts.R")

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
sim <- rts(params, smpl_size = 100000, model)
proc.time() - timer

x <- seq(0.3, .9, .01)

timer <- proc.time()
den_norm <- sapply(x, integrate_density, 
                  choice = choice, sigma = .1, params = as.numeric(params),
                  model = model, tol = 1e-3, maxEval = 0)
proc.time() - timer

den_norm <- den_norm["integral", ] %>% unlist

den_sim <- density(x = filter(sim, choice == 1) %>% 
                       select(rt) %>%
                       unlist, adjust = .7, from = 0.3, to = .9)
prob_sim <- filter(sim, choice == 1) %>% summarise(mu = n() / 100000) %>% unlist

ggplot() +
  geom_line(mapping = aes(x = x, y = den_norm, colour = "green"), size = 1) +
  geom_line(mapping = aes(x = den_sim$x, y = den_sim$y * prob_sim)) +
  geom_vline(x = .1) + xlim(c(.1, 1))


















