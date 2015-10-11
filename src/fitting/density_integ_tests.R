
cd("~/Dropbox/Slava/Masters")
source("src/fitting/load_dependencies.R")
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

load("results/fitting/posterior-chains-norm-fit.RData")
test <- partial_res[1, -24, 3]

# experimental settings
model <- "normal"
source("src/fitting/load_dependencies.R")
smpl_size <- 3
prop <- seq(from = 1, to = 32, by = 7) / 32

# parameter sample from the prior
set.seed(245186)
theta <- initialize_chains(model = model, chain_n = 1) %>% as.numeric
theta_n <- length(theta)

# sample of data from the experiment
train_data <- sample_experiment(theta = theta, model = model,
                                prop = prop, smpl_size = smpl_size)

data_mat <- combine_data(train_data, test)


params <- data_mat[2, -c(1, 2)] %>% as.matrix
params_rt <- params
params_rt[4] <- params[4] / (params[4] + params[5])
params_rt[5] <- params[4] * params[5] / 
  (params[4] + params[5]) ^ 2 /
  (params[4] + params[5] + 1)
params_rt[6] <- params[6] * params[7]
params_rt[7] <- params[6] * params[7] ^ 2
colnames(params_rt) <- names(cbind(combine_param(nu = nu, wiener = wiener,
                                                 rho = rho, omega = omega))[1, -11])


model <- "normal"
timer <- proc.time()
sim <- rts(params_rt, smpl_size = 1, model)
proc.time() - timer

choice = "lower"
x <- seq(0, 3, .05)
params <- as.numeric(params)
timer <- proc.time()
den_norm <- sapply(x, integrate_density, 
                   choice = choice, sigma = 1, alpha = params[1],
                   nu = params[2], eta = params[3], shape1 = params[4],
                   shape2 = params[5],
                   shape = params[6], scale = params[7], rho_db = params[8],
                   rho_dt = params[9], rho_bt = params[10],
                   model = model, tol = 1e-2, maxEval = 1e4)
proc.time() - timer

model <- "normal"
x <- seq(0.2, 3, .01)
choice <- "upper"
params <- cbind(combine_param(nu = nu, wiener = wiener,
                              rho = rho, omega = omega))[6, -11]
timer <- proc.time()
sim <- rts(params, smpl_size = 1e5, model)
proc.time() - timer

shape1 <- as.numeric(((1 - params["lambda"]) / params["gamma"] ^ 2 -
                        1 / params["lambda"]) * params["lambda"] ^ 2)
shape2 <- as.numeric(shape1 * (1 / params["lambda"] - 1))
shape <- (as.numeric(params[1, "chi"]) / as.numeric(params[1, "phi"])) ^ 2
scale <- as.numeric(params[1, "phi"]) ^ 2 / as.numeric(params[1, "chi"])
params <- as.numeric(params)

timer <- proc.time()
den_norm <- sapply(1.09, integrate_density, 
                   choice = choice, sigma = .1, alpha = params[1],
                   nu = params[2], eta = params[3], shape1 = shape1,
                   shape2 = shape2,
                   shape = shape, scale = scale, rho_db = params[8],
                   rho_dt = params[9], rho_bt = params[10],
                   model = model, tol = 1e-2, maxEval = 1e4)
proc.time() - timer

den_sim <- density(x = filter(sim, choice == 1) %>% 
                       dplyr::select(rt) %>%
                       unlist, adjust = .6, from = .2, to = 3)
prob_sim <- filter(sim, choice == 1) %>% summarise(mu = n() / 1e5) %>% unlist

ggplot() +
  geom_line(mapping = aes(x = x, y = den_norm), colour = "green", size = 1) +
  geom_line(mapping = aes(x = den_sim$x, y = den_sim$y * prob_sim)) +
  geom_vline(x = .1) + xlim(c(.1, 3.1))






dwiener(2.402128, 1.5161616528634834111e+01 / .1,
        7.9750192331705271886e-01,
        1.9653784129904643372e-10,
        3.6945669097124578251e+00 / .1, "lower") #%>% print(digits = 16)

test_mat <- bind_rows(data_mat[1, ], data_mat[1, ], data_mat[1, ], data_mat[1, ])
test_mat <- bind_rows(test_mat, test_mat, test_mat, test_mat)
test_mat <- bind_rows(test_mat, test_mat, test_mat)
test_mat[, 1] <- seq(1, 3, length.out = 48)

test_dens <- calc_likelihood(test_mat, model, 3, 1)

plot(test_mat[, 1] %>% unlist,
     test_dens,
     type = "p", pch = 19)


calc_likelihood(data_mat, model, 3, 1)


joint_logdensity(data_mat, test, model, 1, 3, 1)

int_test <- function(x, alpha, tau, beta, delta, resp) {
  #print(x, digits = 16)
  dens <- dwiener(x, alpha, tau, beta, delta, resp)
  if (is.infinite(dens))
    dens <- .Machine$double.xmin
  dens
}

adaptIntegrate(int_test, 0, 10, alpha = 6.76343 / .1,
               tau = 1.372988,
               beta = .4983964,
               delta = -3.866355 / .1, 
               resp = "lower",
               tol = 1e-2, fDim = 1, maxEval = 1e4)

plot(seq(2, 3.8, .01),
dwiener(seq(2, 3.8, .01), 6.76343 / .1, 1.372988,
        .4983964, -3.866355 / .1, rep("lower", 181)), type = "l", xlim = c(2, 3))

