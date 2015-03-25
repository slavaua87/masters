
source("src/fitting/integrate_pdf.R")
source(file = "src/predictions/behavior/wiener_parameters.R")
source(file = "src/predictions/behavior/combine_parameters.R")

params <- combine_param(nu = nu, wiener = wiener,
                           rho = rho, omega = omega)[5, ]

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

microbenchmark(integrate_density(1, 2, .1, params$alpha, params$nu, params$eta,
                                 params$lambda, params$gamma, params$chi, params$phi,
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


params <- combine_param(nu = nu, wiener = wiener,
                        rho = rho, omega = omega)[5, ]

system.time(integrate_density(1, 2, .1, params$alpha, params$nu, params$eta,
                              params$lambda, params$gamma, params$chi, params$phi,
                              model = "independent", 
                              tol = 1e-2, maxEval = 0))

posterior[1, -18, 1]
params <- c(.5, .2, posterior[2, 15, 1], posterior[2, 6, 1],
            posterior[2, 16, 1], posterior[2, 7, 1], posterior[2, 17, 1])

params <- c(4.6, 1,
            0.0272078,
            -0.2877625, 0.4686921,
            0.09243406, 0.1685477,
            0.4990823, 0.02194619)

system.time(integrate_density(params[1], params[2], .1, params[3],
                              params[4], params[5],
                              params[6], params[7], params[8], params[9],
                              model = "independent", 
                              tol = 1e-2, maxEval = 0))

lamda = 0.370787; gamma = 0.02411874
chi = 0.02279598; phi = 0.6375755

shape1 <- ((1 - lambda) / gamma ^ 2 -
             1 / lambda) * lambda ^ 2
shape2 <- shape1 * (1 / lambda - 1)

shape <- (chi / phi) ^ 2
scale <- phi ^ 2 / chi

2.0338934      2 0.1431938 -0.06121770 0.4132354 0.370787 0.02411874 0.02279598

calc_logdensity(0.126984126984127, 0.34375, 2.267939e-312,
                -0.06121770, 0.4132354, 
                shape1, shape2, shape, scale, model = model)

theta_test <- c(0.143193788746827,
                0.0612177023901255,
                0.437761852046756,
                6.53183594902776,
                2.00059544961599,
                0.370787013760467,
                0.0227959781472774,
                0.725608116302457,
                0.181540081553821,
                0.268600458192888,
                19.6119612824863,
                1.66588656097463,
                0.212964985404184,
                1.45613928045427,
                0.413235438168435,
                0.0241187358440099,
                0.637575519885987
)
0.251249792438284
0.197098548797642
0.318353127596461
8.77431683063136
2.37340566751904
0.465367204333156
1.30580501542129
0.283142267008866
0.380229383918215
0.211334106012348
31.6328160402147
0.401884614547948
0.396207065686988
0.391335982064184
0.350615327016505
0.107555062687134
0.627046403937985

dat_mat <- combine_data(train_data, theta_test)
joint_logdensity(train_data, theta_test, model)


1, 2, _, 4, 5, 6, 7, 8, _, 10

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
den_ind <- sapply(x, integrate_density, 
                  choice = choice, sigma = .1, params = as.numeric(params),
                  model = model, tol = 1e-3, maxEval = 0)
proc.time() - timer

den_ind <- den_ind["integral", ] %>% unlist

den_sim <- density(x = filter(sim, choice == 1) %>% 
                       select(rt) %>%
                       unlist, adjust = .7, from = 0.3, to = .9)
prob_sim <- filter(sim, choice == 1) %>% summarise(mu = n() / 100000) %>% unlist

ggplot() +
  geom_line(mapping = aes(x = x, y = den_ind, colour = "green"), size = 1) +
  geom_line(mapping = aes(x = den_sim$x, y = den_sim$y * prob_sim)) +
  geom_vline(x = .1) + xlim(c(.1, 1))


















