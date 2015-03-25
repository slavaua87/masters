
library("dplyr")
source("src/fitting/integrate_pdf.R")
source("src/fitting/calculate_weibull.R")

lkj_logkernel <- function(rho, omega) {
  logkernel <- log((1 - sum(rho ^ 2) + 2 * prod(rho)) ^ (omega - 1))
  return(logkernel)
}

combine_data <- function(behav_data, theta) {
  # Purpose: forms a data matrix with observations and parameters
  # Input: numeric data.frame behav_data, numeric vector theta
  # Output: numeric data.frame dat_mat
  
  data_mat <- bind_rows(
    transmute(filter(behav_data, instr == 0), rt = rt, choice = resp,
              alpha = theta[1], 
              nu = weibull(prop, theta[2], theta[3], theta[4], theta[5]),
              eta = theta[15], lambda = theta[6], gamma = theta[16],
              chi = theta[7], phi = theta[17], rho_db = theta[18],
              rho_dt = theta[19], rho_bt = theta[20]),
    transmute(filter(behav_data, instr == 1), rt = rt, choice = resp,
              alpha = theta[8], 
              nu = weibull(prop, theta[9], theta[10], theta[11], theta[12]),
              eta = theta[15], lambda = theta[13], gamma = theta[16],
              chi = theta[14], phi = theta[17], rho_db = theta[18],
              rho_dt = theta[19], rho_bt = theta[20]))
  return(data_mat)
}

joint_logdensity <- function(behav_data, theta, model) {
  
  data_mat <- combine_data(behav_data, theta)
  
  if (model == "independent") {
    prior_logdensity <- log(c(dunif(x = theta[1], min = 0.001, max = 0.786),
                              dunif(x = theta[2], min = 0, max = 1.172),
                              dunif(x = theta[3], min = 0, max = 1.172),
                              dunif(x = theta[4], min = 0, max = 50),
                              dunif(x = theta[5], min = 0, max = 5),
                              dunif(x = theta[6], min = 0, max = 1),
                              dunif(x = theta[7], min = 0, max = 1.884),
                              dunif(x = theta[8], min = 0.001, max = 0.786),
                              dunif(x = theta[9], min = 0, max = 1.172),
                              dunif(x = theta[10], min = 0, max = 1.172),
                              dunif(x = theta[11], min = 0, max = 50),
                              dunif(x = theta[12], min = 0, max = 5),
                              dunif(x = theta[13], min = 0, max = 1),
                              dunif(x = theta[14], min = 0, max = 1.884),
                              dunif(x = theta[15], min = 0, max = 0.658),
                              dunif(x = theta[16], min = 0, max = 0.5),
                              dunif(x = theta[17], min = 0, max = 1.260)))
  }
  if (model == "normal") {
    prior_logdensity <- c(prior_logdensity, 
                          lkj_logkernel(rho = theta[18:20], omega = 1))
  }
  behav_logdensity <- log(
    integrate_density_vec(rt = data_mat$rt, 
                          choice = data_mat$choice,
                          sigma = .1, 
                          alpha = data_mat$alpha, nu = data_mat$nu,
                          eta = data_mat$eta, lambda = data_mat$lambda,
                          gamma = data_mat$gamma, chi = data_mat$chi,
                          phi = data_mat$phi, rho_db = data_mat$rho_db,
                          rho_dt = data_mat$rho_dt, rho_bt = data_mat$rho_bt,
                          model = model, tol = 1e-2))
  bayes_logdensity <- sum(prior_logdensity, behav_logdensity)
  return(bayes_logdensity)
}
  