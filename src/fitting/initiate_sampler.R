

check_beta <- function(lambda, gamma) {
  # Purpose: checks that mean and variance of a beta variable are proper
  # Inputs: double vector lambda, gamma
  # Output: logical vector
  
  check <- gamma ^ 2 < lambda * (1 - lambda)
  check
}

check_chi_phi <- function(chi, phi) {
  # Purpose: control for unreasonable combinations of chi,
  #          phi that drastically slow down integration
  # Inputs: double vectors chi, phi
  # Output: logical vector
  
  check <- chi / phi >= 1
  check
}

trans_theta_sampler <- function(theta) {
  # Purpose: transform parameters into the sampling space
  # Inputs: double vector theta
  # Output: double vector theta

  beta_sh1_1 <- ((1 - theta[7]) / theta[8] ^ 2 - 1 / theta[7]) * theta[7] ^ 2
  beta_sh2_1 <- beta_sh1_1 * (1 / theta[7] - 1)
  beta_sh1_2 <- ((1 - theta[17]) / theta[18] ^ 2 - 1 / theta[17]) * theta[17] ^ 2
  beta_sh2_2 <- beta_sh1_2 * (1 / theta[17] - 1)
  tau_sh_1 <- (theta[9] / theta[10]) ^ 2
  tau_sc_1 <- theta[10] ^ 2 / theta[9]
  tau_sh_2 <- (theta[19] / theta[20]) ^ 2
  tau_sc_2 <- theta[20] ^ 2 / theta[19]
  theta[c(7, 8, 17, 18, 9, 10, 19, 20)] <- c(beta_sh1_1, beta_sh2_1,
                                             beta_sh1_2, beta_sh2_2,
                                             tau_sh_1, tau_sc_1,
                                             tau_sh_2, tau_sc_2)
  theta[c(1, 4:11, 14:20)] <- log(theta[c(1, 4:11, 14:20)])
  theta
}

sample_prior <- function(model) {
  # Purpose: sample initial vector of parameters for a given model
  # Inputs: character scalar model
  # Output: double vector params
  # Notes: order of parameters follows exposition of the model
  
  repeat {
    prior_theta <- c(runif(n = 1, min = 0.001, max = 0.786), # alpha
                     runif(n = 1, min = -.5, max = 0), # nu_lo 
                     runif(n = 1, min = 0, max = .5), # nu_hi 
                     runif(n = 1, min = 0, max = 25), # nu_sh 
                     runif(n = 1, min = 0, max = 5), # nu_sc
                     runif(n = 1, min = 0, max = 0.658), # eta
                     runif(n = 1, min = .35, max = .65), # beta_mu 
                     runif(n = 1, min = 0, max = 0.5), # beta_sd
                     runif(n = 1, min = 0, max = 1), # tnd_mu
                     runif(n = 1, min = 0, max = .5), # tnd_sd  
                     runif(n = 1, min = 0.001, max = 0.786), # alpha
                     runif(n = 1, min = -.5, max = 0), # nu_lo 
                     runif(n = 1, min = 0, max = .5), # nu_hi 
                     runif(n = 1, min = 0, max = 25), # nu_sh 
                     runif(n = 1, min = 0, max = 5), # nu_sc
                     runif(n = 1, min = 0, max = 0.658), # eta
                     runif(n = 1, min = .35, max = .65), # beta_mu
                     runif(n = 1, min = 0, max = 0.5), # beta_sd
                     runif(n = 1, min = 0, max = 1), # tnd_mu
                     runif(n = 1, min = 0, max = .5)) # tnd_sd
    
    if (initial)
      prior_theta <- c(.221, -.352, .329, 4.413, .526,
                       .127, .464, .065, .279, .041,
                       .05, -.565, .511, 5.227, .521, 
                       .127, .464, .065, .279, .041) 
      
    if (all(check_beta(prior_theta[7], prior_theta[8]),
            check_beta(prior_theta[17], prior_theta[18]),
            check_chi_phi(prior_theta[9], prior_theta[10]),
            check_chi_phi(prior_theta[19], prior_theta[20])))
      break
  }
  prior_theta <- trans_theta_sampler(prior_theta)
  if (model == "normal") {
    prior_theta <- c(prior_theta, 
                     as.numeric(clusterGeneration::rcorrmatrix(d = 3,
                                            alphad = 1)[upper.tri(diag(3))]))
    
     if (initial) 
       prior_theta <- c(prior_theta[1:20], 0, 0, 0)
  }
  prior_theta
}

















