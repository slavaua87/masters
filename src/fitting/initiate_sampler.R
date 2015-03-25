
# Makes available truncated t and correlation matrix samplers
rcorrmatrix <- clusterGeneration::rcorrmatrix

# Checks conditions for the beta density
check_beta <- function(lambda, gamma) {
  # Purpose: checks that mean and variance of a beta variable are proper
  # Input: numerical vectors lambda, gamma
  # Output: logical vector
  
  gamma ^ 2 < lambda * (1 - lambda)
}

# Check relation of chi to phi
check_chi_phi <- function(chi, phi) {
  # Purpose: control for unreasonable combinations of chi, phi that drastically slow down
  # integration
  # Input: numeric vectors chi, phi
  # Output: logical vector
  
  chi / phi >= 1
}

sample_prior <- function(model) {
  # Purpose: sample initial vector of parameters for a given model
  # Input: character scalar model
  # Output: numerical vector params
  # Notes: order of parameters follows exposition of the model
  
  while (TRUE) {
    prior_theta <- c(runif(n = 1, min = 0.001, max = 0.786),
                     runif(n = 1, min = 0, max = .5), #1.172
                     runif(n = 1, min = 0, max = .5), #1.172
                     runif(n = 1, min = 0, max = 50),
                     runif(n = 1, min = 0, max = 5),
                     runif(n = 1, min = 0, max = 1),
                     runif(n = 1, min = 0, max = 1.884),
                     runif(n = 1, min = 0.001, max = 0.786),
                     runif(n = 1, min = 0, max = .5), #1.172
                     runif(n = 1, min = 0, max = .5), #1.172
                     runif(n = 1, min = 0, max = 50),
                     runif(n = 1, min = 0, max = 5),
                     runif(n = 1, min = 0, max = 1),
                     runif(n = 1, min = 0, max = 1.884),
                     runif(n = 1, min = 0, max = 0.658),
                     runif(n = 1, min = 0, max = 0.5),
                     runif(n = 1, min = 0, max = 1.260))
    if (all(check_beta(c(prior_theta[6], prior_theta[13]), prior_theta[16]),
            check_chi_phi(c(prior_theta[7], prior_theta[14]), prior_theta[17])))
      break
  }
  if (model == "normal") {
    prior_theta <- c(prior_theta, 
                     as.numeric(rcorrmatrix(d = 3,
                                            alphad = 1)[upper.tri(diag(3))]))
  }
  return(prior_theta)
}

initialize_chains <- function(model, chain_n) {
  # Purpose: sample initial parameter vectors for all the chains
  # Input: character scalar model, integer scalar chain_n
  # Output: numeric matrix initial_state
  
  initial_state <- t.default(replicate(chain_n, sample_prior(model)))
  return(initial_state)
}


















