
library(package = "doParallel")
library(package = "doRNG")
library(package = "magrittr")
library(package = "dplyr")

source("src/fitting/check_support.R")
source("src/fitting/calculate_joint.R")

transition_demcmc <- function(posterior, chain, chain_n) {
  # Purpose: uses transition kernel based on the DEMCMC
  # proposal to update the state of a Markov chain
  # Input: numeric matrix posterior, integer scalar chain, integer scalar chain_n
  # Output: numeric vector theta_new or theta_old
  
  theta_n <- length(theta_old)
  theta_old <- posterior[chain, ]
  while (TRUE) {
    diff_vectors <- posterior[seq_len(length.out = chain_n)[-chain][
      sample.int(x = chain_n, size = 2)], ]
    theta_new <- theta_old + 
      runif(n = theta_n, min = 0, max = 1) * 
      (diff_vectors[1, ] - diff_vectors[2, ]) +
      runif(theta_n, min = -.01, max = .01)
    if (check_support(theta_new))
      break
  }
  if (log(runif(n = 1)) < joint_density(theta_old) - joint_density(theta_new))
    return(theta_new)
  else
    return(theta_old)
}
  
update_chains <- function(posterior, chain_n) {
  # Purpose: parallelizes simulation of the DEMCMC transition kernel for n chains
  # Input: numeric matrix posterior, integer scalar chain_n
  # Output: numeric matrix posterior_updated
  
  posterior_updated <- foreach(chain = iter(seq_len(chain_n))) %rng% {
    transition_demcmc(posterior, chain)
  } %>% rbind_all
  return(posterior_updated)
}
