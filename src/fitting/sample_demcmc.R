
library(package = "doParallel")
library(package = "doRNG")
library(package = "magrittr")
library(package = "dplyr")

source("src/fitting/check_support.R")
source("src/fitting/calculate_joint.R")

transition_demcmc <- function(train_data, posterior, chain_id, 
                              chain_n, theta_n, model) {
  # Purpose: uses transition kernel based on the DEMCMC
  # proposal to update the state of a Markov chain
  # Input: numeric matrix posterior, integer scalar chain_id, 
  # integer scalars chain_n, theta_n, character scalar model
  # Output: numeric vector theta_new or theta_old
  
  theta_ind <- -(theta_n + 1)
  theta_old <- posterior[chain_id, theta_ind]
  while (TRUE) {
    diff_vectors <- posterior[sample(x = seq_len(chain_n)[-chain_id], size = 2),
                              theta_ind]
    theta_new <- theta_old + 
      runif(n = theta_n, min = .5, max = 1) * 
      (diff_vectors[1, ] - diff_vectors[2, ]) +
      runif(n = theta_n, min = -.01, max = .01)
    if (check_support(theta_new, model))
      break
  }
  joint_new <- joint_logdensity(train_data, theta_new, model)
  joint_old <- posterior[chain_id, -theta_ind]
  if (log(runif(n = 1)) < joint_new - joint_old)
    return(c(theta_new, joint_new))
  else
    return(c(theta_old, joint_old))
}
  
update_chains <- function(train_data, posterior, chain_n, theta_n, cores) {
  # Purpose: parallelizes simulation of the DEMCMC transition kernel for n chains
  # Input: numeric matrix posterior, integer scalars chain_n, theta_n
  # Output: numeric matrix posterior_updated
  
  posterior_updated <- foreach(chain_id = iter(seq_len(chain_n)),
                               .combine = "rbind",
                               .options.multicore = list(cores = cores)) %dorng% {
    timer <- proc.time()
    new_state <- transition_demcmc(train_data, posterior, chain_id, 
                                   chain_n, theta_n, model)
    timer <- proc.time() - timer
    cat(paste(chain_id, timer["elapsed"]), 
        file = "results/fitting/progress-log-test.txt",
        sep = "\n",
        append = TRUE)
    new_state
  }
  return(posterior_updated)
}

