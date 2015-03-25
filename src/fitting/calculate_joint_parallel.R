
source("src/fitting/calculate_joint.R")

joint_parallel <- function(train_data, posterior, chain_n, theta_n, cores) {
  # Purpose: parallelizes calculation of the joint density
  # Input: numeric matrices train_data, posterior, integer scalars chain_n, theta_n
  # Output: numeric vector joint_vector
  
  joint_vector <- foreach(chain_id = iter(obj = seq_len(chain_n)),
                          .combine = "c",
                          .options.multicore = list(cores = cores)) %dorng% {
  timer <- proc.time()
  joint <- joint_logdensity(train_data, 
                            posterior[chain_id, ],
                            model)
  timer <- proc.time() - timer
  cat(paste(chain_id, timer["elapsed"]), 
      file = "~/Masters/results/fitting/progress-log-test.txt",
      sep = "\n",
      append = TRUE)
  joint
  }
  return(joint_vector)
}
