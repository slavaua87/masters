
# Random walk simulation of a bounded Wiener process
# Takes process parameters and number of paths to return a list of sample paths

rand_walk = function(trial = 10, v = .01, s = 1, z = .75, a = 1.3,
                     ter = .2, lb = 0, tau = .0001) {
  
  # Container for data
  results = matrix(0,trial,2)
  
  # Step size
  delta = s*sqrt(tau)
  
  # Probability of moving up
  p.up = .5*(1+v*sqrt(tau)/s)
  if(p.up > 1) p.up = 1
  if(p.up < 0) p.up = 0
  
  # Sampling mechanism
  for (i in 1:trial) {
    
    # Initial position and time
    position = z
    time = 0
    
    # Simulation of a random walk
    while(position < a & position > lb) {
      
      # Make a move
      if(rbinom(1,1,p.up) == 1) position = position + delta
      else position = position - delta
      # Update time
      time = time + tau      
    }
    
    # Add data to the container    
    if(position >= a) results[i,] = c(1,time)
    else results[i,] = c(0,time)
    
  }
  results[,2] <- results[,2] + ter
  # Output data
  return(results)
}


# Vectorize the simulation
rw.vectorized <- Vectorize(FUN = rw.basic,
                           vectorize.args = list('v','z','a','ter'))

dat <- list(a = c(3, 5), b = 4)
do.call(merge, list(x = dat, by = 0, all = TRUE)) 

x <- c(rnorm(10), rep(NA, 10))

