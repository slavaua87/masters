
# Makes available truncated t and correlation matrix samplers
rtmvt <- tmvtnorm::rtmvt
rcorrmatrix <- clusterGeneration::rcorrmatrix

# Tests parameter restrictions are satisfied for the beta density
check_lambda <- function(mu, sig) all(sig ^ 2 < mu * (1 - mu))

# Transforms mean, standard deviation into gamma density shape parameters
shape1 <- function(mu, sig) ((1 - mu) / sig ^ 2 - 1 / mu) * mu ^ 2
shape2 <- function(mu, sig) shape1(mu, sig) * (1 / mu - 1)

sample_prior <- function(model) {
  # Purpose: sample initial vector of parameters for a given model
  # Input: character scalar model
  # Output: numerical vector params
  # Notes: order of parameters follows exposition of the model
  
  while (TRUE) {
    group_param <- c(runif(n = 2, min = 0.001, max = 0.786),
                     runif(n = 4, min = 0, max = 1.172),
                     runif(n = 2, min = 0, max = 50),
                     runif(n = 2, min = 0, max = 5),
                     runif(n = 2, min = 0, max = 1),
                     runif(n = 2, min = 0, max = 1.884),
                     rtmvt(n = 5, sigma = as.matrix(25), df = 1,
                           lower = 0, upper = Inf), 
                     runif(n = 1, min = 0, max = .5),
                     rtmvt(n = 1, sigma = as.matrix(25), df = 1,
                           lower = 0, upper = Inf)) 
    if (check_lambda(mu = group_param[11:12], sig =  group_param[20]))
      break
  }
  
  if (model == "normal") {
      omega <- rtmvt(n = 1, sigma = as.matrix(625), df = 1,
                     lower = 0, upper = Inf) # rho
      group_param <- c(group_param, omega)
  } 
  person_param <- 
    as.numeric(cbind(rgamma(n = 3, 
                            shape = (group_param[1] / group_param[15]) ^ 2, 
                            scale = group_param[15] ^ 2 / group_param[1]),  
                     rgamma(n = 3,            
                            shape = (group_param[2] / group_param[15]) ^ 2,             
                            scale = group_param[15] ^ 2 / group_param[2]),         
                     rgamma(n = 3,        
                            shape = (group_param[3] / group_param[16]) ^ 2, 
                            scale = group_param[16] ^ 2 / group_param[3]),
                     rgamma(n = 3,               
                            shape = (group_param[4] / group_param[16]) ^ 2,   
                            scale = group_param[16] ^ 2 / group_param[4]),         
                     rgamma(n = 3,               
                            shape = (group_param[4] / group_param[16]) ^ 2,                
                            scale = group_param[16] ^ 2 / group_param[4]),         
                     rgamma(n = 3,             
                            shape = (group_param[6] / group_param[17]) ^ 2, 
                            scale = group_param[17] ^ 2 / group_param[6]),       
                     rgamma(n = 3,               
                            shape = (group_param[7] /group_param[18]) ^ 2,   
                            scale = group_param[18] ^ 2 / group_param[7]),         
                     rgamma(n = 3,               
                            shape = (group_param[8] / group_param[18]) ^ 2, 
                            scale = group_param[18] ^ 2 / group_param[8]),         
                     rgamma(n = 3,
                            shape = (group_param[9] / group_param[19]) ^ 2,             
                            scale = group_param[19] ^ 2 / group_param[9]),         
                     rgamma(n = 3,               
                            shape = (group_param[10] / group_param[19]) ^ 2,   
                            scale = group_param[19] ^ 2 / group_param[10]),              
                     rbeta(n = 3,               
                           shape1 = shape1(group_param[11], group_param[20]),
                           shape2 = shape2(group_param[11], group_param[20])),
                     rbeta(n = 3,               
                           shape1 = shape1(group_param[12], group_param[20]),
                           shape2 = shape2(group_param[12], group_param[20])),
                     rgamma(n = 3,               
                            shape = (group_param[13] / group_param[21]) ^ 2, 
                            scale = group_param[21] ^ 2 / group_param[13]),       
                     rgamma(n = 3,               
                            shape = (group_param[14] / group_param[21]) ^ 2, 
                            scale = group_param[21] ^ 2 / group_param[14]),             
                     runif(n = 3, min = 0, max = 0.658),
                     runif(n = 3, min = 0, max = 0.860),
                     runif(n = 3, min = 0, max = 1.260))) 
  if (model == "normal") {
    person_rho <- as.numeric(
      replicate(n = 3,
                expr = 
                  rcorrmatrix(d = 3,
                              alphad = group_param[22])[upper.tri(diag(3))]))
    person_param <- c(person_param, person_rho)
  }
  params <- c(group_param, person_param)
  params <- params + runif(n = length(params), min = .005, max = .015)
  return(params)
}

initialize_chains <- function(model, chain_n) {
  # Purpose: sample initial parameter vectors for all the chains
  # Input: character scalar model, integer scalar chain_n
  # Output: numeric matrix initial_state
  
  initial_state <- t(replicate(chain_n, sample_prior(model)))
  return(initial_state)
}


















