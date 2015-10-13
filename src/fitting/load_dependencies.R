
# load packages
library("copula")
library("mvtnorm")
library("magrittr")
library("dplyr")
library("RWiener")
library("Rcpp")
library("RcppArmadillo")
library("integral")

# source all the custom function relative to the working directory
source("src/fitting/calculate_weibull.R") 
source("src/fitting/reflect.R") 
source("src/fitting/check_support.R")  
source("src/fitting/simulate_wiener_parameters.R") 
source("src/fitting/simulate_rndwalk_rts.R") 
source("src/fitting/sample_experiment.R") 
source("src/fitting/initiate_sampler.R") 
source("src/fitting/update_chains.R") 
source("src/fitting/sample_posterior.R") 
source("src/fitting/sample_am.R") 
