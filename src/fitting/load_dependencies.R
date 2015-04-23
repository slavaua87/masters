
# load (and install) libraries
if (!require("copula")) {
  install.packages("copula")
  library("copula")
}
if (!require("mvtnorm")) {
  install.packages("mvtnorm")
  library("mvtnorm")
}
if (!require("doParallel")) {
  install.packages("doParallel")
  library("doParallel")
}
if (!require("doRNG")) {
  install.packages("doRNG")
  library("doRNG")
}
if (!require("magrittr")) {
  install.packages("magrittr")
  library("magrittr")
}
if (!require("dplyr")) {
  install.packages("dplyr")
  library("dplyr")
}
if (!require("cubature")) {
  install.packages("cubature")
  library("cubature")
}
if (!require("RWiener")) {
  install.packages("RWiener")
  library("RWiener")
}
if (!require("clusterGeneration")) {
  install.packages("clusterGeneration")
  library("clusterGeneration")
}

# Source all the custom function relative to the working directory
source("src/fitting/calculate_weibull.R")
source("src/fitting/reflect.R")
source("src/fitting/check_support.R")
source("src/fitting/calc_copula_dens.R")
source("src/fitting/simulate_wiener_parameters.R")
source("src/fitting/simulate_rndwalk_rts.R")
source("src/fitting/calculate_joint.R")
source("src/fitting/integrate_pdf.R")
source("src/fitting/sample_experiment.R")
source("src/fitting/initiate_sampler.R")
source("src/fitting/update_chains.R")
source("src/fitting/calculate_joint_parallel.R")
source("src/fitting/sample_posterior.R")
source("src/fitting/sample_am.R")
