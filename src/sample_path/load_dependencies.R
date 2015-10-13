
# packages
library("compiler")
library("dplyr")
library("magrittr")
library("doParallel")
library("doRNG")

# support functions
source("src/sample_path/calculate_weibull.R") 
source("src/sample_path/calculate_path_stats.R") 
source("src/sample_path/wiener_parameters.R") 
source("src/sample_path/combine_parameters.R")  
source("src/sample_path/simulate_parameters.R") 
source("src/sample_path/simulate_rndwalk.R") 
source("src/sample_path/simulate_sample_paths.R") 
