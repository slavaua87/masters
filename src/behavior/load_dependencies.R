
# packages
library("dplyr")
library("magrittr")
library("doParallel")
library("doRNG")
library("copula")
library("RWiener")
library("ggplot2")
library("ggthemes")

# custom functions
source("src/behavior/calculate_weibull.R")
source("src/behavior/wiener_parameters.R")
source("src/behavior/combine_parameters.R")
source("src/behavior/simulate_wiener_parameters.R")
source("src/behavior/simulate_rndwalk_rts.R")
source("src/behavior/simulate_behavior.R")
source("src/behavior/calculate_summary.R")
source("src/behavior/plot_behavior.R")

