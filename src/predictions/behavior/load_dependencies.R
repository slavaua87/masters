
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
source("src/predictions/behavior/calculate_weibull.R")
source("src/predictions/behavior/wiener_parameters.R")
source("src/predictions/behavior/combine_parameters.R")
source("src/predictions/behavior/simulate_wiener_parameters.R")
source("src/predictions/behavior/simulate_rndwalk_rts.R")
source("src/predictions/behavior/calculate_summary.R")
source("src/predictions/behavior/simulate_behavior.R")
source("src/predictions/behavior/plot_behavior.R")

