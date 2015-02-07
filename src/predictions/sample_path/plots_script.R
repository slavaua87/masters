
# Sets the root directory and then uses relative paths
setwd("~/Dropbox/Slava/Masters/")
source("src/predictions/sample_path/plot_mean_paths_full.R")
load(file = "results/sample_path/paths2-2015-02-01.RData")

# Specifies plotting settings
settings <- list(cores = 3,
                 nrow = 6,
                 ncol = 4,
                 nplots = 22)

# Generates nplots .pdf images in the results directory for sample paths
# Ignore Error: Results...
with(data = settings, expr = plot_mean_paths(paths = paths, cores = cores, 
                                             nrow = nrow, ncol = ncol, 
                                             nplots = nplots))
 
  
