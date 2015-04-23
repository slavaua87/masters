
# Sets the root directory and then uses relative paths
setwd("~/Dropbox/Slava/Masters/")
source("src/predictions/sample_path/plot_mean_paths.R")
load(file = "results/sample_path/paths-2015-04-02.RData")

# Specifies plotting settings
settings <- list(cores = 1,
                 rows = 6,
                 cols = 2)

# Generates nplots .pdf images in the results directory for sample paths
# Ignore Error: Results...
with(data = settings, expr = plot_mean_paths(paths = paths, cores = cores, 
                                             rows = rows, cols = cols))
 
  
