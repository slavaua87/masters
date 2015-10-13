
# Sets the root directory and then uses relative paths
setwd("~/Masters")
source("src/predictions/sample_path/plot_mean_paths.R")
load(file = "results/sample_path/paths-2015-06-20.RData")

# Specifies plotting settings 
settings <- list(cores = 3,
                 rows = 6,
                 cols = 2)

# Generates nplots .pdf images in the results directory for sample paths
# Ignore Error: Results...
timer <- proc.time()
with(data = settings, expr = plot_mean_paths(paths = paths, cores = cores, 
                                             rows = rows, cols = cols))
timer <- proc.time() - timer
# Captures time, settings and data used in a metadata file
capture.output(timer["elapsed"], settings, print("data = paths-2015-06-20.RData"),
               file = paste0("results/sample_path/plots-metadata-",
                             Sys.Date(), ".txt"))


