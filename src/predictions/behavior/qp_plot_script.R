
# Sets the root directory and then uses relative paths
setwd("~/Masters/")
source("src/predictions/behavior/load_dependencies.R")
behavior_sum <- 
  read.table(file = "results/behavior/behavior-summary-2015-05-12.txt",
             header = TRUE)

# Specifies plotting settings
settings <- list(rows = 3,
                 cols = 2)

# Generates nplots .pdf images in the results directory for sample paths
# Ignore Error: Results...
with(data = settings, expr = plot_behavior_sum(behavior_sum = behavior_sum, 
                                               rows = rows, cols = cols))
 
  
