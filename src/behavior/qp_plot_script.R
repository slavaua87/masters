
# sets the root directory and then uses relative paths
setwd("~/Masters/")
source("src/predictions/behavior/load_dependencies.R")
behavior_sum <- 
  read.table(file = "results/behavior/behavior-summary-2015-06-15.txt",
             header = TRUE)

# specifies plotting settings
settings <- list(rows = 3, cols = 2)

# generates nplots .pdf images in the results directory for sample paths
# ignore Error: Results...
with(data = settings, expr = plot_behavior_sum(behavior_sum = behavior_sum, 
                                               rows = rows, cols = cols))
 
  
