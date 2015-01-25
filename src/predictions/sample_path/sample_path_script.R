
# Sets the root directory
setwd("~/Dropbox/Slava/Masters/")
source("src/predictions/sample_path/simulate_sample_paths.R")

# Specifies simulation settings
settings <- list(model = c("independent", "normal", "t"),
                 smpl_size = 5,
                 seeds = c(1316048320, -1572737661, 195896225),
                 sigma = 1,
                 time_unit = 1e-4,
                 cores = 2)

# Simulates paths for the three models
paths <- mapply(FUN = simul_paths, model = settings$model, seed = settings$seeds, 
                MoreArgs = settings[c("smpl_size", "sigma", "time_unit", "cores")], 
                SIMPLIFY = FALSE)

# Saves a compressed nested list of paths in binary format
save(paths, file = "results/sample_path/paths.RData", compress = "gzip")
# Saves machine/software/settings configuration for reproducibility
capture.output(sessionInfo(), settings,
               file = "src/predictions/sample_path/reproduce_config.txt")
