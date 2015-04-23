
# Sets the root directory and then uses relative paths
setwd("~/Dropbox/Slava/Masters/")
source("src/predictions/sample_path/simulate_sample_paths.R")

# Specifies simulation settings
settings <- list(models = c("independent", "normal", "t"),
                 smpl_size = 50,
                 seeds = c(1316048320, -1572737661, 195896225),
                 sigma = .1,
                 time_unit = 1e-3,
                 cores = 2)

# Simulates paths for the three models
timer <- proc.time()
paths <- with(data = settings, 
              expr = mapply(FUN = simul_paths, model = models, seed = seeds, 
                            MoreArgs = list(smpl_size, sigma, time_unit, cores), 
                            SIMPLIFY = FALSE))
timer <- proc.time() - timer

# Saves computational time
capture.output(timer, 
               file = paste0("results/sample_path/computational-time-",
                             Sys.Date(), ".txt"))

# Saves a compressed nested list of paths in binary format
save(paths, file = paste0("results/sample_path/paths-", Sys.Date(),
                          ".RData"), compress = "gzip")

# Saves machine/software/settings configuration for reproducibility
capture.output(sessionInfo(),
               file = paste0("results/sample_path/machine-software-config-",
                             Sys.Date(), ".txt"))
capture.output(settings,
               file = paste0("results/sample_path/simulation-config-",
                             Sys.Date(), ".txt"))

