
# Sets the root directory and load necessary functions
setwd("~/Masters/")
source("src/sample_path/load_dependencies.R")

# catch errors
errors <- file("results/sample_path/errors.Rout", open = "wt")
sink(file = errors, append = TRUE, type = "message")

# Specifies simulation settings
settings <- list(models = c("independent", "normal", "t"),
                 smpl_size = 1e5,
                 sim_size = 1e3,
                 seeds = c(1316048320, -1572737661, 195896225),
                 sigma = .1,
                 time_unit = 1e-3,
                 cores = 50,
                 group_n = 72)

# Simulates paths for the three models
timer <- proc.time()
paths <- with(data = settings, 
              expr = mapply(FUN = simul_paths, model = models, seed = seeds, 
                            MoreArgs = list(smpl_size, sim_size, sigma,
                                            time_unit, cores, group_n), 
                            SIMPLIFY = FALSE))
timer <- proc.time() - timer

# Saves computational time  
capture.output(timer, 
               file = paste0("results/sample_path/paths-computational-time-",
                             Sys.Date(), ".txt"))
# Saves a compressed nested list of paths in binary format
save(paths, file = paste0("results/sample_path/paths-", Sys.Date(),
                          ".RData"))
# Saves machine/software/settings configuration for reproducibility
capture.output(sessionInfo(),
               file = paste0("results/sample_path/machine-software-config-",
                             Sys.Date(), ".txt"))
capture.output(settings, 
               file = paste0("results/sample_path/simulation-config-",
                             Sys.Date(), ".txt"))

