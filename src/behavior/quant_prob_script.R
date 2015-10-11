
# Sets the root directory; loads required functions
setwd("~/Masters/")
source("src/behavior/load_dependencies.R")

# catch errors
errors <- file(paste0("results/behavior/errors-report-", Sys.Date(), ".Rout"),
                      open = "wt")
sink(file = errors, append = TRUE, type = "message")

# Specifies simulation settings
settings <- list(models = c("independent", "normal", "t"),
                 smpl_size = 1e7,
                 sim_size = 1e3,
                 seeds = c(-1605742457, -1480525591, 1868731723),
                 sigma = 0.1,                 
                 cores = 35,
                 group_n = 54)

# Simulates behavior for the three models 
timer <- proc.time()
behavior_sum <- with(data = settings, 
                     expr = mapply(FUN = simul_behavior, 
                                   model = models, 
                                   seed = seeds, 
                                   MoreArgs = list(smpl_size, sim_size,
                                                   sigma, cores, group_n), 
                                   SIMPLIFY = FALSE)) %>% bind_rows
timer <- proc.time() - timer

# Saves a table of behavior summaries in plain text format
write.table(x = behavior_sum, 
            file = paste0("results/behavior/behavior-summary-", Sys.Date(),
                          ".txt"))

# Saves computational time
capture.output(timer, 
               file = paste0("results/behavior/computational-time-",
                             Sys.Date(), ".txt"))

# Saves machine/software/settings configuration for reproducibility
capture.output(sessionInfo(),
               file = paste0("results/behavior/machine-software-config-",
                             Sys.Date(), ".txt"))
capture.output(settings,
               file = paste0("results/behavior/simulation-config-",
                             Sys.Date(), ".txt"))









