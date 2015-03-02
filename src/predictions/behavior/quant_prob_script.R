
# Sets the root directory and then loads behavioral data sampler
setwd("~/Dropbox/Slava/Masters/")
source("src/predictions/behavior/simulate_behavior.R")

# Specifies simulation settings
settings <- list(models = c("independent", "normal", "t"),
                 smpl_size = 1e3,
                 seeds = c(-1605742457, -1480525591,  1868731723),
                 sigma = .1,                 
                 cores = 1)

# Simulates behavior for the three models
timer <- proc.time()
behavior_sum <- with(data = settings, 
                     expr = mapply(FUN = simul_behavior, model = models, 
                                   seed = seeds, 
                                   MoreArgs = list(smpl_size, sigma, cores), 
                                   SIMPLIFY = FALSE)) %>% rbind_all
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








