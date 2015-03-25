
# Sets the root directory and then loads posterior sampler
setwd("~/Dropbox/Slava/Masters/")
source("src/fitting/sample_posterior.R")

# Load training data and rescale RT to seconds
train_data <- read.table("data/train_data.txt") %>% 
  mutate(rt = rt / 1000, prop = prop / 32)

# Posterior sampler settings
# Note: theta_n = 17 for model = independent, theta_n = 20 for model = normal
settings <- list(train_data = train_data,
                 theta_n = 17,
                 chain_n = 20,
                 draw_n = 3,
                 model = "independent",
                 cores = 1,
                 seed = 1800968452)

# Obtains posterior sample for a model
timer <- proc.time()
res <- with(data = settings,
            expr = sample_posterior(train_data = train_data,
                                    theta_n = theta_n, 
                                    chain_n = chain_n, 
                                    draw_n = draw_n,
                                    model = model,
                                    cores = cores,
                                    seed = seed))
timer <- proc.time() - timer

# Saves a array of posterior draws in plain text format
# for combiningapply(res, 2, c)
save(x = res, 
     file = paste0("results/fitting/posterior-chains", Sys.Date(),
                   ".txt"))

# Saves computational time
capture.output(timer, 
               file = paste0("results/fitting/computational-time-",
                             Sys.Date(), ".txt"))

# Saves machine/software/settings configuration for reproducibility
capture.output(sessionInfo(),
               file = paste0("results/fitting/machine-software-config-",
                             Sys.Date(), ".txt"))
capture.output(settings[[-1]],
               file = paste0("results/fitting/simulation-config-",
                             Sys.Date(), ".txt"))
