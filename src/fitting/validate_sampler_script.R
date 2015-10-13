
# root directory - all files are sourced based on relative paths
setwd("~/Masters/")

# load dependencies
source("src/fitting/load_dependencies.R")

# catch errors
errors <- file("results/fitting/errors.Rout", open = "wt")
sink(file = errors, append = TRUE, type = "message")

# load chains to use last state as initial state
#load("results/fitting/posterior-chains-norm-fit.RData")

# experimental settings
model <- c("independent", "normal")[2]
smpl_size <- 1e4
prop <- c(.5) #seq(from = 0, to = 32, by = 4) / 32
initial <- TRUE

# parameter sample from the prior
set.seed(245166)
theta <- sample_prior(model)
theta_n <- length(theta)

# sample of data from the experiment
train_data <- sample_experiment(theta, model,
                                prop, smpl_size)
initial <- TRUE

# sampler settings
settings <- list(train_data = train_data,
                 theta_n = theta_n,
                 draw_n = 1e5,
                 model = model,
                 alpha = .10,
                 thread_n = 60, 
                 chunk_n = 1,
                 tol = 1e-3,
                 maxEvals = 5e5,
                 seed = 1800968452)

# Saves machine/software/settings configuration for reproducibility
capture.output(sessionInfo(),
               file = paste0("results/fitting/machine-software-config-norm-fit-",
                             Sys.Date(), ".txt"))
capture.output(c(settings[-1], smpl_size = smpl_size, prop = list(prop)),
               file = paste0("results/fitting/simulation-config-norm-fit-",
                             Sys.Date(), ".txt"))

# Obtains posterior sample for a model
attach(settings, 2, warn.conflicts = FALSE)
timer <- proc.time()
res <- sample_posterior()
timer <- proc.time() - timer

# Saves a array of posterior draws in plain text format
write.table(x = res, 
            file = "results/fitting/posterior-chains-norm-fit.txt")

# Saves computational time
capture.output(timer, 
               file = paste0("results/fitting/computational-time-norm-fit-",
                             Sys.Date(), ".txt"))


