
# root directory - all files are sourced based on relative paths
setwd("~/Masters/")

# load dependencies
source("src/fitting/load_dependencies.R")

# load chains to use last state as initial state
# load(file = "results/fitting/posterior-chains-test.RData")
# chain_last <- res[, , dim(res)[3]]

# experimental settings
model <- "normal"
smpl_size <- 1
prop <- seq(from = 1, to = 32, by = 7) / 32

# parameter sample from the prior
set.seed(245186)
theta <- initialize_chains(model = model, chain_n = 1) %>% as.numeric
theta_n <- length(theta)

# sample of data from the experiment
train_data <- sample_experiment(theta = theta, model = model,
                                prop = prop, smpl_size = smpl_size)

# sampler settings
settings <- list(train_data = train_data,
                 theta_n = theta_n,
                 chain_n = 3, 
                 draw_n = 10000,
                 model = model,
                 cores = 60, 
                 chunk = 1,
                 seed = 1800968452,
                 continue = FALSE,
                 chain_last = NULL)

# Obtains posterior sample for a model
timer <- proc.time()
res <- with(data = settings,
            expr = sample_posterior(train_data = train_data,
                                    theta_n = theta_n,
                                    chain_n = chain_n, 
                                    draw_n = draw_n,
                                    model = model,
                                    cores = cores,
                                    chunk = chunk,
                                    seed = seed,
                                    continue = continue,
                                    chain_last = chain_last))
timer <- proc.time() - timer

# Saves a array of posterior draws in plain text format
write.table(x = res, 
            file = "results/fitting/posterior-chains-norm-fit.txt")

# Saves computational time
capture.output(timer, 
               file = paste0("results/fitting/computational-time-norm-fit-",
                             Sys.Date(), ".txt"))

# Saves machine/software/settings configuration for reproducibility
capture.output(sessionInfo(),
               file = paste0("results/fitting/machine-software-config-norm-fit-",
                             Sys.Date(), ".txt"))
capture.output(settings[-1],
               file = paste0("results/fitting/simulation-config-norm-fit-",
                             Sys.Date(), ".txt"))



