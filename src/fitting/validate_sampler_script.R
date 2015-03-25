
setwd("~/Dropbox/Slava/Masters/")

source(file = "src/fitting/sample_experiment.R")
source(file = "src/fitting/sample_posterior.R")

# experimental settings
model <- "independent"
smpl_size <- 10
prop <- seq(from = 1, to = 32, by = 5) / 32

# parameter sample from the prior
set.seed(245187)
theta <- initialize_chains(model = model, chain_n = 1) %>% as.numeric
theta_n <- length(theta)

# sample of data from the experiment
train_data <- sample_experiment(theta = theta, model = model,
                                prop = prop, smpl_size = smpl_size)

# sampler settings
settings <- list(train_data = train_data,
                 theta_n = theta_n,
                 chain_n = 20,
                 draw_n = 3,
                 model = model,
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
save(x = res, 
     file = paste0("results/fitting/posterior-chains-test", Sys.Date(),
                   ".RData"))

# Saves computational time
capture.output(timer, 
               file = paste0("results/fitting/computational-time-test-",
                             Sys.Date(), ".txt"))

# Saves machine/software/settings configuration for reproducibility
capture.output(sessionInfo(),
               file = paste0("results/fitting/machine-software-config-test-",
                             Sys.Date(), ".txt"))
capture.output(settings[-1],
               file = paste0("results/fitting/simulation-config-test-",
                             Sys.Date(), ".txt"))



