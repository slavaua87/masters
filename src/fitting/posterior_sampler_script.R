
# Sets the root directory and then loads posterior sampler
setwd("~/Masters/")
source("src/fitting/load_dependencies.R")
#load("results/fitting/posterior-chains-norm-fit.RData")

#res <- partial_res[seq_len(9100), ]
# Load training data and rescale RT to seconds
train_data <- read.table("data/train_data.txt") %>% 
  mutate(rt = rt / 1000, prop = prop / 32)

# Posterior sampler settings
# Note: theta_n = 17 for model = independent, theta_n = 20 for model = normal
model <- "normal"
theta_n <- 20
# sampler settings
settings <- list(train_data = train_data,
                 theta_n = theta_n,
                 init_n = 5 * theta_n,
                 chain_n = 3, 
                 draw_n = 10000,
                 model = model,
                 cores_init = 60,
                 cores_update = 60, 
                 chunk_init = 1,
                 chunk_update = 1,
                 seed = 1800968452,
                 continue = FALSE,
                 history_last = NULL)

# Obtains posterior sample for a model
timer <- proc.time()
res <- with(data = settings,
            expr = sample_posterior(train_data = train_data,
                                    theta_n = theta_n,
                                    init_n = init_n,
                                    chain_n = chain_n, 
                                    draw_n = draw_n,
                                    model = model,
                                    cores_init = cores_init,
                                    cores_update = cores_update,
                                    chunk_init = chunk_init,
                                    chunk_update = chunk_update,
                                    seed = seed,
                                    continue = continue,
                                    history_last = history_last))
timer <- proc.time() - timer

# Saves a array of posterior draws in plain text format
# for combiningapply(res, 2, c)
# Saves a array of posterior draws in plain text format
save(x = res, 
     file = paste0("results/fitting/posterior-chains-norm-fit",
                   ".RData"))

# Saves computational time
capture.output(timer, 
               file = paste0("results/fitting/computational-time-norm-fit-",
                             Sys.Date(), ".txt"))

# Saves machine/software/settings configuration for reproducibility
capture.output(sessionInfo(),
               file = paste0("results/fitting/machine-software-config-norm-fit-",
                             Sys.Date(), ".txt"))
capture.output(settings[-c(1, 13)],
               file = paste0("results/fitting/simulation-config-norm-fit-",
                             Sys.Date(), ".txt"))
