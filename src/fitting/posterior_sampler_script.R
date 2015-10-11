
# Sets the root directory and then loads posterior sampler
setwd("~/Dropbox/Slava/Masters/")
source("src/fitting/load_dependencies.R")

# Load training data and rescale RT to seconds
train_data <- read.table("data/train_data.txt") %>% 
  filter(subj == "kr", prop %in% c(0, 1)) %>%  # nh may be better
  dplyr::select(rt, resp, prop, instr) %>%
    mutate(rt = rt / 1000, prop = prop / 32)
# levels(train_data$resp) <- c("upper", "lower")
# train_data %<>% as_data_frame

# catch errors
errors <- file("results/fitting/kr-errors.Rout", "wt")
sink(errors, TRUE, "message")

# load chains to use last state as initial state
load("results/fitting/kr-posterior-chains-norm-fit.RData")

# Posterior sampler settings
# Note: theta_n = 20 for model = independent,
# theta_n = 23 for model = normal
model <- "normal"
theta_n <- 23
initial <- TRUE

# sampler settings
settings <- list(train_data = train_data,
                 theta_n = theta_n,
                 draw_n = 5,
                 model = model,
                 alpha = .15,
                 thread_n = 3, 
                 chunk_n = 1,
                 tol = 1e-3,
                 maxEvals = 5e6,
                 seed = 1800968452,
                 continue = FALSE)

# Saves machine/software/settings configuration for reproducibility
capture.output(sessionInfo(),
               file = paste0("results/fitting/kr-machine-software-config-norm-fit-",
                             Sys.Date(), ".txt"))
capture.output(settings[-1],
               file = paste0("results/fitting/kr-simulation-config-norm-fit-",
                             Sys.Date(), ".txt"))

# Obtains posterior sample for a model
attach(settings, 2, warn.conflicts = FALSE)
timer <- proc.time()
res <- sample_posterior()
timer <- proc.time() - timer

# Saves a array of posterior draws in plain text format
write.table(x = res, 
            file = "results/fitting/kr-posterior-chains-norm-fit.txt")

# Saves computational time
capture.output(timer, 
               file = paste0("results/fitting/kr-computational-time-norm-fit-",
                             Sys.Date(), ".txt"))


