
cd("~/Dropbox/Slava/Masters/")
qp <- read.table("results/behavior/behavior-summary-2015-06-15.txt")
qp %<>% mutate(resp = rep(rep(c("bright", "dark"), each = 5), times = 3 * 108))

# examine probabilities
ind_diff <- function(x) round(x - x[1], 3)

filter(qp, quant_rank == 1) %>% 
  group_by(condition, model) %>% 
  summarise(prob = max(prob)) %>%
  mutate(diffs = ind_diff(prob)) %>% 
  ungroup %>%
  arrange(diffs) %>%
  as.data.frame

filter(qp, quant_rank == 1) %>% 
  group_by(condition, model) %>% 
  summarise(prob = max(prob)) %>%
  mutate(diffs = ind_diff(prob)) %>% 
  ungroup %>%
  arrange(condition, model, diffs) %>% 
  filter(model == "normal", condition %in% 73:108) %>%
  as.data.frame


# examine quantiles
quant_diff <- function(dat) {
  diffs <- dplyr::select(dat, quant) - 
    filter(dat, model == "independent") %>% 
    dplyr::select(quant) %>%
    unlist %>%
    rep(times = 3)
  unlist(diffs %>% multiply_by(1000) %>% round(0))
}

filter(qp, condition %in% (61:72)) %>%
  mutate(., diffs = quant_diff(.)) %>%
  dplyr::select(-1, -2) %>%
  filter(model %in% "normal") %>%
  arrange(condition, instruction, model)

# examine asymmetry
asymm_diff <- function(dat) {
  bright <- filter(dat, resp == "bright")
  dark <- filter(dat, resp == "dark")
  diffs <- (dplyr::select(bright, quant) - dplyr::select(dark, quant)) %>%
    multiply_by(1000) %>% round(0)
  return(bind_cols(dplyr::select(bright, -c(1, 2, 7)), asym = diffs))
}

filter(qp, condition %in% (61:72), model != "t") %>% 
  asymm_diff(.)

################ examine starting points

shape1 <- as.numeric(((1 - params["lambda"]) / params["gamma"] ^ 2 -
                        1 / params["lambda"]) * params["lambda"] ^ 2)
shape2 <- as.numeric(shape1 * (1 / params["lambda"] - 1))

x_sd <- sqrt(shape1 * shape2 / (shape1 + shape2) ^ 2 / (shape1 + shape2 + 1))

c(.221, .05) * .464 = 0.102544 0.023200

bright_spd <- data.frame(prop = c(.505, .520, .540, .565, .605, .740))

bright_acc <- data.frame(prop = c(.496, .509, .522, .541, .565, .600))

sp <- read.table("results/sample_path/initial_states.txt", header = TRUE)

acc <- .221
spd <- .05

start <- function(sp, inst) {
  diffs <- dplyr::select(sp, value) / inst
  unlist(diffs %>% round(3))
}

filter(sp, condition %in% c(67:72)) %>%
  mutate(., diffs = start(., spd)) %>%
  dplyr::select(-1, -6, -7) %>%
  filter(model %in% c("normal")) %>%
  mutate(diffs = round(round(diffs, 3) / .464 * 100 - 100, 2))

start_diff <- function(sp, inst) {
  diffs <- dplyr::select(sp, value) / inst - 
    filter(sp, model == "independent") %>% 
    dplyr::select(value) %>%
    divide_by(inst) %>%
    unlist %>%
    rep(times = 3) 
  unlist(diffs %>% round(4))
}

filter(sp, condition %in% c(49:54)) %>%
  mutate(., diffs = start_diff(., acc)) %>%
  dplyr::select(-1, -6, -7) %>%
  filter(model %in% c("normal")) 


# %>%
#   mutate(diffs = round(diffs / .464 * 100, 2)) %>% 
#   as.data.frame


# arrange by value to figure out largest effects
group_by(sp, model, response, instruction, graph, condition, sum_stat) %>% 
  filter(condition %in% c(61:72), model == "normal", instruction == 1,
         sum_stat %in% c("fast_mean", "middle_mean", "slow_mean")) %>%
  mutate(value = round(round(value, 4) / 0.0232 * 100 - 100, 2)) %>% 
  as.data.frame
