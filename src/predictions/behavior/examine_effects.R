

qp <- read.table("results/behavior/behavior-summary-2015-03-30.txt") 

# examine probabilities
ind_diff <- function(x) round(x - x[1], 3)

filter(qp, quant_rank == 1) %>% 
  group_by(condition, model) %>% 
  summarise(prob = max(prob)) %>%
  mutate(diffs = ind_diff(prob)) %>% 
  ungroup %>%
  arrange(condition, diffs) %>%
  as.data.frame

filter(qp, quant_rank == 1) %>% 
  group_by(condition, model) %>% 
  summarise(prob = max(prob)) %>%
  mutate(diffs = ind_diff(prob)) %>% 
  ungroup %>%
  arrange(condition, diffs) %>%
  as.data.frame

# examine quantiles
quant_diff <- function(dat) {
  diffs <- dplyr::select(dat, quant) - 
    filter(dat, model == "independent") %>% 
    dplyr::select(quant) %>%
    unlist %>%
    rep(times = 3)
  unlist(diffs)
}

filter(qp, condition %in% (1:12)) %>%
  mutate(., diffs = quant_diff(.) * 1000) %>%
  dplyr::select(-1, -2)

# examine starting points
c(.221, .05) * .464 = 0.102544 0.023200

bright_spd <- data.frame(prop = c(.505, .520, .540, .565, .605, .740))

bright_acc <- data.frame(prop = c(.496, .509, .522, .541, .565, .600))

sp <- read.table("results/sample_path/initial_states.txt", header = TRUE)

group_by(sp, model, response, instruction, graph, condition, sum_stat) %>% 
  filter(model %in% c("independent", "normal", "t"),
         sum_stat %in% c("fast_mean", "middle_mean", "slow_mean"),
         condition == 32) %>%
  mutate(value = round(round(value, 4) / 0.0232 * 100 - 100, 2)) %>% 
  as.data.frame


# arrange by value to figure out largest effects
group_by(sp, model, response, instruction, graph, condition, sum_stat) %>% 
  filter(condition %in% c(61:72), model == "normal", instruction == 1,
         sum_stat %in% c("fast_mean", "middle_mean", "slow_mean")) %>%
  mutate(value = round(round(value, 4) / 0.0232 * 100 - 100, 2)) %>% 
  as.data.frame

