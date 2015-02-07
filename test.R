

rndwalk <- function(smpl_size = 1, delta = .01, sigma = 1, beta = .75, 
                    alpha = 1.3, low_bound = 0, time_unit = .0001) {
  # Random walk simulation of a bounded Wiener process
  

test <- ind_param[1, ]

with(test, rndwalk_cmp(1, delta = nu, sigma = .1, .5, alpha,
            low_bound = 0, time_unit = 1e-3) %>% unlist %>% 
  data.frame(a = .) %T>%
  summarise(print(length(a))) %>%
  ggplot(., aes(x = seq(from = 1, by = 1, length.out = length(a)),
                y = a), ylim = .5) + geom_line() + theme_bw() + xlab("RT (ms)") +
  ylab("Evidence") + geom_hline(yintercept = c(0, alpha)))



sigma = .1
delta = .1
alpha = .2
time_unit = 1e-3

state_unit <- sigma * sqrt(time_unit)
prob.up <- .5 * (1 + delta * sqrt(time_unit) / sigma)
state_unit; prob.up


with(test, rndwalk_cmp(1, time_unit = 1e-3)) %>% unlist %>% 
  data.frame(a = .) %T>% 
  summarise(print(length(a))) %>%
  ggplot(., aes(x = seq(from = 1, by = 1, length.out = length(a)),
                y = a)) + geom_line() + theme_bw() + xlab("RT (ms)") +
  ylab("Evidence")

rw.basic.cmp <- cmpfun(f = rw.basic)

with(test, rw.basic.cmp(5000, v = nu, s = .1, z = alpha / 2,
                        a = alpha, ter = 0, tau = 1e-3) %>% as.data.frame %>% 
       rename(r = V1, rt = V2) %>%
  ggplot(, aes(group = r)) + facet_wrap(facets = ~r, nrow = 2) +
  geom_histogram(aes(x = rt), binwidth = .1))












