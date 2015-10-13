

setwd("~/Dropbox/Slava/Masters/")

library("gridExtra")
obs <- read.table(file = "~/Dropbox/Slava/Masters/data/rr_processed.txt",
                  sep = " ")
time_series <- ggplot() + 
  geom_line(aes(x = 1:length(rt), y = rt),
            filter(obs, subj == "nh", sess %in% c(3, 5, 7, 9)) %>%
            arrange(sess, block, tr_id)) + 
  theme_solarized_2() +
  xlab("trial") + 
  ylab("reaction time (ms)")
rt_prop <- ggplot() +
  geom_point(aes(x = prop / 32, y = rt),
             filter(obs, subj == "nh")) +
  theme_solarized_2() +
  xlab("brightness proportion") + 
  ylab("reaction time (ms)")

late_prop1 <- ggplot() + 
  stat_summary(aes(x = prop / 32, y = rt, group = prop),
              filter(obs, subj == "nh", instr == 1, resp == 1), fun.y = mean,
              fun.ymin = mean, fun.ymax = mean) + 
  theme_solarized_2() +
  xlab("brightness proportion") + 
  ylab("mean reaction time (ms)")

late_prop2 <- ggplot() + 
  stat_summary(aes(x = prop / 32, y = rt, group = prop),
               filter(obs, subj == "nh", instr == 0, resp == 1), fun.y = mean,
               fun.ymin = mean, fun.ymax = mean) + 
  theme_solarized_2() +
  xlab("brightness proportion") +
  ylab("mean reaction time (ms)")
pdf("results/exploratory/data_structure.pdf")
gridExtra::grid.arrange(time_series, rt_prop, late_prop1, late_prop2)
dev.off()
