

library(dplyr)
library(magrittr)
setwd("~/Dropbox/Slava/Masters/")
source(file = "src/fitting/clean_contaminants.R")

# Read in, remove non-data, rename columns
rr_processed <- read.table('data/rr98_data.txt', sep = ' ', header = F) %>% 
  select(-6, -10)
colnames(rr_processed) <- c('subj', 'sess', 'dist', 'prop', 
                            'instr', 'resp', 'rt', 'tr_id')

# Determine accuracy of choices
rr_processed %<>% mutate(resp = resp + 1) %>% 
  mutate(acc = as.integer(dist == resp))

# Remove duplicate session 11 from nh
norm <- (nrow(rr_processed) - 823):nrow(rr_processed)
rr_processed %<>% slice(-norm)

# Add a block variable
norm <- which(select(rr_processed, rt) == 32751)
rr_processed <- bind_rows(slice(rr_processed, 1:(norm - 1)),
                          slice(rr_processed, norm),
                          slice(rr_processed, -(1:(norm - 1))))

blocks_1 <- rep(x = seq_len(8), each = 103) %>% rep(times = 11)
blocks_2 <- c(rep(x = seq_len(8), each = 103) %>% rep(times = 7),
              rep(x = seq_len(7), each = 103),
              rep(x = seq_len(8), each = 103) %>% rep(times = 3))
blocks_3 <- rep(x = seq_len(8), each = 103) %>% rep(times = 11)
rr_processed %<>% mutate(block = c(blocks_1, blocks_2, blocks_3)) %>% 
  filter(tr_id != 0)

# Remove contaminants
rr_processed %<>% remove_slow(center_meas = median, n_sd = 3) %>%
  remove_fast(window_size = 33, p = 0.53)

# Remove non-responses
rr_processed %<>% filter(resp != 3)

# Save processed data
write.table(x = rr_processed, file = "data/rr_processed.txt", sep = " ",
            col.names = TRUE)







