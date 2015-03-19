
obs <- read.table(file = "data/rr_processed.txt", sep = " ")

train_data <- bind_rows(filter(obs, block == 1), filter(obs, block == 3),
                        filter(obs, block == 5), filter(obs, block == 7))
test_data <- bind_rows(filter(obs, block == 2), filter(obs, block == 4),
                       filter(obs, block == 6), filter(obs, block == 8))
write.table(x = train_data, file = "data/train_data.txt",
            sep = " ", col.names = TRUE)
write.table(x = test_data, file = "data/test_data.txt",
            sep = " ", col.names = TRUE)

