
# weibull psychometric model
weibull_param <- data.frame(lower = c(-.352, -.565), upper = c(.329, .511),
                            scale = c(.526, .521), shape = c(4.413, 5.227),
                            row.names = c("acc", "spd"))

# brightness covariates 
bright_spd <- data.frame(prop = c(.512, .527, .547, .570, .605, .740))

bright_acc <- data.frame(prop = c(.503, .515, .530, .545, .565, .590))

# drift rate (reflects effect of s/a and brightness)
nu <- as.data.frame(t(data.frame(
  acc = as.vector(weibull(bright_acc, 
                          weibull_param["acc", "lower"],
                          weibull_param["acc", "upper"], 
                          weibull_param["acc", "scale"], 
                          weibull_param["acc", "shape"])),
  spd = as.vector(weibull(bright_spd, 
                          weibull_param["spd", "lower"],
                          weibull_param["spd", "upper"],
                          weibull_param["spd", "scale"], 
                          weibull_param["spd", "shape"])))))
colnames(nu) <- paste("prop", as.character(1:6), sep = "")

# wiener process (reflects effect of s/a)
wiener <- data.frame(alpha = c(.221, .05), eta = c(.127, .127),
                     lambda = c(.464, .464), gamma = c(.065, .065),
                     chi = c(.279, .279), phi = c(.041, .041), 
                     row.names = c("acc", "spd"))

# copula degrees of freedom and correlations
omega <- 3
rho <- data.frame(rho_db = c(.15, .5, .85, -.15, -.5, -.85, .15, .5, .85),
                  rho_dt = c(-.15, -.5, -.85, -.15, -.5, -.85, .15, .5, .85), 
                  rho_bt = c(-.15, -.5, -.85, .15, .5, .85, .15, .5, .85))
