
### Stores and generates parameters used for predictions in data.frames
# weibull psychometric model
weibull_param <- data.frame(lower = c(-.352, -.565), upper = c(.329, .511),
                            scale = c(.526, .521), shape = c(4.413, 5.227),
                            row.names = c("acc", "spd"))

# brightness covariate
bright <- data.frame(prop = c(.188, .25, .313, .375, .438, .5, 
                              .563, .625, .688, .750, .813))

# drift rate (reflects effect of s/a and brightness)
source("src/predictions/sample_path/calculate_weibull.R")
nu <- as.data.frame(t(data.frame(
                   acc = as.vector(weibull(bright, weibull_param["acc", "lower"],
                                 weibull_param["acc", "upper"], 
                                 weibull_param["acc", "scale"], 
                                 weibull_param["acc", "shape"])),
                   spd = as.vector(weibull(bright, weibull_param["spd", "lower"],
                                 weibull_param["spd", "upper"],
                                 weibull_param["spd", "scale"], 
                                 weibull_param["spd", "shape"])))))
colnames(nu) <- paste("prop", as.character(1:11), sep = "")

# wiener process (reflects effect of s/a)
wiener <- data.frame(alpha = c(.221, .05), eta = c(.127, .127),
                     lambda = c(.464, .464), gamma = c(.065, .065),
                     chi = c(.279, .279), phi = c(.041, .041), 
                     row.names = c("acc", "spd"))

# copula 
omega <- 5
rho <- data.frame(rho_db = c(.15, .5, .85, -.15, -.5, -.85, -.15, -.5, -.85),
                  rho_dt = c(-.15, -.5, -.85, .15, .5, .85, -.15, -.5, -.85), 
                  rho_bt = c(-.15, -.5, -.85, .15, .5, .85, .15, .5, .85))

