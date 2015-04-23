



# Overall block performances over subjects for some session
rr_data %<>% group_by(subj, )
filter(rr_data, tr_id != 0, sess == 1) %>%
  ggplot(data = ., aes(x = rt, 
                       group = interaction(block, sess, subj), 
                       colour = subj)) + 
  facet_wrap(facets = ~block, nrow = 3, ncol = 3, scale = "free") + 
  geom_density() + 
  xlim(0, 1300)


est <- c(-16.6, 0,  17.7,  -32.9,  -50.0,  -66.9)  

exp(17.7) / sum(exp(est))



########### Data structure 
# 3 subjects, sess 1 is practice, instr 1 is acc

# Add labels

# Adjust response and add accuracy
rt.data[, 'resp'] <- rt.data[, 'resp'] + 1
rt.data[, 'acc'] <- as.numeric(rt.data[, 'dist'] == rt.data[, 'resp'])
# View structure
str(rt.data)
head(rt.data)

jf kr nh
########### Clean out non-trials
rt.data2 <- rt.data[rt.data[, 'tr_id'] != 0, ]
head(rt.data2)
########### Clean out extremely fast trials by subject 
rt.data3 <- rbind(rt.data2[rt.data2[, 'subj'] == 'jf' & 
                             rt.data2[, 'rt'] >= 200, ],
                  rt.data2[rt.data2[, 'subj'] == 'kr' & 
                             rt.data2[, 'rt'] >= 200, ],
                  rt.data2[rt.data2[, 'subj'] == 'nh' & 
                             rt.data2[, 'rt'] >= 200, ])

########### Clean out extremely slow trials
# Take data for a condition and reomves slow outliers 
slow.cut <- function(norm, rt.data3, cntr.tend, n.sd) {
  subj <- c('jf', 'kr', 'nh')
  rt.cond <- rt.data3[rt.data3[, 'subj'] == subj[norm[1]] & 
                        rt.data3[, 'instr'] == norm[2] & 
                        rt.data3[, 'prop'] == norm[3], ]
  loc <- mean(rt.cond[, 'rt'])
  spr <- sd(rt.cond[, 'rt']) * 4
  cut.off <- loc + spr
  rt.clean <- rt.cond[rt.cond[, 'rt'] < cut.off, ]
  return(rt.clean)
}
# Takes the whole dataset and removes slow outliers by condition
slow.clean <- function(rt.data3, cntr.tend, n.sd) {
  subj <- rep(1:3, each = 66)
  prop <- rep(c(0:32, 0:32),  3)
  instr <- rep(c(numeric(33), numeric(33) + 1), 3)
  norm <- matrix(c(subj, instr, prop), 198, 3)
  rt.data4 <- apply(X = norm, FUN = slow.cut, 1, 
                    rt.data3, cntr.tend, n.sd)
  return(rt.data4)
}
# Remove slow outliers by condition
rt.data4 <- slow.clean(rt.data3, mean, 3)
rt.final <- ldply(rt.data4)

########### Density of all data
windows()
plot(density(rt.final$rt),lwd=2)
rug(rt.final$rt)

########## Plots for masters
### Conditional distribution plots
# Take RTs given instructions and stimulus quality
Get.Rt.Plot1 <- function(exp.data, part, instr, ratio) {
  # Take out correct and error rts for 2 subjects
  # with some difficulty/instruction
  rts <- list(s1.cor = exp.data[exp.data$subj == part[1]  & 
                                  exp.data$instr == instr   &
                                  exp.data$prop == ratio &
                                  exp.data$acc == 1, 'rt'],
              s1.err = -exp.data[exp.data$subj == part[1]  & 
                                   exp.data$instr == instr   &
                                   exp.data$prop == ratio &
                                   exp.data$acc == 0, 'rt'],
              s2.cor = exp.data[exp.data$subj == part[2]  & 
                                  exp.data$instr == instr   &
                                  exp.data$prop == ratio &
                                  exp.data$acc == 1, 'rt'],
              s2.err = -exp.data[exp.data$subj == part[2]  & 
                                   exp.data$instr == instr   &
                                   exp.data$prop == ratio &
                                   exp.data$acc == 0, 'rt'])
  return(rts)
}
# Graph conditional distributions
Graph.Cond.Dist <- function(exp.data, part, instr, ratio) {
  # Get relevant rts
  rts <- Get.Rt.Plot1(exp.data, part, instr, ratio)
  # Initiate a plot and add lines
  plot(NULL, lwd = 2, col = 'black', xlim = c(-2500, 2500), ylim = c(0, .002),
       xlab = 'Response Times', ylab = 'Probability Density', main = NULL)
  lines(density(rts$s1.cor), lwd = 2, col = 'black')
  lines(density(rts$s2.cor), lwd = 2, col = 'blue')
  lines(density(rts$s1.err), lwd = 2, col = 'black')
  lines(density(rts$s2.err), lwd = 2, col = 'blue')
}
# Create a conditional distribution plot
part <- c('nh', 'jf')
instr <- 1
ratio <- 16
jpeg('Cond_distr_RRdata.jpg')
Graph.Cond.Dist(rt.final, part, instr, ratio)
graphics.off()

### Quantile probability plots
# Get data
Get.Rt.Plot2 <- Vectorize(FUN = function(exp.data, part, inst, ratio) {
  # Mirror proportions of pixels
  norm <- seq(0, 32, 1)[c(1 + ratio, 33 - ratio)]
  # Get full dataset for a particular accuracy level and instruction
  Cond.data <- rbind(exp.data[exp.data$subj == part &
                                exp.data$inst == inst &
                                exp.data$prop == norm[1],
                              c('acc', 'rt')],
                     exp.data[exp.data$subj == part &
                                exp.data$inst == inst &
                                exp.data$prop == norm[2],
                              c('acc', 'rt')])
  return(Cond.data)
}, vectorize.args = c('ratio'), SIMPLIFY = F)

# Calculate quantiles
Calc.Cond.Quants <- function(Cond.data, qs) {
  ### Description of arguments
  # cond.data - matrix
  # qs - vector
  
  # Upper responses index
  norm <- Cond.data[, 1] == 1
  # Quantiles for responses
  qs.upp <- quantile(x = Cond.data[norm, 2], probs = qs)
  qs.low <- quantile(x = Cond.data[!norm, 2], probs = qs)
  return(c(qs.upp, qs.low))
}

# Create quantile probability graph
Calc.Cond.Graph <- function(Cond.data, qs) {
  ### Description of arguments
  # exp.data - array
  # qs - vector
  
  # Quantile Counts and index
  qs.n <- length(qs)
  cond.n <- length(Cond.data)
  qs.norm <- c((cond.n + 1):(cond.n * 2), cond.n:1)
  
  
  # Calculate probability of a correct response for each condition
  cor.probs <- 0
  for (i in 1:cond.n) cor.probs[i] <- mean(Cond.data[[i]][, 1])
  # Get correct and error rates
  all.probs <- c(cor.probs, 1 - cor.probs)
  # Calculate quantiles for each condition
  quants <- matrix(0, 2 * qs.n, cond.n)
  for (i in 1:cond.n) quants[, i] <- Calc.Cond.Quants(Cond.data[[i]], qs)
  # Restructure quantiles matrix
  quants <- cbind(quants[1:qs.n, ], quants[-(1:qs.n), ])
  quants <- rbind(quants, all.probs)
  quants <- quants[, qs.norm]  
  return(quants)
}

# Graph the QP plot
Graph.QP <- function(exp.data, part, inst, ratio, 
                     qs, plot.main) {
  ### Description of arguments
  # exp.data - array
  # qs - vector
  
  # Get RT data
  Cond.data <- Get.Rt.Plot2(exp.data, part, inst, ratio)
  
  # Lines for 2 conditions
  rt.graph <- Calc.Cond.Graph(Cond.data, qs)
  
  # Plotting options
  qs.n <- length(qs)
  min.rt <- min(rt.graph[-(qs.n + 1), ])
  max.rt <- max(rt.graph[-(qs.n + 1), ])
  #c(min.rt - 50, max.rt + 50)
  # Qpp plot
  plot(NULL, NULL, main = plot.main,
       ylab = 'Reaction time (ms)', xlab = 'Response Probability', 
       xlim = c(0, 1), ylim = c(min.rt - 50, max.rt + 50))
  for (i in 1:qs.n) { 
    lines(rt.graph[qs.n + 1, ], rt.graph[i, ], type = 'o', pch = 20)
  }
}

# Make a QP graph
ratio <- c(5, 9, 11, 13, 14)
qs <- c(.1, .3, .5, .7, .9) 
part <- 'jf'
windows()
par(mfrow = c(1, 2))
Graph.QP(rt.final, part, 0, ratio, qs, 'Speed')
Graph.QP(rt.final, part, 1, ratio, qs, 'Accuracy')

### Auto-correlation plots
# Graph combining block sequence and acf plot
Graph.Ser.Acf <- function(exp.data, part, sess, block) {
  # Block index
  block.norm <- c(((block[1] - 1) * 100 + block[1]):(block[1] * 100 + block[1]),
                 ((block[2] - 1) * 100 + block[2]):(block[2] * 100 + block[2]))
  # Get a block of RTs   
  rts <- exp.data[exp.data$subj == part  & 
                    exp.data$sess == sess, 'rt'][block.norm]
  plot(acf(rts, plot = F), ylab = 'Auto-correlation', 
       main = 'Auto-correlation within a Block')
  plot.ts(rts, xlab = 'Trial', ylab = 'Response Time', 
          main = 'Trial Series of Response Times', xlim = c(0, 35))
  print(rts)
}
# Plot settings
part <- 'jf'
sess <- 2
block <- c(7, 8)
windows()
par(mfrow = c(1, 2))
Graph.Ser.Acf(rt.data2, part, sess, block)
########## Check variability among participants across conditions
# Extract relevant response times
Ind.Rt <- function(rt.final, part, instr, ratio) {
  # Create 4 matrices for each combination of conditions
  rt <- list(subj1.cor = c(rt.final[rt.final$subj == part[1]  & 
                                      rt.final$instr == instr       &
                                      rt.final$prop == ratio[1] &
                                      rt.final$acc == 1, 7],
                           rt.final[rt.final$subj == part[1]  & 
                                      rt.final$instr == instr       &
                                      rt.final$prop == ratio[2] &
                                      rt.final$acc == 1, 7]),
             subj1.err = c(rt.final[rt.final$subj == part[1]  & 
                                      rt.final$instr == instr       &
                                      rt.final$prop == ratio[1] &
                                      rt.final$acc == 0, 7],
                           rt.final[rt.final$subj == part[1]  & 
                                      rt.final$instr == instr      &
                                      rt.final$prop == ratio[2] &
                                      rt.final$acc == 0, 7]),
             subj2.cor = c(rt.final[rt.final$subj == part[2]  & 
                                      rt.final$instr == instr       &
                                      rt.final$prop == ratio[1] &
                                      rt.final$acc == 1, 7],
                           rt.final[rt.final$subj == part[2]  & 
                                      rt.final$instr == instr       &
                                      rt.final$prop == ratio[2] &
                                      rt.final$acc == 1, 7]),
             subj2.err = c(rt.final[rt.final$subj == part[2]  & 
                                      rt.final$instr == instr       &
                                      rt.final$prop == ratio[1] &
                                      rt.final$acc == 0, 7],
                           rt.final[rt.final$subj == part[2]  & 
                                      rt.final$instr == instr       &
                                      rt.final$prop == ratio[2] &
                                      rt.final$acc == 0, 7])
  )
  return(rt)                                               
}

# Graph a single data set
Graph.Rt <- function(rt.cond, lim) {
  plot(density(rt.cond$subj1.cor), xlab = 'Reaction Time', ylab = 'Density',
       lwd = 2, main = '', ylim = c(0, lim), xlim = c(0, 2000))
  lines(density(rt.cond$subj1.err), lty = 2, lwd = 2)
  lines(density(rt.cond$subj2.cor), lwd = 2, col = 'blue')
  lines(density(rt.cond$subj2.err), lty = 2, lwd = 2, col = 'blue')
}

# Graph response times
Graph.Var <- function(rt.final, ratio, part, lim) {
  # Obtain relevant data
  extreme.speed <- Ind.Rt(rt.final, part, 1, ratio[3:4])
  extreme.acc <- Ind.Rt(rt.final, part, 0, ratio[3:4])
  neutral.speed <- Ind.Rt(rt.final, part, 1, ratio[1:2])
  neutral.acc <- Ind.Rt(rt.final, part, 0, ratio[1:2])
  # Create a matrix of plots
  windows()
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
  Graph.Rt(extreme.speed, lim[1])
  title('Extreme - Accuracy')
  Graph.Rt(extreme.acc, lim[2])
  title('Extreme - Speed')
  Graph.Rt(neutral.speed, lim[3])
  title('Neutral - Accuracy')
  Graph.Rt(neutral.acc, lim[4])
  title('Neutral - Speed')
  title(main = 'Figure 1 - Variability in Reaction Times', 
        outer = T, cex.main = 1.5)
}

# Graphing settings
ratio <- c(16, 18, 3, 30)
part <- c('jf', 'kr')
lim <- c(.006, .013, .0018, .010)
# Graph
Graph.Var(rt.final, ratio, part, lim)

### Plot time series
windows()
plot.ts(rt.final[rt.final$subj == 'jf' & rt.final$sess == 2, 7], 
        ylab = 'Reaction Time', xlab = 'Trial', 
        main = 'Figure 2 - Across-trial Reaction Time Series')
lines(ts(rt.final[rt.final$subj == 'kr' & rt.final$sess == 2, 7]), col = 'blue')
plot(ts(rt.final[rt.final$subj == 'nh' & rt.final$sess == 3, 7]), col = 'pink')
acf(rt.final[rt.final$subj == 'nh' & rt.final$sess == 2, 7][1:360], lag = 100)
ccf(x = rt.final[rt.final$subj == 'nh' & rt.final$sess == 2, 7][1:180],
    y = rt.final[rt.final$subj == 'nh' & rt.final$sess == 2, 7][1:180])
plot(density(-rt.final[rt.final$subj == 'jf' & rt.final$sess == 2, 7][1:360]),
     xlim = c(-2000, 2000), ylim = c(0, .003))
lines(density(rt.final[rt.final$subj == 'jf' & rt.final$sess == 3, 7][1:360]))
### Plot autocorrelation
windows()
plot.ts(cbind(rt.final[rt.data2$sess == 2 & rt.data2$subj == 'jf', 'rt'],
              rt.final[rt.data2$sess == 2 & rt.data2$subj == 'jf', 'acc']),
        plot.type = 'single')

acf(cbind(rt.final[rt.data2$sess == 2 & rt.data2$subj == 'jf', 'rt'],
          rt.final[rt.data2$sess == 2 & rt.data2$subj == 'jf', 'acc']), 
    lag.max = 100)

ccf(rt.final[rt.data2$sess == 2 & rt.data2$subj == 'jf', 'rt'],
    rt.final[rt.data2$sess == 2 & rt.data2$subj == 'jf', 'acc'])

cor.test(rt.data2[rt.data2$sess == 2 & rt.data2$subj == 'jf', 'rt'],
         rt.data2[rt.data2$sess == 2 & rt.data2$subj == 'jf', 'acc'],
         method = 'pearson')


#
#
#### Latency - probability plots ####
#
#

# Explore stimulus proportions


# Proportions of white to black
seq(0,1,length.out=33)*100
sort(unique(rt.final$prop))

# Means
h.mu <- .625
l.mu <- .375
h.sd <- l.sd <- .1875

########### Features of the processed data
trials <- table(subj)
trials

########### Sampling algorithm functions


# subjects overlaid for speed/accuracy and 2 stimulus conditions combinations
#- 4 panes



require(tmvtnorm)

dom <- seq(0, 1, .01)

windows()
plot(dom, dtmvnorm.marginal(xn = dom, n = 1, mean = .25, sigma = .035, 
                            lower = .15, upper = 1), 
     type = 'l', lwd = 2, xlim = c(0, 1))


vioplot(rnorm(1000), rgamma(1000, 1, 1), col = 'yellow', lwd = 1.5, 
        rectCol = 'purple')


####### Full experiment plot

Extract.Responses <- function(rt.data2) {
  resp.set <- 1:2
  n <- length(rt.data2$rt)
  resps <- numeric(n)
  for (i in 1:n) {
    if (rt.data2$acc[i] == 0) resps[i] <- resp.set[-rt.data2$dist[i]]
    else resps[i] <- resp.set[rt.data2$dist[i]]
  }
  return(resps)
}

all.data <- rt.data2[, c(6, 7)]
all.data[all.data$acc == 0, 1] <- 2
all.data[all.data$acc == 1, 1] <- 1

exp.series <- cumsum(all.data$rt)
exp.responses <- Extract.Responses(rt.data2)
windows()
# Responses and response times
plot(exp.series[1:101], exp.responses[1:101], ylim = c(-.5, 3.5), pch = 19)
# Stimuli and response times
points(exp.series[1:101], rt.data2$dist[1:101], col = 'blue')
# Accuracy and response times
points(exp.series[1:100], all.data$acc[1:100], col = 'blue')



# Explore correlations between response times and responses
norm <- 502:601
cor(all.data$acc[norm], all.data$rt[norm])





#### probability integration


calc_prob_integrand <- function(dimensions, diffusion, params,
                                lower, upper, model) {
  # Purpose: Calculates integrand for the probability of choice associated
  # with the lower bound
  # Input: numeric scalars drift, bias, nondec, diffusion, 
  # numeric vector params, character scalar model
  # Output: numeric scalar integrand
  
  drift <- dimensions[1]
  bias <- dimensions[2]
  nondec <- dimensions[3]
  
  drift_trans <- drift / (1 - drift ^ 2)
  nondec_trans <- nondec / (1 - nondec)
  
  choice_prob <- pwiener(q = .Machine$double.xmax, 
                         alpha = params[1, "alpha"] / diffusion, 
                         tau = nondec_trans, beta = bias, 
                         delta = drift_trans / diffusion,
                         resp = "lower")
  param_dens <- calc_density(drift = drift_trans, bias = bias,
                             nondec = nondec_trans,
                             params = params, lower = lower,
                             upper = upper, model = model)
  jacobian <- (1 + drift ^ 2) / (1 - drift ^ 2) ^ 2 / (1 - nondec) ^ 2
  integrand <- choice_prob * param_dens * jacobian
  return(integrand)
}

integrate_prob <- function(diffusion, params, model, tol, maxEval) {
  # Purpose: Integrates choice probability with respect to parameters' density
  # Input: numerical scalar diffusion, numerical vector params, 
  # character scalar model, numerical scalar tol
  # Output: numerical scalar integral_out
  
  integral_out <- adaptIntegrate(f = calc_prob_integrand, 
                                 lowerLimit = c(-1, 0, 0), 
                                 upperLimit = c(1, 1, 1), 
                                 diffusion, choice, params,
                                 model, tol = tol, fDim = 1,
                                 maxEval = maxEval)
  return(integral_out)
}


