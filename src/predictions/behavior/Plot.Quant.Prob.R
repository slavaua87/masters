
#### Quantile-probability plot to examine joint distribution of response times
### and responses

# Required package
require(scales)

# Calculate quantiles
Calc.Cond.Quants <- function(Cond.data, qs) {
  ### Description of arguments
  # cond.data - matrix
  # qs - vector
  
  # Upper responses index
  ind <- Cond.data[, 1] == 1
  # Quantiles for responses
  qs.upp <- quantile(x = Cond.data[ind, 2], probs = qs)
  qs.low <- quantile(x = Cond.data[!ind, 2], probs = qs)
  return(c(qs.upp, qs.low))
}

# Graph quantile probability plots
Calc.Cond.Graph <- function(Cond.data, qs) {
  ### Description of arguments
  # exp.data - array
  # qs - vector
  
  # Quantile and Counts
  qs.n <- length(qs)
  cond.n <- dim(Cond.data)[3]
  
  # Calculate probability of a correct response for each condition
  cor.probs <- colMeans(Cond.data[, 1, ])
  all.probs <- c(cor.probs, 1 - cor.probs)
  all.probs <- all.probs[c(8:5, 1:4)]
  # Calculate quantiles for each condition
  quant <- apply(Cond.data, 3, Calc.Cond.Quants, qs = qs)
  # Restructure quantiles matrix
  quant <- cbind(quant[1:5, ], quant[-(1:5), ])
  quant <- quant[, c(8:5, 1:4)]  
  return(rbind(quant * 1000, sort(all.probs)))
}

Plot.Cond.Graph <- function(Cond.data.i, Cond.data.n, Cond.data.t, 
                            qs, plot.main) {
  ### Description of arguments
  # exp.data - array
  # qs - vector
  
  # Row number for probabilities
  x.ind <- (1:3) * length(qs) + 1:3
  # Lines for 3 models
  Lines.xy <- rbind(Calc.Cond.Graph(Cond.data.i, qs),
                    Calc.Cond.Graph(Cond.data.n, qs),
                    Calc.Cond.Graph(Cond.data.t, qs))
  X.coord <- Lines.xy[x.ind, ]
  Y.coord <- Lines.xy[-x.ind, ]
  # Plot settings
  min.rt <- min(Y.coord)
  max.rt <- max(Y.coord)
  xs <- rep(1:3, each = 5)
  
  # Qpp plot
  plot(NULL, NULL, main = plot.main,
       ylab = 'Reaction time (ms)', xlab = 'Response Probability', 
       xlim = c(0, 1), ylim = c(min.rt - 50, max.rt + 50))
  for (i in 1:(length(qs) * 3)) { 
    lines(X.coord[xs[i], ], Y.coord[i, ], type = 'o', pch = 20,
          col = c('black', alpha(c('red', 'blue'), 0.5))[xs[i]])
  }
  legend(x = 'topright', legend = c('i', 'n', 't'), 
         col = c('black', 'red', 'blue'), lty = 1)
}

Plot.Exp.Graph <- function(Theta.a, Pred.i, Pred.n, Pred.t, exp.n, qs) {
  ### Description of arguments
  # Theta.a - matrix
  # Pred.x - array
  # expr - vector (scalar)
  # qs - vector
  
  # Indeces
  pred.ind <- Theta.a[Theta.a[, 'expr'] == exp.n, 'index']
  cond.ind <- 1:(length(pred.ind) / 2)
  # Pull out observations relevant to the experiment
  Pred.i <- Pred.i[, , pred.ind]
  Pred.n <- Pred.n[, , pred.ind]
  Pred.t <- Pred.t[, , pred.ind]
  # Produce speed and accuracy conditions graphs
  windows()
  par(mfrow = c(1, 2))
  # Speed
  Plot.Cond.Graph(Pred.i[, , cond.ind], Pred.n[, , cond.ind],
                  Pred.t[, , cond.ind], qs, 'Speed Conditions')
  # Accuracy
  Plot.Cond.Graph(Pred.i[, , -cond.ind], Pred.n[, , -cond.ind],
                  Pred.t[, , -cond.ind], qs, 'Accuracy Conditions')
}






