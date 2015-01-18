
### Prediction simulator
Simul.Pred <- function(theta, n, omega, a, b, model, tau) {
  # theta order: alpha, nu, eta, lambda, gamma, chi, phi, rhos
  # Container
  Preds <- matrix(0, n, 2)
  # Obtain copula points
  Pnts.cop <- Simul.Cop(n, theta, omega, model)
  # Map copula points into parameter values
  Params <- Map.Cop(Pnts.cop, theta, a, b)
  # Preds order: responses, response times
  # Simulate predictions
  Preds <- rw.vectorized(1, Params[, 2], .1, Params[, 3] * Params[1, 1],
                         Params[1, 1], Params[, 4], 0, tau)
  return(t(Preds))
}

############## Prediction results

m <- dim(Theta.a)[1]
n <- 10000
omega <- 5
a <- .2
b <- 1
model <- c('normal', 't', 'ind')
tau <- 1e-4
### Run the simulation
# Normal
Pred.n <- array(0, dim = c(n, 2, m))
prog.bar <- txtProgressBar(min = 0, max = m, style = 3)
for(i in 1:m) {
  Pred.n[, , i] <- Simul.Pred(Theta.a[i, ], n, omega, a, b, model[1], tau)
  setTxtProgressBar(pb = prog.bar, value = i)
}
save(Pred.n, file = 'Predictions_normal.RData')
# t
Pred.t <- array(0, dim = c(n, 2, m))
prog.bar <- txtProgressBar(min = 0, max = m, style = 3)
for(i in 1:m) {
  Pred.t[, , i] <- Simul.Pred(Theta.a[i, ], n, omega, a, b, model[2], tau)         
  setTxtProgressBar(pb = prog.bar, value = i)
}
save(Pred.t, file = 'Predictions_t.RData')
# Independent
#Pred.i <- array(0, dim = c(n, 2, m))
cl <- makeCluster(2)
registerDoParallel(cl)
Pred.i <- array(as.vector(
  foreach(param = iter(Theta.a, by = 'row'), 
          .multicombine = T, 
          .combine = cbind,
          .packages = c('copula', 'tmvtnorm')) %dopar% {
            Simul.Pred(param, n, omega, a, b, model[3], tau)}),
  dim = c(n, 2, m))
stopCluster(cl)

save(Pred.i, file = 'Predictions_independent.RData')




