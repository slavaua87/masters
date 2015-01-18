
# Clean workspace
rm(list=ls())
# Working directory
setwd("C://Users//Liberte//Dropbox//Slava//Masters//R_code//")
# Load relevant packages and functions
require(tmvtnorm) 
require(copula)
require(doParallel)
require('pingr')
source('PDF_FPT_ddm.R')
source('study_B_datasim.R')

############ Bayesian model components ############
### Likelihood
# Alpha
Eval.Like.Alp <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                          eta, lambda, lambda.n, inst, gamma, chi, 
                          phi, rho, omega, model, s, 
                          like.eps, obs.total, block) {
  source("C://Users//Liberte//Dropbox//Slava//Masters//R_code//PDF_FPT_ddm.R")
  # Half of observations
  obs.half <- 1:(obs.total / 2)
  # Index for upper responses for each alpha
  index1 <- Data[obs.half, 1] == 1
  index2 <- Data[-obs.half, 1] == 1
  # Get upper and lower rts for each alpha
  upper.rt1 <- Data[obs.half, 2][index1]
  lower.rt1 <- Data[obs.half, 2][!index1]
  upper.rt2 <- Data[-obs.half, 2][index2]
  lower.rt2 <- Data[-obs.half, 2][!index2]
  # Restructure given parameters
  Params <- matrix(params[-alpha], obs.total, 3, byrow = T)
  # Get upper and lower parameters
  Params.upper1 <- Params[obs.half, ][index1, ]
  Params.lower1 <- Params[obs.half, ][!index1, ]
  Params.upper2 <- Params[-obs.half, ][index2, ]
  Params.lower2 <- Params[-obs.half, ][!index2, ]
  # Calculate joint log densities for for performance data
  l <- log(c(upperRT.dens_vec(upper.rt1, params[alpha][1],
                              Params.upper1[, 2] * params[alpha][1],
                              Params.upper1[, 1], 
                              Params.upper1[, 3], s, like.eps),
             lowerRT.dens_vec(lower.rt1, params[alpha][1],
                              Params.lower1[, 2] * params[alpha][1],
                              Params.lower1[, 1], 
                              Params.lower1[, 3], s, like.eps),
             upperRT.dens_vec(upper.rt2, params[alpha][2],
                              Params.upper2[, 2] * params[alpha][2],
                              Params.upper2[, 1], 
                              Params.upper2[, 3], s, like.eps),
             lowerRT.dens_vec(lower.rt2, params[alpha][2],
                              Params.lower2[, 2] * params[alpha][2],
                              Params.lower2[, 1], 
                              Params.lower2[, 3], s, like.eps)))
  return(l)
} 
# Delta, beta, ter
Eval.Like.Del <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                          eta, lambda, lambda.n, inst, gamma, chi, 
                          phi, rho, omega, model, s, 
                          like.eps, obs.total, block) {
  params <- c(params[alpha], params[delta])
  # Upper bound response
  if(Data[1] == 1) return(log(upperRT.dens_vec(Data[2], params[1], 
                          params[3] * params[1], params[2], params[4],  
                          s, like.eps)))
  # Lower bound response
  else return(log(lowerRT.dens_vec(Data[2], params[1], 
              params[3] * params[1], params[2], params[4],  
              s, like.eps)))
}
### Parameters
# Alphas
Eval.Param.Alp <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                           eta, lambda, lambda.n, inst, gamma, chi, 
                           phi, rho, omega, model, s, 
                           like.eps, obs.total, block) {
  return(log(dunif(params[alpha], .01, .45)))
}
# Dependent
Eval.Param.Dep <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                           eta, lambda, lambda.n, inst, gamma, chi, 
                           phi, rho, omega, model, s, 
                           like.eps, obs.total, block) {
  # Trial updates
  if (length(delta) == 3) {
    Params <- matrix(params[delta], 1, 3)
    Data <- matrix(Data, 1, 2)
  }
  # Extract n parameters for (j,k)th condition
  else { 
    Params <- matrix(params[-alpha], length(params[-alpha]) / 3, 3) 
  }
  # Construct elliptical pdf
  copula.pdf <- ellipCopula(family = model, param = hparams[rho],
                            dim = 3, dispstr = 'un', 
                            df = hparams[omega])
  # Construction of the distribution
  copula.part <- dCopula(cbind(pnorm(Params[, 1], 
                                   rep(hparams[nu], each = nu.n),
                                   hparams[eta]),
                               pbeta(Params[, 2], 
                                     rep(Repam.Shape1(hparams[lambda[inst]],
                                                hparams[gamma[inst]]),
                                         each = lambda.n),
                                     rep(Repam.Shape2(hparams[gamma[inst]],
                                                hparams[gamma[inst]]),
                                         each = lambda.n)),
                               ptmvnorm(Params[, 3], 
                                        hparams[chi],
                                        hparams[phi],
                                        .1, Data[, 2])), copula.pdf)
  
  product.part <- c(dnorm(Params[, 1], 
                          rep(hparams[nu], 
                              each = nu.n),
                          hparams[eta]),
                    dbeta(Params[, 2], 
                          rep(Repam.Shape1(hparams[lambda[inst]],
                                     hparams[gamma[inst]]),
                              each = lambda.n),
                          rep(Repam.Shape2(hparams[lambda[inst]],
                                     hparams[gamma[inst]]),
                              each = lambda.n)),
                    dtmvnorm.vec(Params[, 3], 
                             hparams[chi],
                             hparams[phi], 
                             .1, Data[, 2]))
  return(log(c(copula.part, product.part)))
}
# Independent
Eval.Param.Ind <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                           eta, lambda, lambda.n, inst, gamma, chi, 
                           phi, rho, omega, model, s, 
                           like.eps, obs.total, block) {
  # Trial updates
  if (length(delta) == 3) {
    Params <- matrix(params[delta], 1, 3)
    Data <- matrix(Data, 1, 2)
  }
  # Extract n parameters for (j,k)th condition
  else {
    Params <- matrix(data = params[-alpha][delta], 
                   nrow = length(params[-alpha][delta]) / 3,
                   ncol = 3,
                   byrow = T)
  }
  # Select appropriate prior
  if (block %in% c('nu', 'eta')) {
  # Calculate density for each variable parameter
  drifts <- dnorm(Params[, 1], 
                  rep(hparams[nu], 
                      each = nu.n),
                  hparams[eta])
  return(log(drifts))
  }
  if (block == 'lambda') {
  biases <- dbeta(Params[, 2], 
                  rep(Repam.Shape1(hparams[lambda[inst]],
                             hparams[gamma[inst]]),
                      each = lambda.n),
                  rep(Repam.Shape2(hparams[lambda[inst]],
                             hparams[gamma[inst]]),
                      each = lambda.n))
  return(log(biases))
  }
  if (block == 'chi') {
  ters <- dtmvnorm.vec(xn = Params[, 3], n = 1,
                   mean = hparams[chi],
                   sigma = hparams[phi], 
                   lower = .1, upper = Data[, 2])
  return(log(ters))
  }
  # trial level update
  else {
    drifts <- dnorm(Params[, 1], 
                    rep(hparams[nu], 
                        each = nu.n),
                    hparams[eta])
    biases <- dbeta(Params[, 2], 
                    rep(Repam.Shape1(hparams[lambda[inst]],
                               hparams[gamma[inst]]),
                        each = lambda.n),
                    rep(Repam.Shape2(hparams[lambda[inst]],
                               hparams[gamma[inst]]),
                        each = lambda.n))
    ters <- dtmvnorm.vec(xn = Params[, 3], n = 1,
                         mean = hparams[chi],
                         sigma = hparams[phi], 
                         lower = .1, upper = Data[, 2])
    return(log(c(drifts, biases, ters)))
  }
}


### Hyperparameters
# Drifts
Eval.Hypparam.Nu <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                              eta, lambda, lambda.n, inst, gamma, chi, 
                              phi, rho, omega, model, s, 
                              like.eps, obs.total, block) {
  l <- log(dunif(hparams[nu], -1, 1))
  return(l)
}
# Drifts and sd
Eval.Hypparam.Eta <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                             eta, lambda, lambda.n, inst, gamma, chi, 
                             phi, rho, omega, model, s, 
                             like.eps, obs.total, block) {
  l <- log(dunif(hparams[eta], 0, .3))
  return(l)
}
# Nondecision time and sd
Eval.Hypparam.Chi <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                              eta, lambda, lambda.n, inst, gamma, chi, 
                              phi, rho, omega, model, s, 
                              like.eps, obs.total, block) {
  l <- log(c(dunif(hparams[chi], .1, 1),
             dunif(hparams[phi], 0, .35)))
  return(l)
}
# Biases and sd
Eval.Hypparam.Lam <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                              eta, lambda, lambda.n, inst, gamma, chi, 
                              phi, rho, omega, model, s, 
                              like.eps, obs.total, block) {
  l <- log(c(dunif(hparams[lambda[inst]], .1, .9),
             dunif(hparams[gamma[inst]], 0, .4)))
  return(l)
}
# Rhos
Eval.Hypparam.Rho <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                              eta, lambda, lambda.n, inst, gamma, chi, 
                              phi, rho, omega, model, s, 
                              like.eps, obs.total, block) {
  l <- log(dunif(hparams[rho], -1, 1))
  return(l)
}
# Omega
Eval.Hypparam.Ome <- function(Data, params, hparams, delta, alpha, nu, nu.n, 
                              eta, lambda, lambda.n, inst, gamma, chi, 
                              phi, rho, omega, model, s, 
                              like.eps, obs.total, block) {
  l <- log(c(dunif(hparams[rho], -1, 1), dgamma(hparams[omega], 1, 1)))
  return(l)
}

### Condition checker
# Takes a proposal and block information to choose relevant support limits
Check.Cond <- function(new.params, block, Data) {
  # nus
  if(block == 'nu') return(all(new.params > -1) &
                           all(new.params < 1))
  # eta
  if (block == 'eta') return(new.params > 0 &
                             new.params < .3)
  # lambda and gamma
  if(block == 'lambda') return(new.params[1] > .1 &
                               new.params[1] < .9 &
                               new.params[2] > 0 &
                               new.params[2] < .4)
  # chi and phi
  if(block == 'chi') return(new.params[1] > .1 &
                            new.params[1] < 1 &
                            new.params[2] > 0 &
                            new.params[2] < .35)
  # rho
  if(block == 'rho') return(new.params > -1 &
                            new.params < 1)
  # rho and omega
  if(block == 'omega') return(all(new.params[-4] > -1) &
                              all(new.params[-4] < 1) &
                              new.params[4] > 0 &
                              new.params[4] < Inf)
  # alpha
  if(block == 'alpha') return(all(new.params > 0.01) &
                              all(new.params < .45))
  # delta and beta and ter
  if(block == 'delta') return(new.params[2] > 0 &
                              new.params[2] < 1 &
                              new.params[3] > .1 &
                              new.params[3] < Data[2])
}

### Test ratios
Calc.Ratio <- function(densities, Data, params1, params2, hparams1, 
                       hparams2, delta, alpha, nu, nu.n, eta, lambda, 
                       lambda.n, inst, gamma, chi, phi, rho, 
                       omega, model, s, like.eps, obs.total, block) {
  # arbitrary parameter and hyperparameter densities
  l.new <- c(densities[[1]](Data, params = params1, hparams = hparams1, delta, 
                            alpha, nu, nu.n, eta, lambda, lambda.n, inst, 
                            gamma, chi, phi, rho, omega, model, s, 
                            like.eps, obs.total, block), 
             densities[[2]](Data, params = params1, hparams = hparams1, delta,
                            alpha, nu, nu.n, eta, lambda, lambda.n, inst, 
                            gamma, chi, phi, rho, omega, model, s, 
                            like.eps, obs.total, block))
  l.old <- c(densities[[1]](Data, params = params2, hparams = hparams2, delta, 
                            alpha, nu, nu.n, eta, lambda, lambda.n, inst, 
                            gamma, chi, phi, rho, omega, model, s, 
                            like.eps, obs.total, block),
             densities[[2]](Data, params = params2, hparams = hparams2, delta,
                            alpha, nu, nu.n, eta, lambda, lambda.n, inst,
                            gamma, chi, phi, rho, omega, model, s, 
                            like.eps, obs.total, block))
  
  # Combine and remove logs
  t.r <- exp(sum(l.new) - sum(l.old))
  # Test for NaN or infinity - skip if detected
  if (is.nan(t.r)) return(0)
  return(t.r)
}

### Proposal and acceptance functionrm
# Block updater returning a submatrix of estimates
Update.Param <- function(Updated, Fixed, Data, chain.n, gamma.low,
                         gamma.upp, eps.low, eps.upp, block, delta,
                         nu, nu.n, eta, lambda, lambda.n, alpha,
                         inst, gamma, chi, phi, rho, omega, model, 
                         s, like.eps, obs.total, density1, density2,
                         target, update.index) {  
  # Distribute chains across workers
  Updates <- foreach(chain = 1:chain.n, .combine = rbind, .multicombine = T,  
                     .packages = c('copula', 'tmvtnorm'),
                     .export = c('Calc.Eps', 'Check.Cond', 
                                 'Calc.Ratio', 'Repam.Shape1',
                                 'Repam.Shape2', 'dtmvnorm.vec')) %dopar% {
     # Proposal constant
     # runif(1, gamma.low, gamma.upp)                                 
     gama <- 2.38 / sqrt(2 * length(update.index) )                                
     # Sample a proposal satisfying some restrictions
     while (1 > 0) {
       # Get DE pieces
       old.params <- Updated[chain, ][update.index]
       index <- sample(c(1:chain.n)[-chain], 2, replace = F)
       m.params <- Updated[index[1], ][update.index]
       n.params <- Updated[index[2], ][update.index]
       # Propose a new theta
       new.params <- old.params + gama * 
       (m.params - n.params) + Calc.Eps(1, eps.low, eps.upp)
       # Check that proposal satisfies restrictions
       if (Check.Cond(new.params, block, Data)) break
     }
     # Create full old and new upddated parameters
     old <- Updated[chain, ]
     new <- replace(old, update.index, new.params)
     # Calculate the ratio of posteriors given the target of update
     if (target == 'param') {
       t.ratio <- Calc.Ratio(c(density1, density2), Data,
                             params1 = new,
                             params2 = old,
                             hparams1 = Fixed[chain, ],
                             hparams2 = Fixed[chain, ],
                             delta, alpha, nu, nu.n,
                             eta, lambda, lambda.n, inst, gamma, chi,
                             phi, rho, omega, model, s,
                             like.eps, obs.total, block)
     }
     else { 
       t.ratio <- Calc.Ratio(c(density1, density2), Data,
                             params1 = Fixed[chain, ],
                             params2 = Fixed[chain, ],
                             hparams1 = new,
                             hparams2 = old,
                             delta, alpha, nu, nu.n,
                             eta, lambda, lambda.n, inst, gamma, chi,
                             phi, rho, omega, model, s,
                             like.eps, obs.total, block)
     }
     #print(c('t.ratio', t.ratio))
     #if (t.ratio > .1) print(t.ratio)
     # calculate acceptance probability
     acc.prob <- min(1, t.ratio)
     #if (any(block == c('nu', 'eta'))) print(acc.prob)
     #if (acc.prob > .1 & block != 'delta') print(acc.prob)
     # Accept new params/hparams
     if (runif(1) < acc.prob) return(c(new, 1))
     # Reject new params/hparams
     else return(c(old, 0))
  }
  return(Updates)
}
### Support functions
# Uniform error added to proposal
Calc.Eps <- function(theta.n, low, upp) {
  x <- runif(theta.n, low, upp)
  return(x)
}
# Beta shape parameters
Repam.Shape1 <- function(mu, sigma) {
  return(((1 - mu) / sigma ^ 2 - 1 / mu) * mu ^ 2)
}
Repam.Shape2 <- function(mu, sigma) {
  return((((1 - mu) / sigma ^ 2 - 1 / mu) * mu ^ 2) * (1 / mu - 1))
}
# Beta transformation condition
Calc.Beta.Cond <- function(mu, sigma) {
  return((1 - mu) * mu > sigma ^ 2)
}
# Vectorize truncated normal
dtmvnorm.vec <- Vectorize(FUN = dtmvnorm.marginal, 
                          vectorize.args = list('xn', 'upper'))


# Initial value simulator
Simul.Inits <- function(obs.total, obs.cond, model, Data) {
  # Get hyperparameters
  while (1 > 0) {
  hparams <- c(runif(8, -.3, .3), runif(1, .05, .1), runif(1, .3, .7),
               runif(1, .04, .08), runif(1, .3, .7), runif(1, .04, .08),
               runif(1, .1, .5), runif(1, .05, .15))
  if (Calc.Beta.Cond(hparams[10], hparams[11]) &
      Calc.Beta.Cond(hparams[12], hparams[13])) break
  }
  # Add association parameters
  if (model == 'normal') hparams <- c(hparams, runif(3, -1, 1))
  if (model == 't') hparams <- c(hparams, runif(3, -1, 1), rgamma(1, 1, 1))
  # Get variable parameters
  params <- as.vector(t(cbind(
    rnorm(obs.total, rep(hparams[1:8], each = obs.cond), hparams[9]),
    rbeta(obs.total, rep(Repam.Shape1(hparams[c(10, 12)],
                                hparams[c(11, 13)]),
                                each = obs.cond * 4), 
                     rep(Repam.Shape2(hparams[c(10, 12)],
                                hparams[c(11, 13)]), 
                                each = obs.cond * 4)),
    apply(X = as.matrix(Data[, 2], obs.total, 1), MARGIN = 1, 
          FUN = rtmvnorm, n = 1, mean = hparams[14], 
          sigma = hparams[15], lower = .1, algorithm = 'gibbs', 
          start.value = .15))))
  # Add alpha 
  params <- c(params, runif(2, .01, .45))
  return(c(params, hparams))
}

# Data and intertrial parameter simulator
Simul.Data.Pars <- function(Thetas, n, omega, a, b, model, tau) {
  t.rt <- foreach(theta = iter(Thetas, by = 'row'), 
                  .combine = rbind, .multicombine = T,
                  .packages = c('copula', 'tmvtnorm'),
                  .export = c('Simul.Data', 'Simul.Cop', 'Map.Cop',
                              'qtmvnorm.vec', 'rw.vectorized'
                              )) %dopar% {
            Simul.Data(theta, n, omega, a, b, model, tau)
          }
  return(t.rt)
}

# Progress report
Report.Prog <- function(Posterior, acc.ind, nu, eta, lambda, 
                        gamma, chi, phi, alpha, delta, i, 
                        chain.n) {
  theta.ind <- c(alpha + c(nu, eta, lambda, gamma, chi, phi), 
                 alpha, delta)
  windows()
  plot.ts(t(Posterior[1, theta.ind, 1:i]), plot.type = 'multiple')
  acc.rate <- colSums(matrix(Posterior[, acc.ind, 1:i],
                      chain.n * i, length(acc.ind), byrow = T)) / chain.n / i
  return(acc.rate)
}


### Algorithms
# Blocked DE-mCMC -- Independent
Sample.Ind <- function(Posterior, Data, iter, hparam.ind, param.ind,
                       block.ind, chain.n, gamma.low,
                       gamma.upp, eps.low, eps.upp, delta, nu, 
                       nu.n, eta, lambda, lambda.n, inst, alpha,
                       gamma, chi, phi, rho, omega, model, 
                       s, like.eps, obs.total) {
  # Progress bar
  prog.bar <- txtProgressBar(min = 2, max = iter, style = 3)
  # Update the chains by block
  for (i in 2:iter) {
  #print('nu1')
    # nu1
    Posterior[, c(hparam.ind, block.ind[1]), i] <- Update.Param(
      Updated = Posterior[, hparam.ind, i - 1], 
      Fixed = Posterior[, param.ind, i - 1],
      Data = Data[1:(nu.n * 2), ], chain.n = chain.n, gamma.low = gamma.low, 
      gamma.upp = gamma.upp,
      eps.low = eps.low, eps.upp = eps.upp, block = 'nu', 
      delta = 1:(nu.n * 6), nu = 1:2, 
      nu.n = nu.n, eta = eta, lambda = lambda, lambda.n = lambda.n, 
      alpha = alpha, inst = inst, gamma = gamma, chi = chi, 
      phi = phi, rho = 0, omega = 0, model = model,
      s = 0, like.eps = 0, obs.total = 0, density1 = Eval.Param.Ind, 
      density2 = Eval.Hypparam.Nu, target = 'hparam', update.index = 1:2)
  #print('nu2')
  
  # nu2
  Posterior[, c(hparam.ind, block.ind[2]), i] <- Update.Param(
    Updated = Posterior[, hparam.ind, i], 
    Fixed = Posterior[, param.ind, i - 1],
    Data = Data[(nu.n * 2 + 1):(nu.n * 4), ], chain.n = chain.n, 
    gamma.low = gamma.low, gamma.upp = gamma.upp,
    eps.low = eps.low, eps.upp = eps.upp, block = 'nu', 
    delta = (nu.n * 6 + 1):(nu.n * 12), nu = 3:4, 
    nu.n = nu.n, eta = eta, lambda = lambda, lambda.n = lambda.n, 
    alpha = alpha, inst = inst, gamma = gamma, chi = chi, 
    phi = phi, rho = 0, omega = 0, model = model,
    s = 0, like.eps = 0, obs.total = 0, density1 = Eval.Param.Ind, 
    density2 = Eval.Hypparam.Nu, target = 'hparam', update.index = 3:4)
  
  #print('nu3')
  # nu3
  Posterior[, c(hparam.ind, block.ind[3]), i] <- Update.Param(
    Updated = Posterior[, hparam.ind, i], 
    Fixed = Posterior[, param.ind, i - 1],
    Data = Data[(nu.n * 4 + 1):(nu.n * 6), ], chain.n = chain.n, 
    gamma.low = gamma.low, gamma.upp = gamma.upp,
    eps.low = eps.low, eps.upp = eps.upp, block = 'nu', 
    delta = (nu.n * 12 + 1):(nu.n * 18), nu = 5:6, 
    nu.n = nu.n, eta = eta, lambda = lambda, lambda.n = lambda.n, 
    alpha = alpha, inst = inst, gamma = gamma, chi = chi, 
    phi = phi, rho = 0, omega = 0, model = model,
    s = 0, like.eps = 0, obs.total = 0, density1 = Eval.Param.Ind, 
    density2 = Eval.Hypparam.Nu, target = 'hparam', update.index = 5:6)
  
  #print('nu4')
  # nu4
  Posterior[, c(hparam.ind, block.ind[4]), i] <- Update.Param(
    Updated = Posterior[, hparam.ind, i], 
    Fixed = Posterior[, param.ind, i - 1],
    Data = Data[(nu.n * 6 + 1):(nu.n * 8), ], chain.n = chain.n, 
    gamma.low = gamma.low, gamma.upp = gamma.upp,
    eps.low = eps.low, eps.upp = eps.upp, block = 'nu', 
    delta = (nu.n * 18 + 1):(nu.n * 24), nu = 7:8, 
    nu.n = nu.n, eta = eta, lambda = lambda, lambda.n = lambda.n, 
    alpha = alpha, inst = inst, gamma = gamma, chi = chi, 
    phi = phi, rho = 0, omega = 0, model = model,
    s = 0, like.eps = 0, obs.total = 0, density1 = Eval.Param.Ind, 
    density2 = Eval.Hypparam.Nu, target = 'hparam', update.index = 7:8)
  
  #print('eta')
  # Drift and nu
  Posterior[, c(hparam.ind, block.ind[5]), i] <- Update.Param(
    Updated = Posterior[, hparam.ind, i], 
    Fixed = Posterior[, param.ind, i - 1],
    Data = Data, chain.n = chain.n, gamma.low = gamma.low, gamma.upp = gamma.upp,
    eps.low = eps.low, eps.upp = eps.upp, block = 'eta', 
    delta = -(obs.total * 3 + 1), nu = nu, 
    nu.n = nu.n, eta = eta, lambda = lambda, lambda.n = lambda.n, 
    alpha = alpha, inst = inst, gamma = gamma, chi = chi, 
    phi = phi, rho = 0, omega = 0, model = model,
    s = 0, like.eps = 0, obs.total = 0, density1 = Eval.Param.Ind, 
    density2 = Eval.Hypparam.Eta, target = 'hparam', update.index = eta)
  
  #print('lambda1')
  
    # Lambda and gamma 1
    Posterior[, c(hparam.ind, block.ind[6]), i] <- Update.Param(
      Updated = Posterior[, hparam.ind, i], 
      Fixed = Posterior[, param.ind, i - 1],
      Data = Data[1:(nu.n * 4), ], chain.n = chain.n, gamma.low = gamma.low, 
      gamma.upp = gamma.upp,
      eps.low = eps.low, eps.upp = eps.upp, block = 'lambda', 
      delta = 1:(nu.n * 12), 
      nu = 1:4, nu.n = nu.n, eta = eta, lambda = lambda, 
      lambda.n = lambda.n, alpha = alpha, inst = 1, gamma = gamma, chi = chi, 
      phi = phi, rho = 0, omega = 0, model = model,
      s = 0, like.eps = 0, obs.total = 0, density1 = Eval.Param.Ind, 
      density2 = Eval.Hypparam.Lam, target = 'hparam', update.index = 10:11)
  #print('lambda2')
  
    # Lambda and gamma 2
    Posterior[, c(hparam.ind, block.ind[7]), i] <- Update.Param(
      Updated = Posterior[, hparam.ind, i], 
      Fixed = Posterior[, param.ind, i - 1],
      Data = Data[(nu.n * 4 + 1):(nu.n * 8), ], chain.n = chain.n, 
      gamma.low = gamma.low, gamma.upp = gamma.upp,
      eps.low = eps.low, eps.upp = eps.upp, block = 'lambda', 
      delta = (nu.n * 12 + 1):(nu.n * 24), 
      nu = 5:8, nu.n = nu.n, eta = eta, lambda = lambda, 
      lambda.n = lambda.n, alpha = alpha, inst = 2, gamma = gamma, chi = chi, 
      phi = phi, rho = 0, omega = 0, model = model,
      s = 0, like.eps = 0, obs.total = 0, density1 = Eval.Param.Ind, 
      density2 = Eval.Hypparam.Lam, target = 'hparam', update.index = 12:13)
  #print('chi')
  
    # Chi and phi
    Posterior[, c(hparam.ind, block.ind[8]), i] <- Update.Param(
      Updated = Posterior[, hparam.ind, i], 
      Fixed = Posterior[, param.ind, i - 1],
      Data = Data, chain.n = chain.n, gamma.low = gamma.low, 
      gamma.upp = gamma.upp,
      eps.low = eps.low, eps.upp = eps.upp, block = 'chi', 
      delta = -(obs.total * 3 + 1), nu = nu, 
      nu.n = nu.n, eta = eta, lambda = lambda, 
      lambda.n = lambda.n, alpha = alpha, inst = inst, gamma = gamma, chi = chi, 
      phi = phi, rho = 0, omega = 0, model = model,
      s = 0, like.eps = 0, obs.total = 0, density1 = Eval.Param.Ind, 
      density2 = Eval.Hypparam.Chi, target = 'hparam', update.index = 14:15)
  #print('alpha')
  
    # Alpha 1, 2
    Posterior[, c(param.ind, block.ind[9]), i] <- Update.Param(
      Updated = Posterior[, param.ind, i - 1], 
      Fixed = Posterior[, hparam.ind, i],
      Data = Data, chain.n = chain.n, gamma.low = gamma.low, 
      gamma.upp = gamma.upp,
      eps.low = eps.low, eps.upp = eps.upp, block = 'alpha', delta = 0, nu = 0, 
      nu.n = 0, eta = 0, lambda = 0, 
      lambda.n = 0, alpha = alpha, inst = 0, gamma = 0, chi = 0, 
      phi = 0, rho = 0, omega = 0, model = model,
      s = s, like.eps = like.eps, obs.total = obs.total, 
      density1 = Eval.Like.Alp, density2 = Eval.Param.Alp, target = 'param', 
      update.index = alpha)
  #print('delta')
  
    # Delta, beta, tau
    for (j in 1:obs.total) {
      Posterior[, c(param.ind, block.ind[9 + j]), i] <- Update.Param(
        Updated = Posterior[, param.ind, i], 
        Fixed = Posterior[, hparam.ind, i],
        Data = Data[j, ], chain.n = chain.n, gamma.low = gamma.low, 
        gamma.upp = gamma.upp, eps.low = eps.low, eps.upp = eps.upp, 
        block = 'delta', delta = delta[j, ], nu = rep(nu, each = nu.n)[j], 
        nu.n = 1, eta = eta, lambda = lambda, lambda.n = 1, 
        alpha = alpha[1], 
        inst = rep(1:2, each = lambda.n)[j],
        gamma = gamma, chi = chi, phi = phi, rho = 0, omega = 0,
        model = model, s = s, like.eps = like.eps,
        obs.total = obs.total, density1 = Eval.Param.Ind, 
        density2 = Eval.Like.Del, target = 'param', update.index = delta[j, ])
    }
  # Progress report
  # Progress report
  if (i %in% seq(50, iter - 1, 50)) {
    Report.Prog(Posterior, block.ind[c(1, 5, 6, 8, 9, 10)], alpha[2], nu[1], eta, 
                lambda[1], gamma[1], chi, phi, delta[1, ], i, chain.n)
}
  #print('cycle done')
  # Update the progress bar
  setTxtProgressBar(pb = prog.bar, value = i)
  }
  return(Posterior)
}




########### Testing 
# Initial values
chain.n <- 25
obs.total <- 400
param.n <- obs.total * 3 + 17
block.n <- obs.total + 9
iter <- 250
param.ind <- 1:(obs.total * 3 + 2)
hparam.ind <- (max(param.ind) + 1):(max(param.ind) + 15)
block.ind <- (max(hparam.ind) + 1):(max(hparam.ind) + block.n)
gamma.low <- .5
gamma.upp <- 1
eps.low <- -.001
eps.upp <- .001
delta <- matrix(1:(3 * obs.total), obs.total, 3, byrow = T)
nu <- 1:8
nu.n <- obs.total / 8
eta <- 9
lambda <- c(10, 12)
lambda.n <- obs.total / 2
inst <- c(1, 2)
gamma <- c(11, 13)
chi <- 14
phi <- 15
rho = 0
omega = 0
alpha <- c(obs.total * 3 + 1, obs.total * 3 + 2)
s <- .1
like.eps <- 1e-8
model <- 'ind'
Thetas.test <- Theta.b[1:8, 1:7]

##### Get data and initial values
cl <- makeCluster(4)
registerDoParallel(cl)


Data.pars <- Simul.Data.Pars(Thetas.test, n = nu.n, omega = 5, 
              a = .2, b = 1, model = model, tau = 1e-4)
t.rt <- Data.pars[, 1:2]
Posterior <- array(0, dim = c(chain.n, param.n + block.n, iter))
Posterior[, 1:param.n, 1] <- t(replicate(chain.n, 
                               Simul.Inits(obs.total, nu.n, model, t.rt),
                               simplify = 'array'))



#Results.2 <- Results
#save(Results.1, Results.2, file = 'Results_1_2.RData')
#load('Results1.RData')
#Posterior[, , 1] <- Results


########## Run the posterior estimates


bench <- proc.time()
Results <- Sample.Ind(Posterior, t.rt, iter, hparam.ind, param.ind, block.ind,
                      chain.n, gamma.low,
                      gamma.upp, eps.low, eps.upp, delta, nu, 
                      nu.n, eta, lambda, lambda.n, inst, alpha,
                      gamma, chi, phi, rho, omega, model, 
                      s, like.eps, obs.total)
bench <- proc.time() - bench
ping(4)
stopCluster(cl)

########## Evaluate updating

### mixing
windows()
plot.ts(t(Results[1:10, 1201, ]), plot.type = 'multiple')
acf(t(Results[1:10, 242, ]), lag.max = 100)
plot(density(Results[, 1217, ]), lwd = 2)
abline(v = .05)

### acceptance rate
for (i in 1:length(block.ind)) 
  print(sum(Results[, block.ind[i], ]) / chain.n / iter)



# alpha
plot(seq(.01, .45, .005), dunif(seq(.01, .45, .005), .01, .45), lwd = 2,
     type = 'l', ylim = c(0, 3))
abline(v = .1, col = 'red', lwd = 2)
lines(density(Results[1:10, 241, 150:500]), lwd = 2)
# nu1
plot(seq(-1, 1, .01), dunif(seq(-1, 1, .01), -1, 1), lwd = 2,
     type = 'l', ylim = c(0, 3))
abline(v = -.2, col = 'red', lwd = 2)
lines(density(Results[1:10, 243, ]), lwd = 2)
#nu2
plot(seq(-1, 1, .01), dunif(seq(-1, 1, .01), -1, 1), lwd = 2,
     type = 'l', ylim = c(0, 3))
abline(v = .2, col = 'red', lwd = 2)
lines(density(Results[1:10, 246, ]), lwd = 2)
#nu3
plot(seq(-1, 1, .01), dunif(seq(-1, 1, .01), -1, 1), lwd = 2,
     type = 'l', ylim = c(0, 3))
abline(v = .025, col = 'red', lwd = 2)
lines(density(Results[1:10, 249, ]), lwd = 2)
#eta
plot(seq(0, .3, .01), dunif(seq(0, .3, .01), 0, .3), lwd = 2,
     type = 'l', ylim = c(0, 6))
abline(v = .1, col = 'red', lwd = 2)
lines(density(Results[1:10, 251, ]), lwd = 2)
#chi
plot(seq(.1, 1, .01), dunif(seq(.1, 1, .01), .1, 1), lwd = 2,
     type = 'l', ylim = c(0, 6))
abline(v = .60, col = 'red', lwd = 2)
lines(density(Results[1:10, 254, ]), lwd = 2)
#delta
plot(seq(-1, 1, .01), dunif(seq(-1, 1, .01), -1, 1), lwd = 2,
     type = 'l', ylim = c(0, 6))
abline(v = -.42, col = 'red', lwd = 2)
lines(density(Results[1:10, 1, ]), lwd = 2)
# beta
plot(seq(0, 1, .01), dunif(seq(0, 1, .01), 0, 1), lwd = 2,
     type = 'l', ylim = c(0, 6))
abline(v = .57, col = 'red', lwd = 2)
lines(density(Results[1:10, 2, ]), lwd = 2)
# tau
plot(seq(.1, .37, .01), dunif(seq(.1, .37, .01), .1, .37), lwd = 2,
     type = 'l', ylim = c(0, 6))
abline(v = .28, col = 'red', lwd = 2)
lines(density(Results[1:10, 3, ]), lwd = 2)




###################### COMPUTATIONS ########################

# Sampling settings
chain.n <- 20
iter <- 1500
# True parameter values
s <- 1
like.eps <- 1e-8
# Prior parameters
gamma.low <- .5
gamma.upp <- 1
eps.low <- -.0001
eps.upp <- .0001
p.means <- c(.5,0,0,-1.5)
p.sds <- c(.05,.3,.03,.05)
p.shapes <- c(.5,.5,.5,.5)
p.rates <- c(2,2,2,15)
# Data simulation
part.n <- 10
theta.n <- 4*part.n + 8
obs.n <- 150
time.step <- 1e-4
# Simulate parameters and data
grp_param <- group.sim(p.means,p.sds,p.shapes,p.rates)
ind_param <- ind.sim(grp_param,part.n)
obs <- dif.data(part.n,obs.n,ind_param,time.step,s)
obs <- cbind(obs,rep(1:part.n,each=obs.n))
tau_bound <- min(obs[,2])
all_true <- list(g.theta = grp_param,i.theta = ind_param,data = obs)
save(all_true,file='blocked_DE_1group_truevalues.RData')

# Initial values
param <- array(0,dim=c(chain.n,theta.n,iter))
init_group <- matrix(replicate(chain.n,exp=group.sim(c(.5,0,0,1.5*log(tau_bound)),
                                                     c(.05,.30,.03,.001),p.shapes,p.rates)),chain.n,8,byrow=T)
param[,,1] <- cbind(t(apply(init_group,1,FUN=t_ind.sim,part.n=part.n,min.rt=log(tau_bound))),
                    init_group)
# Full sampling mechanism
block.de <- function(proposals,chain.n,iter,part.n,gamma.low,gamma.upp,eps.low,eps.upp,p.means,
                     p.sds,p.shapes,p.rates,like.eps) {
  # Indexing hyper and individual parameters
  indxs.h <- 1:8+4*part.n
  indxs.i <- seq(0,4*part.n,4)
  # Iterations
  for(i in 2:iter) {
    # Hyperparameters update
    for(j in 1:4) {
      proposals[,indxs.h[c(j,4+j)],i] <- hyper.update(chain.n,proposals[,indxs.h[c(j,4+j)],i-1],
                                                      proposals[,seq(j,4*(part.n-1)+j,4),i-1],
                                                      gamma.low,gamma.upp,eps.low,eps.upp,p.means[j],
                                                      p.sds[j],p.shapes[j],p.rates[j])
    }
    # Individuals update
    for(k in 1:part.n) {
      proposals[,1:4+indxs.i[k],i] <- ind.update(chain.n,obs.n,obs,proposals[,1:4+indxs.i[k],i-1],
                                                 proposals[,indxs.h,i],gamma.low,gamma.upp,eps.low,eps.upp,like.eps,k)
    }
  }
  return(proposals)
}

wtc <- proc.time()
post.samples <- block.de(param,chain.n,iter,part.n,gamma.low,gamma.upp,eps.low,eps.upp,p.means,
                         p.sds,p.shapes,p.rates,like.eps)
proc.time() - wtc

save(post.samples,file='blocked_DE_1group_samples.RData')

load('blocked_DE_1group_samples.RData')
load('blocked_DE_1group_truevalues.RData')



for (i in 1:chain.n) {
print(sum(dtmvnorm.vec(xn = Posterior[i, seq(3, 240, 3), 1], n = 1,
             mean = Posterior[i, 256, 1],
             sigma = Posterior[i, 257, 1],
             lower = .1, upper = t.rt[, 2]) == 0))}


dtmvnorm.vec(x = Posterior[1, seq(3, 240, 3), 1], n = 1,
                  mean = Posterior[1, 256, 1],
                  sigma = Posterior[1, 257, 1],
                  lower = .1, upper = t.rt[1, 2])



x <- seq(.1, t.rt[2, 2], .01)
fx <- dtmvnorm.marginal(x = x, n = 1,
                        mean = Posterior[1, 256, 1],
                        sigma = Posterior[1, 257, 1], 
                        lower = .1, upper = t.rt[2, 2])

windows()
plot(x, fx, type = 'l')
abline(v = Posterior[1, 6, 1])


scatter.smooth(Results[, 1216, ], Results[, 1217, ])





params <- Results[1, 1:1202, 1]
params[1201:1202] <- c(.02, .02)

xy <- expand.grid(seq(.01, 1, .01), seq(.01, 1, .01))

alpha.post <- function(x, y, params) {
  params[1201:1202] <- cbind(x, y)
  return(sum(Eval.Like.Alp(t.rt, params, Results[1, 1202:1217, 1],
              delta, alpha, nu, nu.n, 
              eta, lambda, lambda.n, inst, gamma, chi, 
              phi, rho, omega, model, s, 
              like.eps, obs.total, block),  
  Eval.Param.Alp(t.rt, params, Results[1, 1202:1217, 1],
                 delta, alpha, nu, nu.n, 
                 eta, lambda, lambda.n, inst, gamma, chi, 
                 phi, rho, omega, model, s, 
                 like.eps, obs.total, block)))
}
x <- seq(.01, 1, .01)
y <- seq(.01, 1, .01)
z <- outer(x, y, alpha.post, params = params)


image(xy, z = z)



alpha.post(.05, .05, params)









