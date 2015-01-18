

# Algorithm to evaluate density FPT of the Wiener process

# Probality of the path being absorbed at the lower bound
pr.absorb = function(a, z, v, s) {
  # Rescale parameters
  a <- a / s
  z <- z / s
  v <- v / s
  # Calculate bias
  w = z / a
  # Overall probability of absorbing at a lower bound
  if (abs(v) < .Machine$double.xmin) return(1 - w)
  else return(pmin(1, (1 - exp(-2 * v * a * (1 - w))) / 
             (exp(2 * v * a * w) - exp(-2 * v * a * (1 - w)))))
}

# Small time iterations
small.iter = function(t,a,like.eps) {
  # Rescale time
  t.r = t/a^2
  # Minimum bound
  m.b = ceiling(1/pi/sqrt(t.r))
  # Find appropriate number of iterations
  cond = 2*sqrt(2*pi*t.r)*like.eps < 1
  if ( any(cond) ) 
    # Iterations above the lower bound
    return( pmax(sqrt(t.r)+1, ceiling(2 + sqrt(-2*t.r*log(pmax(.Machine$double.xmin,2*like.eps*sqrt(2*pi*t.r)))))) )
  if (any(!cond) )
    # Iterations equal to lower bound
    return( 2 )
}

# Large time iterations
large.iter = function(rt,a,like.eps) {
  # Rescale time
  rt.r = rt/a^2
  # Lower bound
  l.b = ceiling(1/pi/sqrt(rt.r))
  # Find appropriate number of iterations
  cond = pi*rt.r*like.eps < 1
  if ( any(cond) )
    # Proposed number
    return( pmax(l.b, ceiling(sqrt(-2*log(pmax(.Machine$double.xmin,pi*rt.r*like.eps))/pi^2/rt.r))) )
  if ( any(!cond) ) {
    # Else return lower bound
    return( l.b )
  }
}


# Large pdf component for the lower bound
largepdf.lower = function(t,k,a,z,v) {
  # Bias
  w = z/a
  # Rescaled time
  t.r = t/a^2
  # Summands
  n = 1:k
  # Scaling term
  s.t = 1/a^2*exp(-v*a*w - v^2*t/2)
  # Results 
  return( pmax(0, pi*sum(n*exp(-n^2*pi^2*t.r/2)*sin(n*pi*w))*s.t) )
}

# Small pdf component for the lower bound
smallpdf.lower = function(t,k,a,z,v) {
  # Bias
  w = z/a
  # Rescaled time
  t.r = t/a^2
  # Summands
  n = -floor((k-1)/2):ceiling((k-1)/2)
  # Scaling term
  s.t = 1/a^2*exp(-v*a*w - v^2*t/2)
  # Results 
  return(pmax(0, 1/sqrt(2*pi*t.r^3)*sum((w+2*n)*
         pmax(.Machine$double.xmin,exp(-(w+2*n)^2/2/t.r)))*s.t))
}

# Evaluate lower bound density
lower.dens = function(rt,a,z,v,s,like.eps) { 
  
  # Results container
  dens = rep(0,length(rt))

  # Rescale parameters given variance in diffusion
  a = a/s
  v = v/s
  z = z/s
  
  # Number of summands
  large = large.iter(rt,a,like.eps) 
  small = small.iter(rt,a,like.eps)
  for (i in 1:length(rt)) {
    if (large[i] <= small[i]) dens[i] <- largepdf.lower(rt[i],large[i],a,z,v)
    else dens[i] <- smallpdf.lower(rt[i],small[i],a,z,v)
  }
  
  # Return results
  return(dens)
}

# Evaluate upper bound density
upper.dens = function(rt, a, z, v, s, like.eps) {
  return(lower.dens(rt, a, a-z, -v, s, like.eps))
}

# Evaluate lower bound conditional PDF
lower.PDF = function(rt,a,z,v,s,like.eps) { 
  dens <- lower.dens(rt,a,z,v,s,like.eps)
  prob <- pr.absorb(a, z, v, s)
  return(dens / prob)
}

# Evaluate upper bound conditional PDF
upper.PDF = function(rt,a,z,v,s,like.eps) {
  dens <- upper.dens(rt,a,z,v,s,like.eps)
  prob <- 1 - pr.absorb(a, z, v, s)
  return(dens / prob)
}

# Vectorize FPT conditional pdfs
lower.PDF_vec <- Vectorize(FUN=lower.PDF,vectorize.args=list('a','z','v','rt'))
upper.PDF_vec <- Vectorize(FUN=upper.PDF,vectorize.args=list('a','z','v','rt'))

# Density for lower conditional RT = FPT + nondecision tau
lowerRT.PDF = function(rt,a,z,v,tau,s,like.eps) { 
  # Center with nondecision component
  t = rt - tau
  # Results container
  dens = rep(0,length(t))
  # Probability of the lower bound
  p = pr.absorb(a, z, v, s)
  # Rescale parameters given variance in diffusion
  a = a/s
  v = v/s
  z = z/s
  
  # Number of summands
  large = large.iter(t,a,like.eps) 
  small = small.iter(t,a,like.eps)
  for (i in 1:length(t)) {
    if ( large[i] <= small[i]) dens[i]=largepdf.lower(t[i],large[i],a,z,v)
    else dens[i]=smallpdf.lower(t[i],small[i],a,z,v)
  }
  
  # Return results
  return(dens/p)
}

# Density for upper conditional RT = FPT + nondecision tau
upperRT.PDF = function(rt,a,z,v,tau,s,like.eps) {
  
  # Call lower PDF and change parameters
  return(lowerRT.PDF(rt,a,a-z,-v,tau,s,like.eps))
}

# Vectorize conditional RT densities with nondecision tau
lowerRT.PDF_vec <- Vectorize(FUN = lowerRT.PDF, 
                             vectorize.args = list('a','z','v','tau','rt'))
upperRT.PDF_vec <- Vectorize(FUN = upperRT.PDF,
                             vectorize.args = list('a','z','v','tau','rt'))

# Density for lower part of RT = FPT + nondecision tau
lowerRT.dens <- function(rt, a, z, v, tau, s, like.eps) { 
  # Center data with nondecision parameter
  t <- rt - tau
  # Allocate memory for results
  dens <- vector(length = length(t))
  # Rescale parameters given standard deviation of diffusion
  a <- a / s
  v <- v / s
  z <- z / s
  
  # Number of summands in an integral approximation
  large <- large.iter(t, a, like.eps) 
  small <- small.iter(t, a, like.eps)
  # Use best approximation to calculate density of data
  for(i in 1:length(t)) {
    if(large[i] <= small[i]) dens[i] <- largepdf.lower(t[i], large[i], a, z, v)
    else dens[i] <- smallpdf.lower(t[i], small[i], a, z, v)
  }
  # Return vector of densities
  return(dens)
}

# Density for upper part of RT = FPT + nondecision tau
upperRT.dens <- function(rt, a, z, v, tau, s, like.eps) {  
  # Call lower PDF and change parameters
  return(lowerRT.dens(rt, a, a - z, -v, tau, s, like.eps))
}

# Vectorize RT densities with nondecision tau
lowerRT.dens_vec <- Vectorize(FUN = lowerRT.dens, 
                              vectorize.args = list('a', 'z', 'v', 'tau', 'rt'))
upperRT.dens_vec <- Vectorize(FUN = upperRT.dens, 
                              vectorize.args = list('a', 'z', 'v', 'tau', 'rt'))






