# Simple PFs
lW_bpf <- function(x, obs, sigmaY){-0.5*(x - obs)^2/sigmaY}
lW_pf <- function(x, obs, sigmaY, x_past, opt_mean, opt_var){
  -0.5*(x - obs)^2/sigmaY - 0.5*(x - x_past)^2/sigmaX + 0.5*(x - opt_mean)^2/opt_var
}
lW_mpf <- function(x, obs, sigmaY, x_past, opt_mean, opt_var, W_past){
  lW <- rep(0, length(x))
  for (n in 1:length(x)) {
    lW[n] <- -0.5*(x[n] - obs)^2/sigmaY + log(sum(W_past*exp(-0.5*(x[n] - x_past)^2/sigmaX))) -
      log(sum(W_past*exp(-0.5*(x[n] - opt_mean)^2/opt_var)))
  }
  return(lW)
}

# Auxiliary PFs
W_apf <- function(x, obs, sigmaY, x_past, sigmaX, proposal_mean, proposal_var){
  W <- dnorm(obs, mean = x, sd = sqrt(sigmaY)) * dnorm(x, mean = x_past, sd = sqrt(sigmaX))/
    dnorm(x, mean = proposal_mean, sd = sqrt(proposal_var))
  return(W)
}
W_mapf <- function(x, x_past, obs, obs_future, proposal_mean, proposal_mean_future, sigmaX, sigmaY, proposal_var, W_past){
  W <- rep(0, length(x))
  for (n in 1:length(x)) {
    twisting_future <- optimal_auxiliary_function(x[n], obs_future, proposal_mean_future[n], sigmaX, sigmaY, proposal_var)
    twisting_present <- optimal_auxiliary_function(x_past, obs, proposal_mean, sigmaX, sigmaY, proposal_var)
    W[n] <- dnorm(obs, mean = x[n], sd = sqrt(sigmaY)) * twisting_future
    W[n] <- W[n] * sum(W_past * dnorm(x[n], mean = x_past, sd = sqrt(sigmaX))/twisting_present)/
      sum(W_past * dnorm(x[n], mean = proposal_mean, sd = sqrt(proposal_var)))
  }
  W <- W/sum(W)
  return(W)
}
W_mapf_last <- function(x, x_past, obs, proposal_mean, sigmaX, sigmaY, proposal_var, W_past){
  W <- rep(0, length(x))
  for (n in 1:length(x)) {
    twisting_present <- optimal_auxiliary_function(x_past, obs, proposal_mean, sigmaX, sigmaY, proposal_var)
    W[n] <- dnorm(obs, mean = x[n], sd = sqrt(sigmaY))
    W[n] <- W[n] * sum(W_past * dnorm(x[n], mean = x_past, sd = sqrt(sigmaX))/twisting_present)/
      sum(W_past * dnorm(x[n], mean = proposal_mean, sd = sqrt(proposal_var)))
  }
  W <- W/sum(W)
  return(W)
}
optimal_auxiliary_function <- function(x_past, obs, proposal_mean, sigmaX, sigmaY, proposal_var){
  composite_var <- (2/sigmaX + 2/sigmaY - 1/proposal_var)^(-1)
  constant <- sqrt(proposal_var*composite_var)/(2*pi*sigmaX*sigmaY)
  auxiliary <- exp(-0.5*(2*x_past^2/sigmaX + 2*obs^2/sigmaY-proposal_mean^2/proposal_var)) *
      exp(0.5*(2*x_past/sigmaX + 2*obs/sigmaY-proposal_mean/proposal_var)^2*composite_var)
  auxiliary <- (constant*auxiliary)^(1/2)
  return(auxiliary)
}

# Other useful functions
weights <- function(lW){
  max.lW <- max(lW)
  W <- exp(lW - max(lW))
  W <- W/sum(W)
  return(W)
}
multinomial_resample <- function(W, N, un){
  W <- c(W)
  indices <- rep(0, times = N)
  s <- W[1]
  m <- 1
  for(i in 1:N){
    while(s < un[i]){
      m <- m+1
      s <- s+W[m]
    }
    indices[i] <- m
  }
  return(indices)
}
