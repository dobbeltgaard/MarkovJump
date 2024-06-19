log.lik.cov <- function(m, s1, s2, u, z, beta, states, A_param = NULL, weight = NULL) {
  log.lik <- 0 # init log likelihood
  Cz <- exp(beta[(m):length(beta)] %*% t(z)) # precompute covariate contribution
  lambda <- exp(beta[1:(m - 1)]) # transform to working parameter
  if (is.null(weight)) {
    weight <- rep(1, length(u))
  } # if no weights, weigh all equal
  for (i in 1:length(u)) { # iterate through data
    A <- construct.A(m, lambda * Cz[i], A_param) # construct trans. rate matrix
    Pt <- s1[i] == states # compute state probs
    Ptu <- Pt %*% Matrix::expm(A * u[i]) # compute state probs forward in time u
    log.lik <- log.lik + weight[i] * log(Ptu[which(s2[i] == states)]) # update
  }
  return(-log.lik)
}



rps.cov <- function(m, s1, s2, u, z, beta, states) {
  Cz <- exp(beta[(m):length(beta)] %*% t(z)) # precompute covariate contribution
  lambda <- exp(beta[1:(m - 1)]) # transform to working parameter
  rps = 0
  res = 0
  for (i in 1:length(u)) { # iterate through data
    A <- construct.A(m, lambda * Cz[i]) # construct trans. rate matrix
    Pt <- s1[i] == states # compute state probs
    Ptu <- Pt %*% Matrix::expm(A * u[i]) # compute state probs forward in time u
    foo_state = which(s2[i] == states)
    for(k in 1:m){
      for(j in 1:k){
        rps = rps + Ptu[j] - Pt[j]
      }
      res = res + rps^2
      rps = 0
    }
  }
  return(res/length(u))
}


rps.cov <- function(m, s1, s2, u, z, beta, states) {
  Cz <- exp(beta[(m):length(beta)] %*% t(z)) # precompute covariate contribution
  lambda <- exp(beta[1:(m - 1)]) # transform to working parameter
  rps = 0
  res = 0
  for (i in 1:length(u)) { # iterate through data
    A <- construct.A(m, lambda * Cz[i]) # construct trans. rate matrix
    Pt <- s1[i] == states # compute state probs
    Ptu <- Pt %*% Matrix::expm(A * u[i]) # compute state probs forward in time u
    obs_state = s2[i] == states
    res = res + rps(Ptu, obs_state)
  }
  return(res/length(u))
}

rps.cov(m, d$`s-`, d$s, d$u/365, as.matrix(d[,exo.cols]), beta0, states)
rps.cov(m, d$`s-`, d$s, d$u/365, as.matrix(d[,exo.cols]), res3$lambda, states)

rps <- function(pred, obs){
  m = length(obs)
  rps = 0
  res = 0
  for(k in 1:m){
    for(i in 1:k){
      rps = rps + pred[i] - obs[i]
    }
    res = res + rps^2
    rps = 0
  }
  return(res)
}

library(RcppEigen); library(Rcpp)
sourceCpp("rps.cpp")


discrete_rps_cpp(m, d$`s-`, d$s, d$u/365, as.matrix(d[,exo.cols]), beta0)
discrete_rps_cpp(m, d$`s-`, d$s, d$u/365, as.matrix(d[,exo.cols]), res3$lambda)
discrete_rps_cpp(m, d$`s-`, d$s, d$u/365, as.matrix(d[,exo.cols]), res4$lambda)



rps(rep(1,6)/6, c(0,1,0,0,0,0))
rps_cpp(6, rep(1,6)/6, c(0,1,0,0,0,0))

rps(rep(1,5)/5, c(1,0,0,0,0))


