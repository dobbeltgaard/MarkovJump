

# function for constructing transition rate matrix with r(r-1)/2 free params
construct.A1 <- function(m, lambda) {
  count <- 0 # init counter
  A <- matrix(0, ncol = m, nrow = m) # init A
  for (i in 1:m) { # iterate through rows
    for (j in 1:m) { # iterate through cols
      if (i < j) { # Set upper triangle except diagonal to lambda
        count <- count + 1 # keep track of inputing elements
        A[i, j] <- lambda[count] # assign elements
      }
    }
  }
  diag(A) <- -rowSums(A) # diag elements are equal to negative total rate of transmission
  return(A)
}
# example::: construct.A(m = 5, lambda = 1:10)


construct.A2 <- function(m, lambda) {
  count = 0 # init counter
  A = construct.A3(m,lambda)
  for(i in 1:m){
    for(j in 1:m){
      if(j-i > 1){
        #A[i,j] = 1/(1/A[i,j-1] + 1/A[j-1,j])
        A[i,j] = (A[i,j-1] * A[j-1,j]) / (A[i,j-1] + A[j-1,j])
      }
    }
  }
  diag(A) = 0 # diag elements are equal to negative total rate of transmission
  diag(A) = -rowSums(A) # diag elements are equal to negative total rate of transmission
  
  return(A)
}

# construct.A2(m = 5, c(1,1,1,1))

# function for constructing transition rate matrix with lambda_i,j =0 where j-i>1
construct.A3 <- function(m, lambda) {
  count <- 0 # init counter
  A <- matrix(0, ncol = m, nrow = m) # init A
  for (i in 1:m) { # iterate through rows
    for (j in 1:m) { # iterate through cols
      if (j - i == 1) { # Set upper triangle except diagonal to lambda
        count <- count + 1 # keep track of inputing elements
        A[i, j] <- lambda[count] # assign elements
      }
    }
  }
  diag(A) <- -rowSums(A) # diag elements are equal to negative total rate of transmission
  return(A)
}
# construct.A3(m = 5, c(1,1,1,1))

# function for selecting parameterization of transition rate matrix
construct.A <- function(m, lambda, A_param = NULL) {
  if (is.null(A_param)) {
    A <- construct.A3(m, lambda) # construct transition rate matrix
  } else if (A_param == "A3") {
    A <- construct.A3(m, lambda) # construct transition rate matrix
  } else if (A_param == "A2") {
    A <- construct.A2(m, lambda) # construct transition rate matrix
  } else if (A_param == "A1") {
    A <- construct.A1(m, lambda) # construct transition rate matrix
  }
  return(A)
}


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


rps_r <- function(m, pred, obs) {
  res <- 0.0
  for (k in 1:m) {
    rps <- 0.0
    for (i in 1:k) {
      rps <- rps + (pred[i] - obs[i])
    }
    res <- res + (rps * rps)
  }
  return(res)
}

rps_vectorized <- function(pred, obs) {
  cum_diff <- cumsum(pred - obs)
  return(sum(cum_diff^2))
}

rps.score.cov <- function(m, s1, s2, u, z, beta, states, A_param = NULL) {
  log.lik <- 0 # init log likelihood
  Cz <- exp(beta[(m):length(beta)] %*% t(z)) # precompute covariate contribution
  lambda <- exp(beta[1:(m - 1)]) # transform to working parameter
  for (i in 1:length(u)) { # iterate through data
    A <- construct.A(m, lambda * Cz[i], A_param) # construct trans. rate matrix
    Pt <- s1[i] == states # compute state probs
    Ptu <- Pt %*% Matrix::expm(A * u[i]) # compute state probs forward in time u
    obs <- s2[i] == states
    log.lik <- log.lik + rps_vectorized(Ptu, obs) #weight[i] * log(Ptu[which(s2[i] == states)]) # update
  }
  return(log.lik/length(u))
}

brier_score <- function(m, pred, obs) {
  # Initialize the Brier Score
  brier <- 0.0
  
  # Calculate the Brier Score
  for (i in 1:m) {
    brier <- brier + (pred[i] - obs[i])^2
  }
  
  # Return the average Brier Score
  return(brier / length(pred))
}

brier_score_vectorized <- function(pred, obs) {
  diff = (pred - obs)
  return(sum(diff^2))
}


brier.score.cov <- function(m, s1, s2, u, z, beta, states, A_param = NULL) {
  log.lik <- 0 # init log likelihood
  Cz <- exp(beta[(m):length(beta)] %*% t(z)) # precompute covariate contribution
  lambda <- exp(beta[1:(m - 1)]) # transform to working parameter
  for (i in 1:length(u)) { # iterate through data
    A <- construct.A(m, lambda * Cz[i], A_param) # construct trans. rate matrix
    Pt <- s1[i] == states # compute state probs
    Ptu <- Pt %*% Matrix::expm(A * u[i]) # compute state probs forward in time u
    obs <- s2[i] == states
    log.lik <- log.lik + brier_score_vectorized(Ptu, obs) #weight[i] * log(Ptu[which(s2[i] == states)]) # update
  }
  return(log.lik/length(u))
}


