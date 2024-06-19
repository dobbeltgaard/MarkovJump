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


# very efficient incomplete data log-likelihood computation (all vectorized in case of diagonalizable transition rate matrix)
log.lik.cov2 <- function(m, s1, s2, u, z, beta, states, A_param = NULL, weight = NULL) {
  n <- length(u) # number of observations
  log.lik <- 0 # init log likelihood
  Cz <- exp(beta[(m):length(beta)] %*% t(z)) # precompute covariate contribution
  lambda <- exp(beta[1:(m - 1)]) # transform to working parameter
  if (is.null(weight)) {
    weight <- rep(1, length(u))
  } # if no weights, weigh all equal
  Pt <- outer(s1, states, "==") # compute state probs for all defect obs.
  Ptu <- matrix(NA, nrow = n, ncol = m) # initialize conditioned state distribution after u time
  Ptu.obs <- outer(s2, states, "==") # compute state probs for all defect obs.
  if (length(lambda) == length(unique(lambda))) { # if trans. rates are distinct, A is diagonalizable
    U <- eigenspace.U(m, lambda) # compute eigenspace
    # Uinv <- pracma::inv(U); Uinv[lower.tri(Uinv)] <- F #closed-form expression exists but isn't implemented yet
    Uinv <- eigenspace.U.inv(m, lambda)
    dia <- matrix(outer(Cz * u, c(lambda, 0), "*"), ncol = m) # update diagonal with exponential covariates
    dia <- exp(-dia) # exponential to eigenvalues (i.e. expm(Lambda))
    tpm <- array(0, dim = c(m, m, n)) # initialize 3d transition prob. matrix
    for (i in 1:(m)) { # for each distinct eigenvalue, compute matrix product
      U.foo <- matrix(0, nrow = m, ncol = m) # initialize col space matrix product
      U.foo[, i] <- U[, i] # asign i'th column space from eigenspace
      U.Uinv.foo <- U.foo %*% Uinv # precompute matrix product for chosen column space
      tpm <- tpm + outer(U.Uinv.foo, dia[, i], FUN = "*") # scale by eigenvalue * constant and update tpm
    }
    trows <- which(t(Pt)) - c(0, seq(m, m * n - 1, by = m)) # find which columns is true for all rows in Pt
    trows <- rep(trows, each = m) # stack in long vector
    all.rows <- rep(1:m, n) # choose all rows
    all.mats <- rep(1:n, each = m) # stack in long vector
    index.mat <- cbind(trows, all.rows, all.mats) # construct matrix indicating end points, total corresponding row space of tpm, for all tpm's
    Ptu <- matrix(tpm[index.mat], ncol = m, nrow = n, byrow = T) # extract probability of observed endpoint
    log.lik <- sum(weight * log(Ptu[Ptu.obs] + 10^-10)) # sum probs. of observed end points
  } else { # if A not diagonalizable, then:
    for (i in 1:n) { # iterate through data
      A <- construct.A(m, lambda * Cz[i], A_param) # construct trans. rate matrix
      Ptu[i, ] <- as.matrix(Pt[i, ] %*% Matrix::expm(A * u[i])) # compute state probs forward in time u
    }
    log.lik <- sum(weight * log(Ptu[Ptu.obs])) # sum probs. of observed end points
  }
  return(-log.lik) # return negative log likelihood
}


sourceCpp("funcs_discrete_loglik.cpp")
log.lik.cov3 <- function(m, s1, s2, u, z, beta){
  log.lik <- 0
  if (length(beta) == length(unique(beta))) {
    log.lik <- discrete_loglik_eigen_cpp(m = m, s1 = s1, s2 = s2, u = u, z = z, pars = beta)
  } else {
    log.lik <- discrete_loglik_cpp(m = m, s1 = s1, s2 = s2, u = u, z = z, pars = beta)
  }
  return(log.lik)
}


log.lik.grad.cov3 <- function(m, s1, s2, u, z, beta){
  if (length(beta) == length(unique(beta))) { #
    grad <- discrete_loglik_eigen_grad_cpp(m = m, s1 = s1, s2 = s2, u = u, z = z, pars = beta)
  } else {
    grad <- pracma::grad(f = discrete_loglik_cpp, x0 = beta, m = m, s1 = s1, s2 = s2, u = u, z = z)
  }
  return(grad)
}


# maximum likelihood estimation
mle.cov <- function(m, s1, s2, u, z, beta, states,
                    A_param = NULL, pvalues.bin = F, method = NULL, weight = NULL) {
  start_time <- Sys.time() # Start of run time
  if (is.null(method) | method == "Nelder-Mead") { # Nelder-Mead = standard method
    method <- "Nelder-Mead"
    control <- list(maxit = 20000)
  }
  if (method == "SANN") { # Simulated annealing
    control <- list(maxit = 20000, temp = 10, tmax = 10)
  }
  if (method == "BFGS") { # Simulated annealing
    control <- list()
  }
  mod <- optim(
    par = beta, fn = log.lik.cov3, m = m,
    s1 = s1, s2 = s2, u = u, z = z,
    control = control, hessian = pvalues.bin,
    method = method
  ) # optimize
  end_time <- Sys.time() # End of run time
  if (pvalues.bin) { # check if pvalues should be computed
    hessian <- mod$hessian # retrieve hessian
    pvals <- round(p.values(hessian, mod$par), 5) # compute values
    CI.UP <- mod$par + (1.96 / sqrt(length(u) * diag(pracma::inv(hessian))))
    CI.LO <- mod$par - (1.96 / sqrt(length(u) * diag(pracma::inv(hessian))))
  } else {
    pvals <- NULL
    hessian <- NULL
    CI.UP <- NULL
    CI.LO <- NULL
  }
  return(
    list(
      lambda = mod$par, # retur working params
      pvalues = pvals, # return pvalues
      CI.UP = CI.UP, # return 95 upper CI
      CI.LO = CI.LO, # return 95 lower CI
      mllk = -mod$value, # return log-likelihood value
      AIC = 2 * (mod$value + length(beta)), # return AIC
      convergence = mod$convergence, # return congernece code
      message = mod$message, # return error messages if any
      runtime = end_time - start_time
    ) # return run time
  )
}


mle.no.covs <- function(m, s1, s2, u, beta, states,
                    A_param = NULL, pvalues.bin = F, method = NULL, weight = NULL) {
  start_time <- Sys.time() # Start of run time
  if (is.null(method) | method == "Nelder-Mead") { # Nelder-Mead = standard method
    method <- "Nelder-Mead"
    control <- list(maxit = 20000)
  }
  if (method == "SANN") { # Simulated annealing
    control <- list(maxit = 20000, temp = 10, tmax = 10)
  }
  if (method == "BFGS") { # Simulated annealing
    control <- list()
  }
  mod <- optim(
    par = beta, fn = discrete_loglik_cpp_nocovs, m = m,
    s1 = s1, s2 = s2, u = u,
    control = control, hessian = pvalues.bin,
    method = method
  )
  end_time <- Sys.time() # End of run time
  if (pvalues.bin) { # check if pvalues should be computed
    hessian <- mod$hessian # retrieve hessian
    pvals <- round(p.values(hessian, mod$par), 5) # compute values
    CI.UP <- mod$par + (1.96 / sqrt(length(u) * diag(pracma::inv(hessian))))
    CI.LO <- mod$par - (1.96 / sqrt(length(u) * diag(pracma::inv(hessian))))
  } else {
    pvals <- NULL
    hessian <- NULL
    CI.UP <- NULL
    CI.LO <- NULL
  }
  return(
    list(
      lambda = mod$par, # retur working params
      pvalues = pvals, # return pvalues
      CI.UP = CI.UP, # return 95 upper CI
      CI.LO = CI.LO, # return 95 lower CI
      mllk = -mod$value, # return log-likelihood value
      AIC = 2 * (mod$value + length(beta)), # return AIC
      convergence = mod$convergence, # return congernece code
      message = mod$message, # return error messages if any
      runtime = end_time - start_time
    ) # return run time
  )
}

# CI <- res2$lambda + (1.96 / sqrt( nrow(d)*diag(pracma::inv(res2$hessian)) ) )
# CI2 <- res2$lambda - (1.96 / sqrt( nrow(d)*diag(pracma::inv(res2$hessian)) ) )



gradient.log.lik.cov2 <- function(m, s1, s2, u, z, beta, states, A_param = NULL, weight = NULL){
  gr1 <- rep(0, m-1) #initialize
  gr2 <- rep(0, NCOL(z)) #initialize
  n <- NROW(z) #number of observations
  lambda <- exp(beta[1:(m-1)]) #transform to working parameter
  A <- construct.A(m, lambda,A_param) #construct trans. rate matrix
  covs.effect <- z %*% beta[m:(length(beta))] #precompute for efficiency
  exp.covs.effect <- exp(covs.effect) #precompute for efficiency
  sum.exp.covs.effect <- sum(exp.covs.effect) #precompute for efficiency
  exp.base <- exp(beta[1:(m-1)]) #precompute for efficiency
  sum.z <- colSums(z) #precompute for efficiency
  Cz <- exp(beta[(m):length(beta)] %*% t(z)) #precompute covariate contribution
  lambda <- exp(beta[1:(m-1)]) #transform to working parameter
  if(is.null(weight)){weight <- rep(1, length(u))} #if no weights, weigh all equal
  Pt <- outer(s1, states, "==") #compute state probs for all defect obs.
  Ptu <- matrix(NA, nrow = n, ncol = m) #initialize conditioned state distribution after u time
  Ptu.obs <- outer(s2, states, "==") #compute state probs for all defect obs.
  if( length(lambda) == length(unique(lambda))){ #if trans. rates are distinct, A is diagonalizable
    U <- eigenspace.U(m, lambda) #compute eigenspace
    #Uinv <- pracma::inv(U); Uinv[lower.tri(Uinv)] <- F #closed-form expression exists but isn't implemented yet
    Uinv <- eigenspace.U.inv(m, lambda)
    dia <- matrix(outer(Cz*u, c(lambda,0), "*"), ncol = m) #update diagonal with exponential covariates
    dia <- exp(-dia) #exponential to eigenvalues * k (i.e. expm(Lambda * k))
    #tpm <- array(0,dim = c(m,m,n)) #initialize 3d transition prob. matrix
    #for(i in 1:(m)){ #for each distinct eigenvalue, compute matrix product
    #  U.foo <- matrix(0, nrow = m, ncol = m) #initialize col space matrix product
    #  U.foo[,i] <- U[,i] #asign i'th column space from eigenspace
    #  U.Uinv.foo <- U.foo %*% Uinv #precompute matrix product for chosen column space
    #  tpm <- tpm + outer(U.Uinv.foo, dia[,i],FUN = "*") #scale by eigenvalue * constant and update tpm
    #}
    
    for(i in 1:length(u)){
      P <- diag(dia[i,]) #exp(Lambda*k)
      factor <- Cz[i]*u[i] 
      tpm2 <- (U %*% P %*% Uinv) #transition probability matrix
      
      for(j in 1:(m-1)){
        U.grad <- eigenspace.U.grad(m, lambda, j)
        U.inv.grad <- eigenspace.U.inv.grad(m, lambda, j)
        E <- matrix(F, ncol = m, nrow = m); E[j,j] <- T;
        #E <- diag(rep(T,m)); #E[j,j] <- F;
        P.grad <- factor* P %*% E   #transition probability matrix
        tpm <- (U.grad %*% P %*% Uinv) + (U %*% P.grad %*% Uinv) + (U %*% P %*% U.inv.grad)
        
        gr1[j] <- gr1[j] + ( (Pt[i,] %*% tpm)[Ptu.obs[i,]]) / ((Pt[i,] %*% tpm2)[Ptu.obs[i,]])
        
      }
      for(j in 1:NCOL(z)){ #iterate through covariate coefs.
        tpm3 <- z[i,j]*factor* A %*% tpm2
        gr2[j] <- gr2[j] + ((Pt[i,] %*% tpm3)[Ptu.obs[i,]]) / ( (Pt[i,] %*% tpm2)[Ptu.obs[i,]])   #update gradient
      }
    }
    return(-c(gr1, gr2)) #return negative gradient
    #trows <- which(t(Pt))- c(0,seq(m,m*n-1, by = m)) #find which columns is true for all rows in Pt
    #trows <- rep(trows, each = m) #stack in long vector
    #all.rows <- rep(1:m, n) #choose all rows
    #all.mats <- rep(1:n, each = m) #stack in long vector
    #index.mat <- cbind(trows, all.rows, all.mats) #construct matrix indicating end points, total corresponding row space of tpm, for all tpm's
    #Ptu <- matrix(tpm[index.mat], ncol = m, nrow = n, byrow = T) #extract probability of observed endpoint
    #log.lik <- sum( weight * log(Ptu[Ptu.obs] + 10^-12) ) #sum probs. of observed end points
  } else { #if A not diagonalizable, then:
    return(pracma::grad(f = log.lik.cov2, x0 = beta, m = m, s1 = s1, s2 = s2, u = u, states = states,  z = z))
  }
}
























