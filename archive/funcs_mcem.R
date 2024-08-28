

# function to compute transition count matrix
count.transitions <- function(m, vec) {
  trans <- table(vec[-length(vec)], vec[-1]) # counting
  Tab <- as.table(matrix(0, nrow = m, ncol = m)) # empty large matrix
  colnames(Tab) <- 1:m
  rownames(Tab) <- 1:m # set rownames and colnames
  
  # Find the row and column indices for embedding
  row_indices <- match(rownames(trans), rownames(Tab))
  col_indices <- match(colnames(trans), colnames(Tab))
  
  # Embed the smaller table into the bigger table
  Tab[row_indices, col_indices] <- trans
  return(matrix(Tab, ncol = m, nrow = m))
}



# function to compute holding time for all states
holding.time <- function(m, paths, horizon) {
  R <- rep(0, m) # init. holding values
  for (j in 1:m) { # iterate through states
    R[j] <- sum(paths == j) # sum over all indeces for which the state is j
  }
  res <- horizon / (NROW(paths) - 1) # compute resolution in yrs
  return(R * res) # return scaled holding times
}

sim.Markov.cov <- function(m, states, beta, A_param = NULL, s1, u, z, n.samples, step.size) {
  Cz <- exp(sum(beta[(m):length(beta)] * (z))) # precompute covariate contribution
  lambda <- exp(beta[1:(m - 1)]) # transform to working parameter
  # if(is.null(weight)){weight <- rep(1, length(u))} #if no weights, weigh all equal
  A <- construct.A(m = m, lambda = lambda * Cz, A_param = A_param) # construct transition matrix
  
  Time <- seq(0, u, by = step.size)
  PATHS <- matrix(NA, nrow = length(Time), ncol = n.samples) # init Markov trajectories
  PATHS[1, ] <- s1 # save initial
  
  P <- Matrix::expm(A * step.size) # compute probability matrix
  obs.state <- s1 == states # find state distribution
  dist <- obs.state %*% P # compute distribution at next time step
  s <- sample(x = states, size = n.samples, replace = T, prob = as.matrix(dist)) # sample from the state distribution
  
  for (i in 2:length(Time)) { # iterate through time
    obs.state <- outer(s, states, `==`) # compute path for total sample, i.e. compute state vectors for new sample
    dist <- as.matrix(obs.state %*% P) # state dependent distributions for sample, (carried forward to next time step)
    s <- apply(dist, 1, function(prob_row) {
      sample(x = states, size = 1, replace = TRUE, prob = prob_row)
    }) # sample from state-dependent distribution
    PATHS[i, ] <- s # store sample
  }
  return(PATHS)
}


# function to simulate Markov trajectories removing
sim.Markov.cov_from_i_to_j <- function(m, states, beta, A_param = NULL, 
                                       s1, s2, u, z, n.samples, step.size) {
  Cz <- exp(sum(beta[(m):length(beta)] * (z))) #  covariate contribution
  lambda <- exp(beta[1:(m - 1)]) # transform to working parameter
  # if(is.null(weight)){weight <- rep(1, length(u))} #if no weights
  A <- construct.A(m = m, lambda = lambda * Cz, 
                   A_param = A_param) # construct transition matrix
  
  Time <- seq(0, u, by = step.size)
  
  if (s1 == s2) { # if states are equal then it simplifies
    PATHS <- matrix(s1, nrow = length(Time), ncol = n.samples) 
  } else {
    P <- Matrix::expm(A * step.size) # compute probability matrix
    
    for (i in 1:m) { # find a better way to simplify the matrix exponential
      for (j in 1:m) {
        if (j - i > 1) P[i, j] <- 0
      }
    }
    
    max_time_seconds <- 60 # Set the maximum allowed time in seconds
    start_time <- Sys.time() # Capture the start time
    found <- FALSE
    while (!found) {
      PATHS <- matrix(NA, nrow = length(Time), ncol = n.samples) # init traject.
      PATHS[1, ] <- s1 # save initial
      
      obs.state <- s1 == states # find state distribution
      dist <- obs.state %*% P # compute distribution at next time step
      s <- sample(x = states, size = n.samples, replace = T, 
                  prob = as.matrix(dist)) # sample from the state distribution
      
      for (i in 2:length(Time)) { # iterate through time
        obs.state <- outer(s, states, `==`) # compute path for total sample
        dist <- as.matrix(obs.state %*% P) # state dependent distributions 
        s <- apply(dist, 1, function(prob_row) {
          sample(x = states, size = 1, replace = TRUE, prob = prob_row)
        }) # sample from state-dependent distribution
        if (sum(s <= s2) == 0) { # check if no samples could hit end point
          break # if so, break out of loop...
        } else { # if not, then ...
          PATHS <- PATHS[, s <= s2, drop = FALSE] # remove irrelevant samples 
          s <- s[s <= s2] # same...
          PATHS[i, ] <- s # store samples with potential of hitting end point
        }
        
        elapsed_time <- as.numeric(
          difftime(Sys.time(), start_time, units = "secs")) # Check elapsed time
        # Terminate the loop if it runs for more than a minute
        if (elapsed_time > max_time_seconds) {
          print(paste0("Terminating while loop in obs. = ", i))
          return(NULL)
          # break
        }
      }
      # print(PATHS[NROW(PATHS), , drop = FALSE])
      path.ends <- PATHS[NROW(PATHS), , drop = FALSE] == s2 # check paths 
      found <- sum(path.ends) > 0 # termine while loop if any found
      if (is.na(found)) found <- FALSE # if no proper path, then redo
      PATHS <- PATHS[, path.ends, drop = FALSE]
    }
  }
  return(PATHS) # return only samples that hit the end point
}


MCEM <- function(m, s1, s2, u, z, beta0, states,
                 iters = NULL, step.size = NULL, n.samples = NULL) {
  # Default options
  if (is.null(iters)) {
    iters <- 10
  }
  if (is.null(step.size)) {
    step.size <- 1 / 365 * 5
  }
  if (is.null(n.samples)) {
    n.samples <- 10
  }
  
  # Initialization
  ESTIM <- matrix(NA, ncol = m - 1 + NCOL(z), nrow = iters + 1)
  ESTIM[1, ] <- beta0
  N <- list()
  R <- list()
  
  start_time <- Sys.time()
  bin.foo <- F
  for (j in 1:iters) { # iterate by chosen number of iterations
    print(beta0) # current estimate
    idx.foo <- rep(T, NROW(u)) # init well defined rows from rejection sample
    Nrel <- matrix(NA, ncol = m - 1, nrow = NROW(z)) # Nrel for better storage
    Rrel <- matrix(NA, ncol = m, nrow = NROW(z)) # Rrel for better storage
    for (i in 1:NROW(u)) { # iterate through data points
      ## Step 2 ##
      paths <-
        sim.Markov.cov_from_i_to_j(
          m = m, states = states, beta = beta0, A_param = NULL,
          s1 = s1[i], s2 = s2[i], u = u[i], z = z[i, ],
          n.samples = n.samples, step.size = step.size
        ) # make sample paths over time
      if (is.null(paths)) {
        idx.foo[i] <- F
        next
      } else {
        max.jump <- apply(diff(paths), 2, max) # find maximum jump of all paths
        
        ## Step 3: Estimate N and R ##
        Nrel[i, ] <- t(count.transitions(
          m = m, vec = paths[, which.min(max.jump)]))[
            c((1:(m - 1) + 1) + seq(0, (m - 2) * m, m))] # store upper diagonal
        Rrel[i, ] <- holding.time(m = m, paths[, which.min(max.jump)], 
                                  horizon = u[i])
      }
    }
    if (j == iters) {
      bin.foo <- T
    } # if last iteration, compute Hessian
    ## Step 4: Update model parameters ##
    mod <- optim(
      par = beta0, fn = full.log.lik.cov3, gr = gradient.full.log.lik.cov3,
      m = m, N = Nrel[idx.foo, ], R = Rrel[idx.foo, ], z = z[idx.foo, ],
      method = "BFGS", hessian = bin.foo
    ) # optimization
    beta0 <- mod$par # retrieve updated parameters
    ESTIM[j + 1, ] <- beta0 # storing solution
    
    ## Step 5: Goto step 2 ##
  }
  end_time <- Sys.time() # End of run time
  hessian <- mod$hessian # retrieve hessian
  pvals <- round(p.values(hessian, mod$par), 5) # compute values
  CI.UP <- mod$par + (1.96 / sqrt(length(u) * diag(pracma::inv(hessian))))
  CI.LO <- mod$par - (1.96 / sqrt(length(u) * diag(pracma::inv(hessian))))
  return(
    list(
      beta = beta0, # retur working params
      beta_iter = ESTIM,
      pvalues = pvals, # return pvalues
      CI.UP = CI.UP, # return 95 upper CI
      CI.LO = CI.LO, # return 95 lower CI
      runtime = end_time - start_time
    ) # return run time
  )
}

# very efficient full data likelihood computation (all vectorized)
full.log.lik.cov3 <- function(m, N, R, z, pars, A_param = NULL, weight = NULL) {
  covs.effect <- z %*% pars[m:(length(pars))] # coefficients * covariates
  exp.covs.effect <- exp(covs.effect) # exponential of covariate effect
  exp.base <- exp(pars[1:(m - 1)]) # base transition rates
  pars.covs.effect <- t(matrix(outer(pars[1:(m - 1)], covs.effect, "+"), nrow = m - 1)) # base coefficients added to covariate effect for all data points
  exp.covs.effect.R <- colSums(as.vector(exp.covs.effect) * R[, 1:(m - 1)]) # exp(covariate effect) multiplied with holding time
  log.lik <- sum(N * pars.covs.effect) - sum(exp.base * exp.covs.effect.R) # log.likelihood computation
  return(-log.lik) # return negative log likehood
}

gradient.full.log.lik.cov3 <- function(m, N, R, z, pars, A_param = NULL, weight = NULL) {
  gr1 <- rep(0, m - 1) # initialize
  gr2 <- rep(0, NCOL(z)) # initialize
  covs.effect <- z %*% pars[m:(length(pars))] # precompute for efficiency
  exp.covs.effect <- exp(covs.effect) # precompute for efficiency
  exp.base <- exp(pars[1:(m - 1)]) # precompute for efficiency
  exp.covs.effect.R <- as.vector(exp.covs.effect) * R[, 1:(m - 1)] # exp(covariate effect) multiplied with holding time
  z.exp.covs.effect <- z * as.vector(exp.covs.effect)
  exp.base.R <- sweep(R[, 1:(m - 1)], MARGIN = 2, exp.base, "*")
  gr1 <- colSums(N)/exp.base - colSums(exp.covs.effect.R)
  gr2 <- colSums(z * rowSums(N))
  gr2 <- gr2 - colSums(z.exp.covs.effect * rowSums(exp.base.R))
  return(-c(gr1, gr2)) # return negative gradient
}