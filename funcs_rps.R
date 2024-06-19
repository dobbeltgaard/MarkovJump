

sourceCpp("rps.cpp")


rps.estim.cov <- function(m, s1, s2, u, z, beta, states, pvalues.bin = F, method = NULL) {
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
    par = beta, fn = discrete_rps_cpp_cov, m = m,
    s1 = s1, s2 = s2, u = u, z = z,
    control = control, hessian = pvalues.bin,
    method = method) # optimize
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
      #pvalues = pvals, # return pvalues
      #CI.UP = CI.UP, # return 95 upper CI
      #CI.LO = CI.LO, # return 95 lower CI
      #mllk = -mod$value, # return log-likelihood value
      RPS = mod$value, # return RPS
      #convergence = mod$convergence, # return congernece code
      #message = mod$message, # return error messages if any
      runtime = end_time - start_time
    ) # return run time
  )
}


rps.loglik.estim.cov <- function(m, s1, s2, u, z, beta, states, pvalues.bin = F, method = NULL) {
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
    par = beta, fn = discrete_rps_loglik_cpp_cov, m = m,
    s1 = s1, s2 = s2, u = u, z = z,
    control = control, hessian = pvalues.bin,
    method = method) # optimize
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
      #pvalues = pvals, # return pvalues
      #CI.UP = CI.UP, # return 95 upper CI
      #CI.LO = CI.LO, # return 95 lower CI
      #mllk = -mod$value, # return log-likelihood value
      RPS = mod$value, # return RPS
      #convergence = mod$convergence, # return congernece code
      #message = mod$message, # return error messages if any
      runtime = end_time - start_time
    ) # return run time
  )
}



rps.estim <- function(m, s1, s2, u, beta, states, pvalues.bin = F, method = NULL) {
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
    par = beta, fn = discrete_rps_cpp, m = m,
    s1 = s1, s2 = s2, u = u,
    control = control, hessian = pvalues.bin,
    method = method) # optimize
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
      #pvalues = pvals, # return pvalues
      #CI.UP = CI.UP, # return 95 upper CI
      #CI.LO = CI.LO, # return 95 lower CI
      #mllk = -mod$value, # return log-likelihood value
      RPS = mod$value, # return RPS
      #convergence = mod$convergence, # return congernece code
      #message = mod$message, # return error messages if any
      runtime = end_time - start_time
    ) # return run time
  )
}