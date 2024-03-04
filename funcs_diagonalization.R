############################
##### Diagonalization ######
############################
eigenspace.U <- function(m, lambda, i.idx = NULL, j.idx = NULL) {
  # lambda <- exp(beta)
  U <- matrix(1, nrow = m, ncol = m)
  U[lower.tri(U)] <- 0
  for (i in 1:(m - 1)) {
    for (j in 1:(m - 1)) {
      if (j > i) {
        for (k in i:(j - 1)) {
          U[i, j] <- U[i, j] * (lambda[k] / (lambda[k] - lambda[j]))
          # U[i,j] <- U[i,j] * (lambda[k] / (lambda[j] * (lambda[k] / lambda[j] - 1)))
        }
      }
    }
  }
  if (is.null(i.idx) | is.null(j.idx)) {
    return(U)
  } else {
    return(U[i, j])
  }
}
# eigenspace.U(m, c(1,3,4,2.5) )


eigenspace.U.grad <- function(m, lambda, h) {
  U <- matrix(1, nrow = m, ncol = m)
  Us <- matrix(0, nrow = m, ncol = m)
  for (i in 1:(m - 1)) {
    for (j in 1:(m - 1)) {
      if (j > i) {
        for (k in i:(j - 1)) {
          if (h == j) {
            # Us[i,j] <- Us[i,j] + fi.eigenspace.dev.frac2(beta[k], beta[j])
            Us[i, j] <- Us[i, j] + 1 / (lambda[k] - lambda[j])
          }
          if (h == k) {
            # Us[i,j] <- Us[i,j] + fi.eigenspace.dev.frac1(beta[k], beta[j])
            Us[i, j] <- Us[i, j] + lambda[j] / (lambda[k] * (-lambda[k] + lambda[j]))
          }
          # U[i,j] <- U[i,j] * fi.eigenspace(beta[k], beta[j])
          U[i, j] <- U[i, j] * lambda[k] / (lambda[k] - lambda[j])
        }
      }
    }
  }
  return(U * Us)
}

# eigenspace.U.grad(m, (c(1,3,4,2.5) ), 3)


eigenspace.U.inv <- function(m, lambda) {
  U <- matrix(1, nrow = m, ncol = m)
  U[lower.tri(U)] <- 0
  for (i in 1:(m - 1)) {
    for (j in 1:(m - 1)) {
      if (j > i) {
        for (k in i:(j - 1)) {
          U[i, j] <- U[i, j] * (lambda[k] / (lambda[i] - lambda[k + 1]))
        }
        U[i, j] <- U[i, j] * (-1)^(j - i)
      }
    }
    if (!any(diff((i + 1):(m - 1)) < 0)) { # non-empty products
      for (k in (i + 1):(m - 1)) { # last column, j = m
        U[i, m] <- U[i, m] * (lambda[k] / (lambda[i] - lambda[k]))
      }
    }
    U[i, m] <- U[i, m] * (-1)^(m - i)
  }
  return(U)
}
# eigenspace.U.inv(m, (c(1,3,4,2.5) ))



eigenspace.U.inv.grad <- function(m, lambda, h) {
  U <- matrix(1, nrow = m, ncol = m)
  Us <- matrix(0, nrow = m, ncol = m)
  for (i in 1:(m - 1)) {
    for (j in 1:(m - 1)) {
      if (j > i) {
        for (k in i:(j - 1)) {
          if (h == i) {
            Us[i, j] <- Us[i, j] - 1 / (lambda[i] - lambda[k + 1])
          }
          if (h == k) {
            Us[i, j] <- Us[i, j] + 1 / lambda[k]
          }
          if (h == (k + 1)) {
            Us[i, j] <- Us[i, j] + 1 / (lambda[i] - lambda[k + 1])
          }
          U[i, j] <- U[i, j] * (lambda[k] / (lambda[i] - lambda[k + 1]))
        }
        U[i, j] <- U[i, j] * (-1)^(j - i)
      }
    }
    if (!any(diff((i + 1):(m - 1)) < 0)) { # non-empty products
      for (k in (i + 1):(m - 1)) { # last column, j = m
        if (h == i) {
          Us[i, m] <- Us[i, m] - 1 / (lambda[i] - lambda[k])
        }
        if (h == k) {
          Us[i, m] <- Us[i, m] + lambda[i] / (lambda[k] * (lambda[i] - lambda[k]))
        }
        U[i, m] <- U[i, m] * (lambda[k] / (lambda[i] - lambda[k]))
      }
    }
    U[i, m] <- U[i, m] * (-1)^(m - i)
  }
  return(U * Us)
}
# eigenspace.U.inv.grad(m, (c(1,3,4,2.5)), 3)

