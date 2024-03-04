#########################
### Helping functions ###
#########################

get_full_name <- function(abbreviation) {
  abbreviation <- tolower(gsub(" ", "_", abbreviation))
  
  rail_names <- list(
    asker = "Askerbanen",
    bergen = "Bergensbanen",
    brat = "Bratsbergbanen",
    brevi = "Brevikbanen",
    dovre = "Dovrebanen",
    dramm = "Drammenbanen",
    garde = "Gardermobanen",
    gjov = "Gjøvikbanen",
    hove = "Hovedbanen",
    kongs = "Kongsvingerbanen",
    mera = "Meråkerbanen",
    nord = "Nordlandsbanen",
    ofot = "Ofotbanen",
    ostfoldbanen_v = "Østfoldbanen V",
    ostfoldbanen_o = "Østfoldbanen Ø",
    rand = "Randsfjordbanen",
    roa = "Roa-Hønefossbanen",
    roro = "Rørosbanen",
    solor = "Solørbanen",
    sorl = "Sørlandsbanen",
    spikke = "Spikkestadbanen",
    vestfold = "Vestfoldbanen"
  )
  
  full_name <- rail_names[[abbreviation]]
  
  if (is.null(full_name)) {
    return("Invalid abbreviation")
  }
  
  return(full_name)
}

convert.to.num <- function(d) {
  x <- d
  x[d == "3" & !is.numeric(d)] <- 1
  x[d == "2B" & !is.numeric(d)] <- 2
  x[d == "2A" & !is.numeric(d)] <- 3
  x[d == "1" & !is.numeric(d)] <- 4
  x[d == "0" & !is.numeric(d)] <- 5
  return(as.numeric(x))
}


convert.to.num.inv <- function(dat) {
  x <- dat
  x[dat == 1] <- "3"
  x[dat == 2] <- "2B"
  x[dat == 3] <- "2A"
  x[dat == 4] <- "1"
  x[dat == 5] <- "0"
  return(x)
}


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

# function for constructing transition rate matrix with parameterized trans rates
construct.A2 <- function(m, lambda) {
  count <- 0 # init counter
  A <- matrix(0, ncol = m, nrow = m) # init A
  for (i in m:1) { # iterate through rows
    for (j in 1:m) { # iterate through cols
      if (j - i == 1) { # Set first upper pseudo diagonal to lambda
        count <- count + 1 # keep track of inputing elements
        A[i, j] <- lambda[m - count] # assign elements
      } else if (j - i > 1) { # position condition for constrained parameters
        foo <- 0 # initialize sum variable
        for (k in i:(j - 1)) { # sum left positioned elements
          foo <- foo + A[k, k + 1] # update sum
        }
        A[i, j] <- 1 / foo # assign value to conditioned parameter
      }
    }
  }
  diag(A) <- -rowSums(A) # diag elements are equal to negative total rate of transmission
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

# construct.A(m = 5, lambda = 1:5, A_param = NULL)

# compute p-values of parameters
p.values <- function(hessian, pars) {
  Fisch <- pracma::inv(hessian) # invert hessian
  se <- sqrt(diag(Fisch)) # Get square root of diagonal for standard errors
  z.score <- pars / se # Compute z-scores
  p.val <- 2 * pnorm(-abs(z.score)) # Get p-values
  return(p.val)
}


