rm(list = ls())
graphics.off()


library(Rcpp)
library(RcppEigen)
library(TMB)

compile("objective_function_erlang.cpp")
dyn.load(dynlib("objective_function_erlang"))

D <- read.csv("defect_data.csv")

states <- c("3", "2B", "2A", "1", "0")

m <- length(states)

params_n <- m-1

data <- list(
  s1 = D$s1,
  s2 = D$s2,
   t =  D$t
)

params_0 <- runif(params_n, 0, 1)

parameters <- list(
  Av = params_0
)

l <- MakeADFun(data = data, parameters = parameters, DLL = "objective_function_erlang")

l$fn(params_0)
l$gr(params_0)

# Fit 1
fit1 <- optim(par=params_0, fn=l$fn, gr=l$gr, method="BFGS")
fit1

# Fit 2
fit2 <- nlminb(l$par, l$fn, l$gr, l$he)
fit2


#############################
rm(list = ls())
library(stats)
library(Rcpp)
library(RcppEigen)
library(TMB)
compile("FUNCS_MJP_with_TMB.cpp")
dyn.load(dynlib("FUNCS_MJP_with_TMB"))

states <- c("3", "2B", "2A", "1", "0")
m <- length(states)
D <- read.csv("defect_data.csv")
data <- list(
  s1 = D$s1,
  s2 = D$s2,
  u = D$t,
  z = as.matrix(D[, c("MBT.norm", "profil.norm", "speed.norm", "steel.norm", "invRad.norm")]),
  m = m,
  generator_type = 2,  
  cov_type = 1
)

params_n <- (m - 1)*m/2 + NCOL(data$z)
#params_n <- (m - 1) 
params_0 <- runif(params_n, 0, 1)
parameters <- list(theta = params_0)
l <- MakeADFun(data = data, parameters = parameters, DLL = "FUNCS_MJP_with_TMB")

l$fn(params_0)
l$gr(params_0)

# Fit 1
fit1 <- optim(par=params_0, fn=l$fn, gr=l$gr, method="BFGS")
fit1

# Fit 2
fit2 <- nlminb(l$par, l$fn, l$gr, l$he)
fit2

##remember to use softmax and not exp()
