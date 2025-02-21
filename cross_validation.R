rm(list = ls()) #clear memory

### Read data and define state space and covariates ###
d = read.csv("defect_data.csv")
states <- c(1,2,3,4,5)
m <- length(states) 
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")


library(MASS)
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")

set.seed(123)
k <- 10
ints = seq(1, NROW(d), ceiling(NROW(d)/k)-1)
ints[length(ints)] = NROW(d)
rand.idx <- sample(x = 1:NROW(d),size = NROW(d),replace = F)
PARS = list()


source("archive/funcs_helping.R")
source("archive/funcs_mcem.R")


#RUN K-FOLD PREDICTION
for(i in 1:k){
  start_time <- Sys.time()
  start = ints[i]; 
  end = ints[i+1]; 
  
  #Estimation and test set splits
  pred.idx = rand.idx[start:end]
  d.test = d[pred.idx, ]
  d.train = d[-pred.idx, ]
  
  beta_base = rep(0.25, (m-1))
  beta = c(beta_base, rep(0.1,length(exo.cols)))
  nam = "direct"
  foo = optim( par = beta, fn = MJP_score, m = m, s1 = d.train$s1, s2 = d.train$s2, u = d.train$t, z = as.matrix(d.train[,exo.cols]), generator = "gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "eigenvalue_decomp", method = "BFGS", control = list(maxit = 1000)) #estimate model 
  if(k == 1){ PARS[[nam]] = foo$par} else { PARS[[nam]] = rbind(PARS[[nam]],foo$par)} #store model estimates

  nam = "mcem"
  foo <- MCEM(m = m, s1 = d.train$s1, s2 = d.train$s2, u = d.train$t, z = as.matrix(d.train[,exo.cols]), 
               beta0 = beta, states = states, iters = 3, step.size = NULL, n.samples = NULL )
  if(k == 1){ PARS[[nam]] = foo$beta} else { PARS[[nam]] = rbind(PARS[[nam]],foo$beta)} #store model estimates
  
  end_time <- Sys.time(); runtime <- end_time - start_time
  cat("Fold ", i, "runtime:", runtime, "\n") #approx 3.5 mins per fold
}


#saveRDS(PARS, file  = "results/cross_validation_estimates.RData")

## Get ranges of parameters in Cross validation procedure
rm(list = ls()) #clear memory
m = 5
P = readRDS("results/cross_validation_estimates.RData")
M = P[["direct"]][, -(1:(m-1))]
apply(M , MARGIN = 2, FUN = min)
apply(M , MARGIN = 2, FUN = max)





