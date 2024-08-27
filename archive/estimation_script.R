rm(list = ls()) #clear memory

library(Rcpp);library(RcppEigen)
source("archive/funcs_diagonalization.R")
source("archive/funcs_discrete_loglik.R")
source("archive/funcs_forecasting.R")
source("archive/funcs_helping.R")
source("archive/funcs_mcem.R")


load("defects_covs_base.RData"); D <- foo; rm(foo) #load dat

#########################
### Data preparation ####
#########################
states <- c(1,2,3,4,5) #define states
m <- length(states) #number of states
track <- unique(D$Track) #define investigated tracks
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", 
              "invRad.norm")#, , "dist.to.station.trans") #covariate space

idx.remove <- 
  (D$`s-` == 4 & D$`s` == 4 & D$u > 300) | 
  (D$`s-` == 4 & D$`s` == 5 & D$u > 300) | 
  (D$`s-` == 5 & D$`s` == 5 & D$u > 300) 


d <- D[!idx.remove, c("pos", "s-", "s", "u", exo.cols, "Track0")] #define data set of interest




#########################
### DIRECT MAXIMATION ###
#########################
beta0 <- c(rep(0.25, (m-1)), rep(0.1,length(exo.cols)))  
res3 <- mle.cov(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, 
                z = as.matrix(d[,exo.cols]), beta = beta0, states = states, 
                A_param = NULL, pvalues.bin =  T, method = "BFGS", weight = NULL)




#########################
## Monte Carlo EM algo ##
#########################
beta0 <- c(rep(0.25, (m-1)), rep(0.1,length(exo.cols)))  
res4 <- MCEM(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = as.matrix(d[,exo.cols]), 
             beta0 = beta0, states = states, iters = 10, step.size = NULL, n.samples = NULL )

