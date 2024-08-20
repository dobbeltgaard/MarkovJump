

#Purpose: Test script for organized estimation with RPS, Brier, and Likelihood.
# - The user should be able to choose which score, she will use for estimation. 
# - The user should be able to choose how she want to compute 
#the transient distribution of the finite space continuous time Markov chain, 
#either by Padé, Uniformization, or eigenvalue decomposition.
# - The user should be able to choose whether or not to include covariates.


rm(list = ls()) #clear memory
load("defects_covs_base.RData"); D <- foo; rm(foo) #load dat


#########################
### Data preparation ####
#########################
states <- c(1,2,3,4,5) #define states
m <- length(states) #number of states
track <- unique(D$Track) #define investigated tracks
exo.cols <- c("speed.norm","profil.norm", "steel.norm", "invRad.norm")

idx.remove <- 
  (D$`s-` == 4 & D$`s` == 4 & D$u > 300) | 
  (D$`s-` == 4 & D$`s` == 5 & D$u > 300) | 
  (D$`s-` == 5 & D$`s` == 5 & D$u > 300) 


d$s1 = d$`s-`
d$s2 = d$`s`
d$t = d$u/365
d$t_mgt = d$u/365 * d$MBT.norm

z = as.matrix(d[, exo.cols])



library(Rcpp)
library(RcppEigen)
#library(Matrix)
#source("funcs_discrete_loglik.R")
#source("funcs_helping.R")

sourceCpp("FUNCS_MJP_with_eigen.cpp")


beta0 <- c(c(0.1,0.2,0.3,0.4),
           #rep(0.25, (m-1)),
           #c(0.1,0.2,0.3,0.4, 0.2,0.2,0.3,0.4,0.1,9),
           rep(0.1,length(exo.cols)))  


#beta0 = c(c(-0.1429,0.5735,-1.0762,-0.1925,-0.0694,-0.6080,-0.8727,-0.2814,-0.0523),rep(0.1,length(exo.cols)))


MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", 
               covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, 
          transient_dist_method = "uniformization")

#discrete_loglik_cpp(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z )
#discrete_loglik_eigen_cpp(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z )

x = optim( par = beta0, fn = MJP_score, m = m, s1 = d$`s-`, s2 = d$s, u = d$t_mgt, z = z, generator="gerlang", 
           covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F,
           transient_dist_method = "eigenvalue_decomp", method = "BFGS", control = list(maxit = 1000)) 

x

#todo
#- make revised data frame with days divided by 365. change name of s.?
#- forecast and eval functions
#- speed investigation of transient distribution methods vs. generator parameterizations
#- investigate if t*MGT is more meaningfull time parameterization



#comments: 
# rps and brier estimation on gerlang relax is meaningful, it just takes more iterations to reach convergence
# sensitivity on intialization seems to be resolved
