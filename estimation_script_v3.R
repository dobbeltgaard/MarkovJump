

#Purpose: Test script for organized estimation with RPS, Brier, and Likelihood.
# - The user should be able to choose which score, she will use for estimation. 
# - The user should be able to choose how she want to compute 
#the transient distribution of the finite space continuous time Markov chain, 
#either by Padé, Uniformization, or eigenvalue decomposition.
# - The user should be able to choose whether or not to include covariates.


rm(list = ls()) #clear memory
d = read.csv(file = "defect_data.csv") #read data
states <- c(1,2,3,4,5)
m = length(states)

exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
z = as.matrix(d[,exo.cols]);


library(Rcpp)
library(RcppArmadillo)
library(Matrix)


source("funcs_discrete_loglik.R")
source("funcs_helping.R")

sourceCpp("FUNCS_MJP.cpp")


beta0 <- c(c(0.1,0.2,0.3,0.4),
           #rep(0.25, (m-1)), 
           rep(0.1,length(exo.cols)))  

log.lik.cov(m = m, s1 = d$`s.`, s2 = d$s, u = d$u/365, z = z, beta = beta0, states)
MJP_score(m = m, s1 = d$`s.`, s2 = d$s, u = d$u/365,pars = beta0, z = z, generator="erlang", 
               covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, 
          transient_dist_method = "uniformization")




#eigenvalue_decomp

#todo
#- see if it works for estimation
#- make revised data frame with days divided by 365. change name of s.?
#- eigenvalue decomp of relaxed g_erlang. 
#- forecast and eval functions
#- speed investigation of transient distribution methods



