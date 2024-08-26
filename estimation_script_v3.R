

#Purpose: Test script for organized estimation with RPS, Brier, and Likelihood.
# - The user should be able to choose which score, she will use for estimation. 
# - The user should be able to choose how she want to compute 
#the transient distribution of the finite space continuous time Markov chain, 
#either by Padé, Uniformization, or eigenvalue decomposition.
# - The user should be able to choose whether or not to include covariates.


rm(list = ls()) #clear memory
d = read.csv("defect_data.csv")

#########################
### Data preparation ####
#########################
states <- c(1,2,3,4,5) #define states
m <- length(states) #number of states
track <- unique(d$Track) #define investigated tracks
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
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

head(d)
MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", 
               covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, 
          transient_dist_method = "uniformization")
x = 
  optim( par = beta0, fn = MJP_score, m = m, s1 = d$s1, s2 = d$s2, u = d$t, z = z, 
         generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, 
         transient_dist_method = "eigenvalue_decomp", 
         method = "BFGS", control = list(maxit = 1000)) 
pred = 
  MJP_predict(m = m, s1 = d$s1, u = d$t, pars = x$par, z = z, generator = "gerlang",
              covs_bin = T, transient_dist_method = "eigenvalue_decomp")
obs = make_Ptu_obs(m, d$s2)
mean(logscore_vectors(m, pred, obs))





#EMPIRICAL DISTRIBUTION PREDICTORS
naive.empirical.pred = function(x){return(as.vector(table(x)/length(x)))}
corr.empirical.pred = function(m, idx, naive){
  if(idx > 1){
    prob = sum(naive[1:(idx-1)])
    naive[1:(idx-1)] = 0
    naive[idx:m] = naive[idx:m] + prob/(m-idx+1) 
  } 
  return(naive)
}



#INITIALIZE
set.seed(123)
k <- 10
ints = seq(1, NROW(d), ceiling(NROW(d)/k)-1)
ints[length(ints)] = NROW(d)
rand.idx <- sample(x = 1:NROW(d),size = NROW(d),replace = F)

beta00 <- rep(0.25, (m-1)) 
beta0 <- c(rep(0.25, (m-1)), rep(0.1,length(exo.cols))) 
beta02 <- c(rep(0.25, (m-1)), rep(0.1,length(exo.cols2))) 


PARS = list() 




#todo
#- investigate if t*MGT is more meaningfull time parameterization



#comments: 
# rps and brier estimation on gerlang relax is meaningful, it just takes more iterations to reach convergence
# sensitivity on intialization seems to be resolved
