rm(list = ls()) #clear memory

setwd("C:/Users/ATBD/OneDrive - Danmarks Tekniske Universitet/MarkovJump")
#setwd("C:/Users/askbi/OneDrive - Danmarks Tekniske Universitet/MarkovJump")
library(RcppEigen); library(Rcpp)

#source("funcs_diagonalization.R")
#source("funcs_discrete_loglik.R")
#source("funcs_forecasting.R")
#source("funcs_helping.R")
#source("funcs_mcem.R")


load("defects_covs_base.RData"); D <- foo; rm(foo) #load dat

#########################
### Data preparation ####
#########################
states <- c(1,2,3,4,5) #define states
m <- length(states) #number of states
track <- unique(D$Track) #define investigated tracks
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
idx.remove <- 
  (D$`s-` == 4 & D$`s` == 4 & D$u > 300) | 
  (D$`s-` == 4 & D$`s` == 5 & D$u > 300) | 
  (D$`s-` == 5 & D$`s` == 5 & D$u > 300) 


d <- D[!idx.remove, c("pos", "s-", "s", "u", exo.cols, "Track0")] #define data set of interest

d$s1 = d$`s-`
d$s2 = d$`s`
d$t = d$u/365
### COMPARISON OF FORECAST ABILIT OF... ###
# Maximum Likelihood Estimation
# Optimal Score Estimation
# Reference Predictors

#We perform k-fold cross testing: 
#for each fold
#   estimate model on all other folds
#   predict on k'th fold and store error scores, resulting table with averaged margins


#What are possible reference predictors?
# - persistance predictor
# - uniform predictor
# - markov jump without covariates

library("MASS")
source("funcs_brier.R")
source("funcs_rps.R")
source("funcs_helping.R")
source("funcs_discrete_loglik.R")
source("funcs_diagonalization.R")

sourceCpp("forecast.cpp")

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

pars_MLE <- matrix(NA, ncol = length(beta00), nrow = k)
pars_MLE_cov <- matrix(NA, ncol = length(beta0), nrow = k)
pars_BRIER <- matrix(NA, ncol = length(beta00), nrow = k)
pars_BRIER_cov <- matrix(NA, ncol = length(beta0), nrow = k)
pars_RPS <- matrix(NA, ncol = length(beta00), nrow = k)
pars_RPS_cov <- matrix(NA, ncol = length(beta0), nrow = k)
pars_olr_cov <- matrix(NA, ncol = length(exo.cols)+m, nrow = k)
pars_olr <- matrix(NA, ncol = m, nrow = k)
pars_RPS_MLE_cov <- matrix(NA, ncol = length(beta0), nrow = k)
pars_BRIER_MLE_cov <- matrix(NA, ncol = length(beta0), nrow = k)

RPS_err = matrix(NA, ncol = k, nrow = 14)
log_err = matrix(NA, ncol = k, nrow = 14)
Brier_e = matrix(NA, ncol = k, nrow = 14)

#RUN K-FOLD PREDICTION
for(i in 1:(length(ints)-1)){
  start = ints[i]; 
  end = ints[i+1]; 
  
  #Estimation and test set splits
  pred.idx = rand.idx[start:end]
  d.test = d[pred.idx, ]
  d.train = d[-pred.idx, ]

  #Estimate on training data
  jump_MLE = mle.no.covs(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t, beta = beta00,  method = "BFGS")
  jump_MLE_cov = mle.cov(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t, z = as.matrix(d.train[,exo.cols]), beta = beta0, method = "BFGS")
  jump_BRIER = brier.estim(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t, beta = beta00, method = "BFGS")
  jump_BRIER_cov = brier.estim.cov(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t, z = as.matrix(d.train[,exo.cols]), beta = beta0, method = "BFGS")
  jump_RPS = rps.estim(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t, beta = beta00, method = "BFGS")
  jump_RPS_cov = rps.estim.cov(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t, z = as.matrix(d.train[,exo.cols]), beta = beta0, method = "BFGS")
  olr = polr(as.factor(s2) ~ as.factor(s1) + t, data = d.train, method = "logistic")
  olr_cov = polr(as.factor(s2) ~ as.factor(s1) + t + MBT.norm + speed.norm + profil.norm + steel.norm + invRad.norm, data = d.train, method = "logistic")
  jump_BRIER_MLE_cov = brier.loglik.estim.cov(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t, z = as.matrix(d.train[,exo.cols]), beta = beta0, method = "BFGS")
  jump_RPS_MLE_cov = rps.loglik.estim.cov(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t, z = as.matrix(d.train[,exo.cols]), beta = beta0, method = "BFGS")
  
  #Store estimates and loss value
  pars_MLE[i,] = jump_MLE$lambda 
  pars_MLE_cov[i,] = jump_MLE_cov$lambda
  pars_BRIER[i,] = jump_BRIER$lambda 
  pars_BRIER_cov[i,] = jump_BRIER_cov$lambda 
  pars_RPS[i,] = jump_RPS$lambda 
  pars_RPS_cov[i,] = jump_RPS_cov$lambda 
  pars_olr[i,] = olr$coefficients
  pars_olr_cov[i,] = olr_cov$coefficients
  pars_RPS_MLE_cov[i,] = jump_RPS_MLE_cov$lambda
  
  #Make predictions
  preds_MLE = jump_prediction(m, s1 = d.test$`s-`, u = d.test$t, jump_MLE$lambda)
  preds_MLE_cov = jump_prediction_cov(m, s1 = d.test$`s-`, u = d.test$t, z = as.matrix(d.test[,exo.cols]), jump_MLE_cov$lambda)
  preds_BRIER = jump_prediction(m, s1 = d.test$`s-`, u = d.test$t, jump_BRIER$lambda)
  preds_BRIER_cov = jump_prediction_cov(m, s1 = d.test$`s-`, u = d.test$t, z = as.matrix(d.test[,exo.cols]), jump_BRIER_cov$lambda)
  preds_RPS = jump_prediction(m, s1 = d.test$`s-`, u = d.test$t, jump_RPS$lambda)
  preds_RPS_cov = jump_prediction_cov(m, s1 = d.test$`s-`, u = d.test$t, z = as.matrix(d.test[,exo.cols]), jump_RPS_cov$lambda)
  preds_uni = uniform_prediction(m,  d.test$`s-`)
  preds_per = make_Ptu_obs(m,  d.test$`s-`)
  empi_foo = naive.empirical.pred(d.train$s2)
  pred_empi = matrix(rep(empi_foo, NROW(d.test)), ncol = m, byrow = T)
  pred_empi_corr = t(sapply(X = d.test$s1,  FUN = corr.empirical.pred, m = m, naive = empi_foo))
  preds_olr = predict(olr, newdata = d.test, type = "p")
  preds_olr_cov = predict(olr_cov, newdata = d.test, type = "p")
  preds_RPS_MLE_cov = jump_prediction_cov(m, s1 = d.test$`s-`, u = d.test$t, z = as.matrix(d.test[,exo.cols]), jump_RPS_MLE_cov$lambda)
  preds_BRIER_MLE_cov = jump_prediction_cov(m, s1 = d.test$`s-`, u = d.test$t, z = as.matrix(d.test[,exo.cols]), jump_BRIER_MLE_cov$lambda)
  
  #Evaluate predictions in terms of RPS and log score and Brier score
  obs = make_Ptu_obs(m, d.test$`s`)
  log_err[1, i] = mean(logscore_vectors(m, preds_MLE, obs))
  log_err[2, i] = mean(logscore_vectors(m, preds_MLE_cov, obs))
  log_err[3, i] = mean(logscore_vectors(m, preds_BRIER, obs))
  log_err[4, i] = mean(logscore_vectors(m, preds_BRIER_cov, obs))
  log_err[5, i] = mean(logscore_vectors(m, preds_RPS, obs))
  log_err[6, i] = mean(logscore_vectors(m, preds_RPS_cov, obs))
  log_err[7, i] = mean(logscore_vectors(m, preds_uni, obs))
  log_err[8, i] = mean(logscore_vectors(m, preds_per, obs))
  log_err[9, i] = mean(logscore_vectors(m, pred_empi, obs))
  log_err[10, i] = mean(logscore_vectors(m, pred_empi_corr, obs))
  log_err[11, i] = mean(logscore_vectors(m, preds_olr, obs))
  log_err[12, i] = mean(logscore_vectors(m, preds_olr_cov, obs))
  log_err[14, i] = mean(logscore_vectors(m, preds_RPS_MLE_cov, obs))
  log_err[13, i] = mean(logscore_vectors(m, preds_BRIER_MLE_cov, obs))
  
  Brier_e[1, i] = mean(BrierScore_vectors( preds_MLE, obs))
  Brier_e[2, i] = mean(BrierScore_vectors( preds_MLE_cov, obs))
  Brier_e[3, i] = mean(BrierScore_vectors( preds_BRIER, obs))
  Brier_e[4, i] = mean(BrierScore_vectors( preds_BRIER_cov, obs))
  Brier_e[5, i] = mean(BrierScore_vectors( preds_RPS, obs))
  Brier_e[6, i] = mean(BrierScore_vectors( preds_RPS_cov, obs))
  Brier_e[7, i] = mean(BrierScore_vectors( preds_uni, obs))
  Brier_e[8, i] = mean(BrierScore_vectors( preds_per, obs))
  Brier_e[9, i] = mean(BrierScore_vectors( pred_empi, obs))
  Brier_e[10, i] = mean(BrierScore_vectors( pred_empi_corr, obs))
  Brier_e[11, i] = mean(BrierScore_vectors( preds_olr, obs))
  Brier_e[12, i] = mean(BrierScore_vectors( preds_olr_cov, obs))
  Brier_e[14, i] = mean(BrierScore_vectors( preds_RPS_MLE_cov, obs))
  Brier_e[13, i] = mean(BrierScore_vectors( preds_BRIER_MLE_cov, obs))
  
  RPS_err[1, i] = mean(rps_vectors(m, preds_MLE, obs))
  RPS_err[2, i] = mean(rps_vectors(m, preds_MLE_cov, obs))
  RPS_err[3, i] = mean(rps_vectors(m, preds_BRIER, obs))
  RPS_err[4, i] = mean(rps_vectors(m, preds_BRIER_cov, obs))
  RPS_err[5, i] = mean(rps_vectors(m, preds_RPS, obs))
  RPS_err[6, i] = mean(rps_vectors(m, preds_RPS_cov, obs))
  RPS_err[7, i] = mean(rps_vectors(m, preds_uni, obs))
  RPS_err[8, i] = mean(rps_vectors(m, preds_per, obs))
  RPS_err[9, i] = mean(rps_vectors(m, pred_empi, obs))
  RPS_err[10, i] = mean(rps_vectors(m, pred_empi_corr, obs))
  RPS_err[11, i] = mean(rps_vectors(m, preds_olr, obs))
  RPS_err[12, i] = mean(rps_vectors(m, preds_olr_cov, obs))
  RPS_err[14, i] = mean(rps_vectors(m, preds_RPS_MLE_cov, obs))
  RPS_err[13, i] = mean(rps_vectors(m, preds_BRIER_MLE_cov, obs))
  
}

names = c("jump_MLE", "jump_MLE_cov", "jump_Brier", "jump_Brier_cov", "jump_RPS","jump_RPS_cov",
          "Uniform", "Persistence", "Empirical", "Empirical_corr", 
          "olr", "olr_cov", "jump_Brier_MLE_cov", "jump_RPS_MLE_cov")
errs = cbind(rowMeans(log_err),rowMeans(Brier_e),rowMeans(RPS_err))
rownames(errs) = names
colnames(errs) = c("Log S", "Brier", "RPS")
errs



