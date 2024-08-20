rm(list = ls()) #clear memory

#setwd("C:/Users/ATBD/OneDrive - Danmarks Tekniske Universitet/MarkovJump")
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
exo.cols2 <- c("speed.norm","profil.norm", "steel.norm", "invRad.norm")

idx.remove <- 
  (D$`s-` == 4 & D$`s` == 4 & D$u > 300) | 
  (D$`s-` == 4 & D$`s` == 5 & D$u > 300) | 
  (D$`s-` == 5 & D$`s` == 5 & D$u > 300) 


d <- D[!idx.remove, c("pos", "s-", "s", "u", exo.cols, "Track0")] #define data set of interest

d$s1 = d$`s-`
d$s2 = d$`s`
d$t = d$u/365
d$t_mgt = d$u/365*d$MBT.norm
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
beta02 <- c(rep(0.25, (m-1)), rep(0.1,length(exo.cols2))) 



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
pars_rps_skill_cov <- matrix(NA, ncol = length(beta0), nrow = k)
pars_MLE_cov_A2 <- matrix(NA, ncol = length(beta0), nrow = k)
pars_MLE_BRIER_RPS_cov_A2 <- matrix(NA, ncol = length(beta0), nrow = k)

pars_MLE_cov_timeMGT <- matrix(NA, ncol = length(beta02), nrow = k)

n_predictors = 17

RPS_err = matrix(NA, ncol = k, nrow = n_predictors)
log_err = matrix(NA, ncol = k, nrow = n_predictors)
Brier_e = matrix(NA, ncol = k, nrow = n_predictors)

pred_list = vector("list",n_predictors + 1)
for (i in 1:n_predictors) {
  pred_list[[i]] <- matrix(, nrow=0, ncol=5)
}


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
  jump_MLE_cov_A2 = mle.cov.A2(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t, z = as.matrix(d.train[,exo.cols]), beta = beta0, method = "BFGS")
  jump_MLE_BRIER_RPS_cov_A2 = loglik.brier.rps.estim.cov.A2(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t, z = as.matrix(d.train[,exo.cols]), beta = beta0, method = "BFGS")
  
  jump_MLE_cov_timeMGT = mle.cov(m = m, s1 = d.train$`s-`, s2 = d.train$s, u = d.train$t_mgt, z = as.matrix(d.train[,exo.cols2]), beta = beta02, method = "BFGS")
  
  
  
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
  pars_MLE_cov_A2[i,] = jump_MLE_cov_A2$lambda
  pars_MLE_BRIER_RPS_cov_A2[i,] = jump_MLE_BRIER_RPS_cov_A2$lambda
  
  pars_MLE_cov_timeMGT[i,] = jump_MLE_cov_timeMGT$lambda
  
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
  preds_MLE_cov_A2 = jump_prediction_cov_A2(m, s1 = d.test$`s-`, u = d.test$t, z = as.matrix(d.test[,exo.cols]), jump_MLE_cov_A2$lambda)
  preds_MLE_BRIER_RPS_cov_A2 = jump_prediction_cov_A2(m, s1 = d.test$`s-`, u = d.test$t, z = as.matrix(d.test[,exo.cols]), jump_MLE_BRIER_RPS_cov_A2$lambda)
  
  preds_MLE_cov_timeMGT = jump_prediction_cov(m, s1 = d.test$`s-`, u = d.test$t_mgt, z = as.matrix(d.test[,exo.cols2]), jump_MLE_cov_timeMGT$lambda)
  
  
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
  log_err[13, i] = mean(logscore_vectors(m, preds_BRIER_MLE_cov, obs))
  log_err[14, i] = mean(logscore_vectors(m, preds_RPS_MLE_cov, obs))
  log_err[15, i] = mean(logscore_vectors(m, preds_MLE_cov_A2, obs))
  log_err[16, i] = mean(logscore_vectors(m, preds_MLE_BRIER_RPS_cov_A2, obs))
  
  log_err[17, i] = mean(logscore_vectors(m, preds_MLE_cov_timeMGT, obs))
  
  
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
  Brier_e[13, i] = mean(BrierScore_vectors( preds_BRIER_MLE_cov, obs))
  Brier_e[14, i] = mean(BrierScore_vectors( preds_RPS_MLE_cov, obs))
  Brier_e[15, i] = mean(BrierScore_vectors( preds_MLE_cov_A2, obs))
  Brier_e[16, i] = mean(BrierScore_vectors( preds_MLE_BRIER_RPS_cov_A2, obs))
  
  Brier_e[17, i] = mean(BrierScore_vectors( preds_MLE_cov_timeMGT, obs))
  
  
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
  RPS_err[13, i] = mean(rps_vectors(m, preds_BRIER_MLE_cov, obs))
  RPS_err[14, i] = mean(rps_vectors(m, preds_RPS_MLE_cov, obs))
  RPS_err[15, i] = mean(rps_vectors(m, preds_MLE_cov_A2, obs))
  RPS_err[16, i] = mean(rps_vectors(m, preds_MLE_BRIER_RPS_cov_A2, obs))
  
  RPS_err[17, i] = mean(rps_vectors(m, preds_MLE_cov_timeMGT, obs))
  
  pred_list[[1]] = rbind(pred_list[[1]], preds_MLE)
  pred_list[[2]] = rbind(pred_list[[2]], preds_MLE_cov)
  pred_list[[3]] = rbind(pred_list[[3]], preds_BRIER)
  pred_list[[4]] = rbind(pred_list[[4]], preds_BRIER_cov)
  pred_list[[5]] = rbind(pred_list[[5]], preds_RPS)
  pred_list[[6]] = rbind(pred_list[[6]], preds_RPS_cov)
  pred_list[[7]] = rbind(pred_list[[7]], preds_uni)
  pred_list[[8]] = rbind(pred_list[[8]], preds_per)
  pred_list[[9]] = rbind(pred_list[[9]], pred_empi)
  pred_list[[10]] = rbind(pred_list[[10]], pred_empi_corr)
  pred_list[[11]] = rbind(pred_list[[11]], preds_olr)
  pred_list[[12]] = rbind(pred_list[[12]], preds_olr_cov)
  pred_list[[13]] = rbind(pred_list[[13]], preds_BRIER_MLE_cov)
  pred_list[[14]] = rbind(pred_list[[14]], preds_RPS_MLE_cov)
  pred_list[[15]] = rbind(pred_list[[15]], preds_MLE_cov_A2)
  pred_list[[16]] = rbind(pred_list[[16]], preds_MLE_BRIER_RPS_cov_A2)
  
  pred_list[[17]] = rbind(pred_list[[17]], preds_MLE_cov_timeMGT)
  
  
  pred_list[[n_predictors+1]] = rbind(pred_list[[n_predictors+1]], obs)
  
}

names = c("jump_MLE", "jump_MLE_cov", "jump_Brier", "jump_Brier_cov", "jump_RPS","jump_RPS_cov",
          "Uniform", "Persistence", "Empirical", "Empirical_corr", 
          "olr", "olr_cov", "jump_Brier_MLE_cov", "jump_RPS_MLE_cov", 
          "jump_MLE_cov_A2", "jump_MLE_BRIER_RPS_cov_A2", "timeMGT")
errs = cbind(rowMeans(log_err),rowMeans(Brier_e),rowMeans(RPS_err))
rownames(errs) = names
colnames(errs) = c("Log S", "Brier", "RPS")
errs

## ESTIMATION COMMENT: There seems to be some irregularity on the Brier score and RPS when using
#A2 parameterization of the generator. I suspect it has something to do with the optimizer. Maybe
#logarithm can be taken in each iteration to stabilize, but null-solutions need to be considered
#in that case.

#How good are the respective models to forecast 0 or 1 events?
bin_predictions = function(col1, col2, preds){ #convert multicategorial predictions into binary for event in col1 to col2
  if(col1 == col2){
    res = preds[,col1]
  } else if(col1 < col2){
    res = rowSums(preds[,col1:col2])
  }
  return(res)
}

zero_bin_obs = function(col1, col2, obs){
  idx = apply(obs, 1, function(row) which(row == 1))
  return(idx >= col1 & idx <= col2)
}




####################
### THE TRIPTYCH ###
####################
#severe stuff
source("reliability_diag.R")
rd = ReliabilityDiagram2(bin_predictions(m-1,m,pred_list[[2]])-0.000001, zero_bin_obs(m-1,m,pred_list[[length(pred_list)]]), plot = T, plot.refin = F, attributes = T)
rd = ReliabilityDiagram2(bin_predictions(m-1,m,pred_list[[4]])-0.000001, zero_bin_obs(m-1,m,pred_list[[length(pred_list)]]), plot = T, plot.refin = F, attributes = T)
rd = ReliabilityDiagram2(bin_predictions(m-1,m,pred_list[[6]])-0.000001, zero_bin_obs(m-1,m,pred_list[[length(pred_list)]]), plot = T, plot.refin = F, attributes = T)
rd = ReliabilityDiagram2(bin_predictions(m-1,m,pred_list[[12]])-0.000001, zero_bin_obs(m-1,m,pred_list[[length(pred_list)]]), plot = T, plot.refin = F, attributes = T)
rd = ReliabilityDiagram2(bin_predictions(m-1,m,pred_list[[14]])-0.000001, zero_bin_obs(m-1,m,pred_list[[length(pred_list)]]), plot = T, plot.refin = F, attributes = T)
rd = ReliabilityDiagram2(bin_predictions(m-1,m,pred_list[[15]])-0.000001, zero_bin_obs(m-1,m,pred_list[[length(pred_list)]]), plot = T, plot.refin = F, attributes = T)


library("ROSE")
roc.curve(zero_bin_obs(m-1,m-1,pred_list[[length(pred_list)]]), bin_predictions(m-1,m-1,pred_list[[2]]))
roc.curve(zero_bin_obs(m-1,m-1,pred_list[[length(pred_list)]]), bin_predictions(m-1,m-1,pred_list[[12]]), add.roc = T)
roc.curve(zero_bin_obs(m-1,m,pred_list[[length(pred_list)]]), bin_predictions(m-1,m,pred_list[[2]]))
roc.curve(zero_bin_obs(m-1,m,pred_list[[length(pred_list)]]), bin_predictions(m-1,m,pred_list[[12]]), add.roc = T)
roc.curve(zero_bin_obs(m-1,m,pred_list[[length(pred_list)]]), bin_predictions(m-1,m,pred_list[[15]]), add.roc = T)


library(murphydiagram)
murphydiagram(bin_predictions(m-1,m,pred_list[[15]]), bin_predictions(m-1,m,pred_list[[12]]), bin_predictions(m-1,m,pred_list[[length(pred_list)]]))



#mild stuff
source("reliability_diag.R")
col1 = 1; col2 = 2;
rd = ReliabilityDiagram2(bin_predictions(col1,col2,pred_list[[2]]), zero_bin_obs(col1,col2,pred_list[[15]]), plot = T, plot.refin = F, attributes = T)
rd = ReliabilityDiagram2(bin_predictions(col1,col2,pred_list[[12]]), zero_bin_obs(col1,col2,pred_list[[15]]), plot = T, plot.refin = F, attributes = T)


library("ROSE")
roc.curve(zero_bin_obs(col1,col2,pred_list[[15]]), bin_predictions(col1,col2,pred_list[[2]]))
roc.curve(zero_bin_obs(col1,col2,pred_list[[15]]), bin_predictions(col1,col2,pred_list[[12]]), add.roc = T)


library(murphydiagram)
murphydiagram(bin_predictions(col1,col2,pred_list[[2]]), bin_predictions(col1,col2,pred_list[[12]]), bin_predictions(col1,col2,pred_list[[15]]))











