

#Purpose: Test script for organized estimation with RPS, Brier, and Likelihood.
# - The user should be able to choose which score, she will use for estimation. 
# - The user should be able to choose how she want to compute 
#the transient distribution of the finite space continuous time Markov chain, 
#either by Padé, Uniformization, or eigenvalue decomposition.
# - The user should be able to choose whether or not to include covariates.


rm(list = ls()) #clear memory

### Read data and define state space and covariates ###
d = read.csv("defect_data.csv")
states <- c(1,2,3,4,5) #define states
m <- length(states) #number of states
track <- unique(d$Track) #define investigated tracks
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

### Read libaries ###
library(MASS)
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")

### EMPIRICAL DISTRIBUTION PREDICTORS ###
naive.empirical.pred = function(x){return(as.vector(table(x)/length(x)))}
corr.empirical.pred = function(m, idx, naive){
  if(idx > 1){
    prob = sum(naive[1:(idx-1)])
    naive[1:(idx-1)] = 0
    naive[idx:m] = naive[idx:m] + prob/(m-idx+1) 
  } 
  return(naive)
}

### INITIALIZE k-fold ###
set.seed(123)
k <- 10
ints = seq(1, NROW(d), ceiling(NROW(d)/k)-1)
ints[length(ints)] = NROW(d)
rand.idx <- sample(x = 1:NROW(d),size = NROW(d),replace = F)
PREDS = list()
PARS = list()
log_bin = F; rps_bin = F; brier_bin = F;
parameterizations = c("gerlang", "gerlang_relax", "free_upper_tri") 
scores = c("log", "rps", "brier", "all")
baselines = c("uniform", "persistence", "empirical_dist", "empirical_dist_corr", "olr", "olr_cov", "obs")
n_predictors = length(parameterizations)*length(scores)*2 + length(baselines) #number of predictors = number of mjp models + reference predictors
log_err = matrix(NA, ncol = k, nrow = n_predictors)
RPS_err = matrix(NA, ncol = k, nrow = n_predictors)
Brier_e = matrix(NA, ncol = k, nrow = n_predictors)

#RUN K-FOLD PREDICTION
for(i in 1:k){
  start_time <- Sys.time()
  start = ints[i]; 
  end = ints[i+1]; 
  
  #Estimation and test set splits
  pred.idx = rand.idx[start:end]
  d.test = d[pred.idx, ]
  d.train = d[-pred.idx, ]
  obs = make_Ptu_obs(m, d.test$s2)
  
  ####################################################
  ### Estimate and predict with all MJP variations ###
  ####################################################
  count = 0
  for(gen in parameterizations){
    if(gen == "gerlang" | gen == "gerlang_relax"){beta_base = rep(0.25, (m-1))}
    if(gen == "free_upper_tri"){beta_base = rep(0.25, m*(m-1)/2)}
    
    for(cov in c(F, T)){
      if(cov ){ beta = c(beta_base, rep(0.1,length(exo.cols))); } 
      if(!cov){ beta = beta_base; }
      
      for(score in scores){
        count = count + 1
        if(score == "log" | score == "all"){log_bin = T}
        if(score == "rps" | score == "all"){rps_bin = T}
        if(score == "brier" | score == "all"){brier_bin = T}
        
        #Estimate and store model pars
        foo = optim( par = beta, fn = MJP_score, m = m, s1 = d.train$s1, s2 = d.train$s2, u = d.train$t, z = as.matrix(d.train[,exo.cols]), generator = gen, covs_bin = cov, likelihood_bin = log_bin, rps_bin = rps_bin, brier_bin = brier_bin, transient_dist_method = "eigenvalue_decomp", method = "BFGS", control = list(maxit = 1000)) #estimate model 
        nam = paste(gen,cov,score,sep = "_") #name of specific estimation
        if(i == 1){ PARS[[nam]] = foo$par} else { PARS[[nam]] = rbind(PARS[[nam]],foo$par)} #store model estimates
        
        #Make and store predictions and compute error scores
        pred = MJP_predict(m = m, s1 = d.test$s1, u = d.test$t, pars = foo$par, z = as.matrix(d.test[,exo.cols]), generator = gen, covs_bin = cov, transient_dist_method = "eigenvalue_decomp")
        if(i == 1){ PREDS[[nam]] = pred} else { PREDS[[nam]] = rbind(PREDS[[nam]],pred)} #store predictions
        log_err[count, i] = mean(logscore_vectors(m, pred, obs))
        RPS_err[count, i] = mean(rps_vectors(m, pred, obs))
        Brier_e[count, i] = mean(BrierScore_vectors(m, pred, obs))
        
        log_bin = F; rps_bin = F; brier_bin = F; #reset score bools
      }
    }
  }
  
  ############################
  ### REFERENCE PREDICTORS ###
  ############################
  #Estimate multinomial ordinal regression
  olr = MASS::polr(as.factor(s2) ~ as.factor(s1) + t, data = d.train, method = "logistic")
  olr_cov = MASS::polr(as.factor(s2) ~ as.factor(s1) + t + MBT.norm + speed.norm + profil.norm + steel.norm + invRad.norm, data = d.train, method = "logistic")
  for(j in baselines){
    count = count + 1;
    #Make predictions
    if(j == "uniform"){ pred = uniform_prediction(m,  d.test$s1); } 
    if(j == "persistence"){ pred = make_Ptu_obs(m,  d.test$s1); }
    if(j == "empirical_dist"){ 
      empi_foo = naive.empirical.pred(d.train$s2); 
      pred = matrix(rep(empi_foo, NROW(d.test)), ncol = m, byrow = T);
    }
    if(j == "empirical_dist_corr"){ 
      empi_foo = naive.empirical.pred(d.train$s2); 
      pred = t(sapply(X = d.test$s1,  FUN = corr.empirical.pred, m = m, naive = empi_foo));
    }
    if(j == "olr"){ 
      if(i == 1){ PARS[[j]] = olr$coefficients} else { PARS[[j]] = rbind(PARS[[j]],olr$coefficients)} #store model estimates
      pred = predict(olr, newdata = d.test, type = "p");
      }
    if(j == "olr_cov"){ 
      if(i == 1){ PARS[[j]] = olr_cov$coefficients} else { PARS[[j]] = rbind(PARS[[j]],olr_cov$coefficients)} #store model estimates
      pred = predict(olr_cov, newdata = d.test, type = "p"); }
    if(j == "obs"){ pred = obs; }
    
    if(i == 1){ PREDS[[j]] = pred} else { PREDS[[j]] = rbind(PREDS[[j]],pred)} 
    log_err[count, i] = mean(logscore_vectors(m, pred, obs))
    RPS_err[count, i] = mean(rps_vectors(m, pred, obs))
    Brier_e[count, i] = mean(BrierScore_vectors(m, pred, obs))
  }
  end_time <- Sys.time(); runtime <- end_time - start_time
  cat("Fold ", i, "runtime:", runtime, "\n") #approx 3.5 mins per fold
}

rownames(log_err) = c(names(PARS)[1:(length(parameterizations)*length(scores)*2)], baselines)
rownames(RPS_err) = c(names(PARS)[1:(length(parameterizations)*length(scores)*2)], baselines)
rownames(Brier_e) = c(names(PARS)[1:(length(parameterizations)*length(scores)*2)], baselines)

errs = cbind(rowMeans(log_err),rowMeans(Brier_e),rowMeans(RPS_err))
errs


# write.csv(PARS, file  = "results/estimated_model_pars.csv")
# write.csv(PREDS, file  = "results/model_preditions.csv")
write.csv(log_err, file  = "results/log_score.csv", row.names = T)
write.csv(Brier_e, file  = "results/brier_score.csv", row.names = T)
write.csv(RPS_err, file  = "results/rps_score.csv", row.names = T)
saveRDS(PARS, file = "results/estimated_model_pars.Rdata")
saveRDS(PREDS, file = "results/model_predictions.Rdata")


#comments: 
# rps and brier estimation on gerlang relax is meaningful, it just takes more iterations to reach convergence
# sensitivity on intialization seems to be resolved







####################################################
####################################################
# #Example code of for using functions from cpp file
# z = as.matrix(d[, exo.cols])
# beta0 <- c(c(0.1,0.2,0.3,0.4),
#            #rep(0.25, (m-1)),
#            #c(0.1,0.2,0.3,0.4, 0.2,0.2,0.3,0.4,0.1,9),
#            rep(0.1,length(exo.cols)))  
# 
# MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", 
#                covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, 
#           transient_dist_method = "uniformization")
# x = 
#   optim( par = beta0, fn = MJP_score, m = m, s1 = d$s1, s2 = d$s2, u = d$t, z = z, 
#          generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, 
#          transient_dist_method = "eigenvalue_decomp", 
#          method = "BFGS", control = list(maxit = 1000)) 
# pred = 
#   MJP_predict(m = m, s1 = d$s1, u = d$t, pars = x$par, z = z, generator = "gerlang",
#               covs_bin = T, transient_dist_method = "eigenvalue_decomp")
# obs = make_Ptu_obs(m, d$s2)
# mean(logscore_vectors(m, pred, obs))

