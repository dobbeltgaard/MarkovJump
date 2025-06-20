

#Purpose: Test script for organized estimation with RPS, Brier, and Likelihood.
# - The user should be able to choose which score, she will use for estimation. 
# - The user should be able to choose how she want to compute 
#the transient distribution of the finite space continuous time Markov chain, 
#either by Padé, Uniformization, or eigenvalue decomposition.
# - The user should be able to choose whether or not to include covariates.


rm(list = ls()) #clear memory

d = read.csv("defect_data.csv")
states <- c(1,2,3,4,5)
m <- length(states) 
track <- unique(d$Track)
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

library(MASS)
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")

#library(TMB)
#compile("FUNCS_MJP_with_TMB_warped.cpp")
#dyn.load(dynlib("FUNCS_MJP_with_TMB_warped"))

library(TMB)
compile("FUNCS_MJP_with_TMB.cpp")
compile("FUNCS_MJP_with_TMB_warped.cpp")



switch_model_dll <- function(to_load) {
  other <- if (to_load == "FUNCS_MJP_with_TMB") "FUNCS_MJP_with_TMB_warped" else "FUNCS_MJP_with_TMB"
  if (other %in% names(getLoadedDLLs())) dyn.unload(dynlib(other))
  if (!(to_load %in% names(getLoadedDLLs()))) dyn.load(dynlib(to_load))
}


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
smart.empirical.pred = function(s1, s2){
  trans_counts <- table(s1, s2)
  trans_matrix <- prop.table(trans_counts, margin = 1)  # normalize by row
  return(as.matrix(trans_matrix))
}


### INITIALIZE k-fold ###
set.seed(123)
k <- 5
ints = seq(1, NROW(d), ceiling(NROW(d)/k)-1)
ints[length(ints)] = NROW(d)
rand.idx <- sample(x = 1:NROW(d),size = NROW(d),replace = F)
PREDS = list()
PARS = list()
log_bin = F; rps_bin = F; brier_bin = F;
parameterizations = c("gerlang", "gerlang_relax", "free_upper_tri") 
scores = c("log", "rps", "all")
links = c("softplus")#c("exp", "softplus", "square")
baselines = c("uniform", "persistence", "empirical_dist", "empirical_dist_corr","empirical_dist_smart" , "olr", "olr_cov","opr", "opr_cov","ocllr","ocllr_cov", "obs")
n_predictors = length(parameterizations)*length(scores)*2*length(links)^2 + length(baselines) #number of predictors = number of mjp models + reference predictors
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
  for(warping in c("no_warp","warp")){
    if(warping == "no_warp"){ 
      warp = F
      xi = c()
      switch_model_dll("FUNCS_MJP_with_TMB")
      }
    if(warping == "warp"){ 
      warp = T
      xi = runif(m-1,0,1)
      switch_model_dll("FUNCS_MJP_with_TMB_warped")
      }
    for(gen in parameterizations){
      if(gen == "gerlang"){generator_type = 0; beta_base = runif(m-1,0,1); }
      if(gen == "gerlang_relax"){generator_type = 1; beta_base = runif(m-1,0,1); }
      if(gen == "free_upper_tri"){generator_type = 2; beta_base = runif(m*(m-1)/2,0,1); }
      for(cov in c(F, T)){
        if(cov ){ beta = c(beta_base, xi, runif(length(exo.cols), 0,1)); } 
        if(!cov){ beta = c(beta_base,xi); }
        for(baselink in links){
          for(covslink in links){
              for(score in scores){
                if(score == "log" | score == "all"){log_bin = T}
                if(score == "rps" | score == "all"){rps_bin = T}
                #if(score == "brier" | score == "all"){brier_bin = T}
                  count = count + 1
                  
                  #TMB implementation can only be used with softplus link for now!!
                  data.train <- list(s1 = d.train$s1,s2 = d.train$s2,u = d.train$t,z = as.matrix(d.train[, exo.cols]),m = m,generator_type = generator_type,cov_type = as.integer(cov),use_log_score = as.integer(log_bin),use_rps_score = as.integer(rps_bin), use_brier_score = 0)
                  parameters <- list(theta = beta)
                  if(warping == "no_warp"){ l <- MakeADFun(data = data.train, parameters = parameters, DLL = "FUNCS_MJP_with_TMB")}
                  if(warping == "warp"){ l <- MakeADFun(data = data.train, parameters = parameters, DLL = "FUNCS_MJP_with_TMB_warped")}
                  
                  #Estimate and store model pars
                  nam = paste(gen,cov,baselink,covslink,score,warping,sep = "_") #name of specific estimation
                  print(nam)
                  #foo = optim( par = beta, fn = MJP_score, m = m, s1 = d.train$s1, s2 = d.train$s2, u = d.train$t, z = as.matrix(d.train[,exo.cols]), generator = gen, link_type_base = baselink, link_type_covs = covslink, covs_bin = cov, likelihood_bin = log_bin, rps_bin = rps_bin, brier_bin = brier_bin, transient_dist_method = "pade", method = "BFGS", control = list(maxit = 1000)) #estimate model 
                  foo <- nlminb(l$par, l$fn, l$gr, l$he)
                  if(i == 1){ PARS[[nam]] = foo$par} else { PARS[[nam]] = rbind(PARS[[nam]],foo$par)} #store model estimates
                  
                  #Make and store predictions and compute error scores
                  pred = MJP_predict(m = m, s1 = d.test$s1, u = d.test$t, pars = foo$par, z = as.matrix(d.test[,exo.cols]), generator = gen, link_type_base = baselink, link_type_covs = covslink, covs_bin = cov, transient_dist_method = "pade", warping = warp)
                  if(i == 1){ PREDS[[nam]] = pred} else { PREDS[[nam]] = rbind(PREDS[[nam]],pred)} #store predictions
                  log_err[count, i] = mean(logscore_vectors(m, pred, obs))
                  RPS_err[count, i] = mean(rps_vectors(m, pred, obs))
                  Brier_e[count, i] = mean(BrierScore_vectors(m, pred, obs))
                  
                  log_bin = F; rps_bin = F; brier_bin = F; #reset score bools
                }
            }
          }
        }
      }
  }
  
  ############################
  ### REFERENCE PREDICTORS ###
  ############################
  #Estimate multinomial ordinal regression
  olr = MASS::polr(as.factor(s2) ~ as.factor(s1) + t, data = d.train, method = "logistic")
  olr_cov = MASS::polr(as.factor(s2) ~ as.factor(s1) + t + MBT.norm + speed.norm + profil.norm + steel.norm + invRad.norm, data = d.train, method = "logistic")
  opr = MASS::polr(as.factor(s2) ~ as.factor(s1) + t, data = d.train, method = "probit")
  opr_cov = MASS::polr(as.factor(s2) ~ as.factor(s1) + t + MBT.norm + speed.norm + profil.norm + steel.norm + invRad.norm, data = d.train, method = "probit")
  ocllr = MASS::polr(as.factor(s2) ~ as.factor(s1) + t, data = d.train, method = "cloglog")
  ocllr_cov = MASS::polr(as.factor(s2) ~ as.factor(s1) + t + MBT.norm + speed.norm + profil.norm + steel.norm + invRad.norm, data = d.train, method = "cloglog")
  
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
    if(j == "empirical_dist_smart"){
      empi_foo = smart.empirical.pred(d.train$s1,d.train$s2)
      pred = empi_foo[d.test$s1, ]
    }
    if(j == "olr"){ 
      if(i == 1){ PARS[[j]] = olr$coefficients} else { PARS[[j]] = rbind(PARS[[j]],olr$coefficients)} #store model estimates
      pred = predict(olr, newdata = d.test, type = "p");
      }
    if(j == "olr_cov"){ 
      if(i == 1){ PARS[[j]] = olr_cov$coefficients} else { PARS[[j]] = rbind(PARS[[j]],olr_cov$coefficients)} #store model estimates
      pred = predict(olr_cov, newdata = d.test, type = "p"); }
    if(j == "opr"){ 
      if(i == 1){ PARS[[j]] = opr$coefficients} else { PARS[[j]] = rbind(PARS[[j]],opr$coefficients)} #store model estimates
      pred = predict(opr, newdata = d.test, type = "p"); }
    if(j == "opr_cov"){ 
      if(i == 1){ PARS[[j]] = opr_cov$coefficients} else { PARS[[j]] = rbind(PARS[[j]],opr_cov$coefficients)} #store model estimates
      pred = predict(opr_cov, newdata = d.test, type = "p"); }
    if(j == "ocllr"){ 
      if(i == 1){ PARS[[j]] = ocllr$coefficients} else { PARS[[j]] = rbind(PARS[[j]],ocllr$coefficients)} #store model estimates
      pred = predict(ocllr, newdata = d.test, type = "p"); }
    if(j == "ocllr_cov"){ 
      if(i == 1){ PARS[[j]] = ocllr_cov$coefficients} else { PARS[[j]] = rbind(PARS[[j]],ocllr_cov$coefficients)} #store model estimates
      pred = predict(ocllr_cov, newdata = d.test, type = "p"); }
    if(j == "obs"){ pred = obs; }
    if(i == 1){ PREDS[[j]] = pred} else { PREDS[[j]] = rbind(PREDS[[j]],pred)} 
    log_err[count, i] = mean(logscore_vectors(m, pred, obs))
    RPS_err[count, i] = mean(rps_vectors(m, pred, obs))
    Brier_e[count, i] = mean(BrierScore_vectors(m, pred, obs))
  }
  end_time <- Sys.time(); runtime <- end_time - start_time
  cat("Fold ", i, "runtime:", runtime, "\n") #approx 3.5 mins per fold
}

rownames(log_err) = c(names(PARS)[1:(length(parameterizations)*length(scores)*2*length(links)^2)], baselines)
rownames(RPS_err) = c(names(PARS)[1:(length(parameterizations)*length(scores)*2*length(links)^2)], baselines)
rownames(Brier_e) = c(names(PARS)[1:(length(parameterizations)*length(scores)*2*length(links)^2)], baselines)

errs = cbind(rowMeans(log_err),rowMeans(Brier_e),rowMeans(RPS_err))
errs


# write.csv(PARS, file  = "results/estimated_model_pars.csv")
# write.csv(PREDS, file  = "results/model_preditions.csv")
write.csv(log_err, file  = "results/log_score_V3.csv", row.names = T)
write.csv(Brier_e, file  = "results/brier_score_V3.csv", row.names = T)
write.csv(RPS_err, file  = "results/rps_score_V3.csv", row.names = T)
saveRDS(PARS, file = "results/estimated_model_pars_V3.Rdata")
saveRDS(PREDS, file = "results/model_predictions_V3.Rdata")


#comments: 
# rps and brier estimation on gerlang relax is meaningful, it just takes more iterations to reach convergence
# sensitivity on intialization seems to be resolved



rm(list = ls()) #clear memory
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states); track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
smart.empirical.pred = function(s1, s2){
  trans_counts <- table(s1, s2)
  trans_matrix <- prop.table(trans_counts, margin = 1)  # normalize by row
  return(as.matrix(trans_matrix))
}

library(MASS)
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")
obs = make_Ptu_obs(m, d$s2)
library(TMB)
compile("FUNCS_MJP_with_TMB.cpp")
dyn.load(dynlib("FUNCS_MJP_with_TMB"))

#NAIVE
empi_foo = smart.empirical.pred(d$s1,d$s2)
pred2 = empi_foo[d$s1, ]
err_naive = rps_vectors(m, pred2, obs)
mean(err_naive) #0.4268079
mean(logscore_vectors(m, pred2, obs)) #0.9689272

#MJP (A3+Covariates)
data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = as.integer(T))
beta = runif(m*(m-1)/2 + length(exo.cols),0,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = "FUNCS_MJP_with_TMB")
foo <- nlminb(l$par, l$fn, l$gr, l$he)
pred1 = MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "free_upper_tri", link_type_base = "softplus", link_type_covs = "softplus", covs_bin = T, transient_dist_method = "pade")
err_mjp = rps_vectors(m, pred1, obs)
mean(err_mjp) #0.4378375 (0.4364708 when purely rps optimized)

#MJP (A3+Covariates+state scaling of time)
u_scale = d$t*exp(1*d$s1)
data <- list(s1 = d$s1,s2 = d$s2,u = u_scale,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = as.integer(T))
beta = runif(m*(m-1)/2 + length(exo.cols),0,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = "FUNCS_MJP_with_TMB")
foo <- nlminb(l$par, l$fn, l$gr, l$he)
pred3 = MJP_predict(m = m, s1 = d$s1, u = u_scale, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "free_upper_tri", link_type_base = "softplus", link_type_covs = "softplus", covs_bin = T, transient_dist_method = "pade")
err_mjp3 = rps_vectors(m, pred3, obs)
mean(err_mjp3) #0.43352

#MJP (A3+Covariates+warped)
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
compile("FUNCS_MJP_with_TMB_warped.cpp")
dyn.load(dynlib("FUNCS_MJP_with_TMB_warped"))
data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = as.integer(T), use_log_score = 1,use_rps_score = 0, use_brier_score = 0)
beta = runif(m*(m-1)/2 + m-1 + length(exo.cols),0,2); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = "FUNCS_MJP_with_TMB_warped", hessian=T)
foo <- nlminb(l$par, l$fn, l$gr, l$he)
pred4 = MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "free_upper_tri", link_type_base = "softplus", link_type_covs = "softplus", covs_bin = T, transient_dist_method = "pade", warping=T)
err_mjp4 = rps_vectors(m, pred4, obs)
mean(err_mjp4) #0.4166109 (when purely RPS optimized)
mean(logscore_vectors(m, (pred4+pred2)/2, obs)) #1.017219 (when rps+log optimized) ( 1.619893 when only rps optimized)


#MJP (A2+Covariates+warped)
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
compile("FUNCS_MJP_with_TMB_warped.cpp")
dyn.load(dynlib("FUNCS_MJP_with_TMB_warped"))
data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(1),cov_type = as.integer(T), use_log_score = 1,use_rps_score = 1, use_brier_score = 0)
beta = runif(m-1 + m-1 + length(exo.cols),0,2); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = "FUNCS_MJP_with_TMB_warped", hessian=T)
foo <- nlminb(l$par, l$fn, l$gr, l$he)
pred4 = MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "gerlang_relax", link_type_base = "softplus", link_type_covs = "softplus", covs_bin = T, transient_dist_method = "pade", warping=T)
err_mjp4 = rps_vectors(m, pred4, obs)
mean(err_mjp4) #0.4166109 (when purely RPS optimized)
mean(logscore_vectors(m, (pred4+pred2)/2, obs)) #1.017219 (when rps+log optimized) ( 1.619893 when only rps optimized)




#test covariates (if all 5 are ok) yes thats better with linespeed as well
#test if warping is necessary. Yes that is better with warping. However, warping is more sensitive to score choice than if no warping
#test scoring (if rps + log S is ok). No, slightly worse than naive.
#test conditional expectation of sojourn times


#IDEA: Use cumulative link models as benchmark.
#Naive predictor is just included to make an ensemble forecast, both for MJP and cumu. links.
#Consider making a third baseline (for instance non-parametric).
#Optimize wrt single metrics rather than their sum


#MJP (A3+warped)
compile("FUNCS_MJP_with_TMB_warped.cpp")
dyn.load(dynlib("FUNCS_MJP_with_TMB_warped"))
data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = as.integer(F))
beta = runif(m*(m-1)/2 + m-1,0,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = "FUNCS_MJP_with_TMB_warped", hessian=F)
foo <- nlminb(l$par, l$fn, l$gr, l$he)
pred6 = MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "free_upper_tri", link_type_base = "softplus", link_type_covs = "softplus", covs_bin = F, transient_dist_method = "pade", warping=T)
err_mjp6 = rps_vectors(m, pred6, obs)
mean(err_mjp6) #0.4411932


#MJP (A1+Covariates+warped)
compile("FUNCS_MJP_with_TMB_warped.cpp")
dyn.load(dynlib("FUNCS_MJP_with_TMB_warped"))
data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(0),cov_type = as.integer(T))
beta = runif(m-1 + m-1 + length(exo.cols),0,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = "FUNCS_MJP_with_TMB_warped", hessian = T)
foo <- nlminb(l$par, l$fn, l$gr, l$he)
pred5 = MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "gerlang", link_type_base = "softplus", link_type_covs = "softplus", covs_bin = T, transient_dist_method = "pade", warping=T)
err_mjp5 = rps_vectors(m, pred5, obs)
mean(err_mjp5) #0.4415872




D <- data.frame(s1 = d$s1, s2 = d$s2, err_mjp = rps_vectors(m, pred4, obs), err_naive = rps_vectors(m, pred2, obs))

library(dplyr)

err_summary <- D %>%
  group_by(s1, s2) %>%
  summarise(
    n = n(),
    mean_err_mjp = mean(err_mjp, na.rm = TRUE),
    mean_err_naive = mean(err_naive, na.rm = TRUE),
    diff = mean_err_naive - mean_err_mjp,
    .groups = 'drop'
  )

library(ggplot2)

ggplot(err_summary, aes(x = factor(s1), y = factor(s2), fill = diff)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(x = "s1", y = "s2", fill = "Naive - MJP") +
  theme_minimal()


idx = which(d$s1 == 4 & d$s2 == 5)
pred4[idx,]
pred2[idx,]



### TEST OF ERROR CONTRIBUTIONS

#Best MJP:
data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = as.integer(T))
beta = runif(m*(m-1)/2 + length(exo.cols),0,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = "FUNCS_MJP_with_TMB")

foo <- nlminb(l$par, l$fn, l$gr, l$he)
pred1 = MJP_predict(m = m, s1 = d$s1, u = d$t, pars = c(foo$par[1:(m*(m-1)/2)],log(exp(c(1,1,1,2))-1) , tail(foo$par,length(exo.cols)) ) , z = as.matrix(d[,exo.cols]), generator = "free_upper_tri", link_type_base = "softplus", link_type_covs = "softplus", covs_bin = T, transient_dist_method = "pade", eps = 2^-52,
                    warping=T)
err_mjp = rps_vectors(m, pred1, obs)
mean(err_mjp)









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

library(Matrix)

s1 = 1; s2 = 4; tau = 5; dt = 0.01; mu = rep(0, m); p1 = matrix(rep(0,m),ncol = m); p1[s1] = 1; p2 = p1; 
A = construct.A1(m, rep(0.1,(m-1)*m/2)); tpm=as.matrix(expm(A*dt))

for(i in 1:(tau/dt)){
  p1 = p1 %*% tpm
  p2 = p2 + p2 %*% A * dt
  mu = mu + p2
}
mu = mu*dt



mu


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

