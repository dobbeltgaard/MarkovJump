
rm(list = ls()); gc()

d = read.csv("defect_data.csv")
states <- c(1,2,3,4,5)
m <- length(states) 
track <- unique(d$Track)
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

library(MASS); library(Rcpp); library(RcppEigen); library(ranger); library(nnet); 
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
smart.empirical.pred = function(s1, s2){
  trans_counts <- table(s1, s2)
  trans_matrix <- prop.table(trans_counts, margin = 1)  # normalize by row
  return(as.matrix(trans_matrix))
}
randomwalk.empirical <- function(s1, s2) {
  trans_counts <- table(s1, s2)
  states <- sort(unique(c(s1, s2)))
  m <- length(states)
  trans_matrix <- matrix(0, m, m)
  
  for (i in 1:(m-1)) {
    stay <- ifelse(i %in% rownames(trans_counts), trans_counts[i, i], 0)
    jump <- ifelse(i %in% rownames(trans_counts), trans_counts[i, i+1], 0)
    total <- stay + jump
    if (total > 0) {
      trans_matrix[i, i] <- stay / total
      trans_matrix[i, i+1] <- jump / total
    } else {
      trans_matrix[i, i] <- 1
    }
  }
  trans_matrix[m, m] <- 1
  rownames(trans_matrix) <- colnames(trans_matrix) <- states
  return(trans_matrix)
}
majority.pred <- function(s1, s2) {
  majority_class <- names(sort(table(s2), decreasing = TRUE))[1]
  states <- sort(unique(c(s1, s2)))
  m <- length(states)
  trans_matrix <- matrix(0, m, m)
  col_idx <- which(states == majority_class)
  trans_matrix[, col_idx] <- 1
  rownames(trans_matrix) <- colnames(trans_matrix) <- states
  return(trans_matrix)
}


### INITIALIZE k-fold ###
set.seed(1); k= 5; bin_km= 0.1
d$pos_bin=floor(d$pos / bin_km); d$group_id=interaction(d$Track0, d$pos_bin, drop = TRUE); 
groups=levels(d$group_id); G=length(groups)
fold_id=sample(rep(1:k, length.out = G)); names(fold_id)=groups; fold_size=integer(k)

PREDS = list(); PARS = list(); CONV <- list(); GRADNORM <- list()
baselines = c("DTMC","majority","randomwalk","uniform", "persistence", "empirical_dist", "empirical_dist_corr","empirical_dist_smart" , 
              "olr", "olr_cov","opr", "opr_cov","ocllr","ocllr_cov", "obs")

if (!dir.exists("estimates_V3")) dir.create("estimates_V3")
if (!dir.exists("predictions_V3")) dir.create("predictions_V3")
if (!dir.exists("diagnostics_V3")) dir.create("diagnostics_V3")

pred = NULL
for(i in 1:5){
  start_time <- Sys.time()
  
  #Estimation and test set splits
  test_grps=names(fold_id)[fold_id == i]
  pred.idx=which(d$group_id %in% test_grps)
  d.test=d[pred.idx, , drop = FALSE]; d.train=d[-pred.idx, , drop = FALSE]
  fold_size[i] <- nrow(d.test)
  
  obs = make_Ptu_obs(m, d.test$s2)
  count = 0
  
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
  
  #Estimate other models
  drf.train = d.train; drf.test = d.test
  
  for(j in baselines){
    npar = NA
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
    if(j == "majority"){
      empi_foo = majority.pred(d.train$s1,d.train$s2)
      pred = empi_foo[d.test$s1, ]
    }
    if(j == "randomwalk.empirical"){
      empi_foo = randomwalk.empirical(d.train$s1,d.train$s2)
      pred = empi_foo[d.test$s1, ]
    }
    
    if(j == "olr"){ 
      if(i == 1){ PARS[[j]] = olr$coefficients} else { PARS[[j]] = rbind(PARS[[j]],olr$coefficients)} #store model estimates
      npar = c(olr$coefficients, olr$zeta)
      pred = predict(olr, newdata = d.test, type = "p");
    }
    if(j == "olr_cov"){ 
      if(i == 1){ PARS[[j]] = olr_cov$coefficients} else { PARS[[j]] = rbind(PARS[[j]],olr_cov$coefficients)} #store model estimates
      npar = c(olr_cov$coefficients, olr_cov$zeta)
      pred = predict(olr_cov, newdata = d.test, type = "p"); }
    if(j == "opr"){ 
      if(i == 1){ PARS[[j]] = opr$coefficients} else { PARS[[j]] = rbind(PARS[[j]],opr$coefficients)} #store model estimates
      npar = c(opr$coefficients, opr$zeta)
      pred = predict(opr, newdata = d.test, type = "p"); }
    if(j == "opr_cov"){ 
      if(i == 1){ PARS[[j]] = opr_cov$coefficients} else { PARS[[j]] = rbind(PARS[[j]],opr_cov$coefficients)} #store model estimates
      npar = c(opr_cov$coefficients, opr_cov$zeta)
      pred = predict(opr_cov, newdata = d.test, type = "p"); }
    if(j == "ocllr"){ 
      if(i == 1){ PARS[[j]] = ocllr$coefficients} else { PARS[[j]] = rbind(PARS[[j]],ocllr$coefficients)} #store model estimates
      npar = c(ocllr$coefficients, ocllr$zeta)
      pred = predict(ocllr, newdata = d.test, type = "p"); }
    if(j == "ocllr_cov"){ 
      if(i == 1){ PARS[[j]] = ocllr_cov$coefficients} else { PARS[[j]] = rbind(PARS[[j]],ocllr_cov$coefficients)} #store model estimates
      npar = c(ocllr_cov$coefficients, ocllr_cov$zeta)
      pred = predict(ocllr_cov, newdata = d.test, type = "p"); }
    if(j == "DTMC"){
      breaks <- c(seq(0,3.5,0.05),Inf)#c(0, 0.25, 0.5, 1, 2, Inf) 
      bin_train <- cut(d.train$t, breaks=breaks, right=FALSE)
      bin_test  <- cut(d.test$t,  breaks=breaks, right=FALSE)
      
      P_by_bin <- lapply(levels(bin_train), function(b) {
        idx <- bin_train == b
        N <- with(d.train[idx,], table(factor(s1, levels=1:m), factor(s2, levels=1:m)))
        alpha <- 0.1
        sweep(N + alpha, 1, rowSums(N + alpha), "/")
      })
      names(P_by_bin) <- levels(bin_train)
      
      pred_test <- t(sapply(seq_len(nrow(d.test)), function(i) {
        P_by_bin[[as.character(bin_test[i])]][as.integer(d.test$s1[i]), ]
      }))
      pred = pred_test
      
      N_by_bin <- lapply(levels(bin_train), function(b) {
        idx <- bin_train == b
        with(d.train[idx,], table(factor(s1, levels=1:m), factor(s2, levels=1:m)))
      })
      names(N_by_bin) <- levels(bin_train)
      bins_with_data <- sapply(N_by_bin, function(N) sum(N) > 0)
      L_eff <- sum(bins_with_data)
      
      npar <- L_eff * m * (m - 1)
    }
    if(j == "obs"){ pred = obs; }
    if(i == 1){ PREDS[[j]] = pred} else { PREDS[[j]] = rbind(PREDS[[j]],pred)} 
    write.csv(pred, file = paste0("predictions_V3/", j, "_fold_", i, ".csv"), row.names = FALSE)
    write.csv(data.frame(par = npar), file = paste0("estimates_V3/", j, "_fold_", i, ".csv"), row.names = FALSE)
    
  }
  end_time <- Sys.time(); runtime <- end_time - start_time
  cat("Fold ", i, "runtime:", runtime, "\n") #approx 1.1 second per fold
}
print(fold_size) #703 725 751 752 727

