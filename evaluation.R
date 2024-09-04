

rm(list = ls()) #clear memory
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")
pred_list = readRDS("results/model_predictions.Rdata")

pred_list[["ensemble"]] = (pred_list[["olr_cov"]] + pred_list[["free_upper_tri_TRUE_all"]] )/2#+ pred_list[["empirical_dist_corr"]])/3
pred_list[["ensemble2"]] = (pred_list[["olr_cov"]] + pred_list[["free_upper_tri_TRUE_all"]] + pred_list[["empirical_dist_corr"]])/3


#logs = read.csv("results/log_score.csv")
#rps = read.csv("results/rps_score.csv")
#brier = read.csv("results/brier_score.csv")

m = 5; k = 10; 
log_err = matrix(NA, ncol = k, nrow = length(names(pred_list)))
RPS_err = matrix(NA, ncol = k, nrow = length(names(pred_list)))
Brier_e = matrix(NA, ncol = k, nrow = length(names(pred_list)))
for(i in 1:k){
  count = 0; 
  for(j in names(pred_list)){
    count = count + 1;
    
    pred = pred_list[[j]]
    obs = pred_list[["obs"]]
    
    log_err[count, i] = mean(logscore_vectors(m, pred, obs))
    RPS_err[count, i] = mean(rps_vectors(m, pred, obs))
    Brier_e[count, i] = mean(BrierScore_vectors(m, pred, obs))
    
  }
}




rownames(log_err) = names(pred_list)
rownames(RPS_err) = names(pred_list)
rownames(Brier_e) = names(pred_list)

errs = cbind(rowMeans(log_err),rowMeans(Brier_e),rowMeans(RPS_err))
errs

#we need to save predictions differently
# err_tab = cbind(rowMeans(logs[,2:11]),rowMeans(rps[,2:11]),rowMeans(brier[,2:11]))
# rownames(err_tab) = logs$X
# err_tab








#How good are the respective models to forecast 0 or 1 events?
bin_predictions = function(col1, col2, preds){ #convert multicategorial predictions into binary for event in col1 to col2
  if(col1 == col2){
    res = preds[,col1]
  } else if(col1 < col2){
    res = rowSums(preds[,col1:col2])
  }
  res[res > 1] = 1;
  res[res < 0] = 0;
  return(res)
}

zero_bin_obs = function(col1, col2, obs){
  idx = apply(obs, 1, function(row) which(row == 1))
  return(idx >= col1 & idx <= col2)
}


#you need to add observations to predlist

####################
### THE TRIPTYCH ###
####################

### Reliability diagrams ###
source("binary_eval_funcs/reliability_diag.R")
m = 5; int1 = m-1; int2 = m; 

rd = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["gerlang_TRUE_all"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T)
rd = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["gerlang_relax_TRUE_all"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T)
rd = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["free_upper_tri_TRUE_all"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T)
rd = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["olr_cov"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T)
rd = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["ensemble"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T)


### ROC curves ###
library("ROSE")
rc = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["empirical_dist_corr"]]))
roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["uniform"]]), add.roc = T)
roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["free_upper_tri_TRUE_all"]]), add.roc = T)
roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["olr_cov"]]), add.roc = T)
roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["ensemble"]]), add.roc = T)
roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["ensemble2"]]), add.roc = T)

plot(rc$false.positive.rate, rc$true.positive.rate, type = "l")


### Murphy diagrams ###
source("binary_eval_funcs/murphy_diag.R")
md = murphydiagram(f1 = bin_predictions(int1,int2,pred_list[["free_upper_tri_TRUE_all"]]), f2 = bin_predictions(int1,int2,pred_list[["olr_cov"]]),  y = bin_predictions(int1,int2,pred_list[["obs"]]))
murphydiagram(f1 = bin_predictions(int1,int2,pred_list[["olr_cov"]]), f2 = bin_predictions(int1,int2,pred_list[["ensemble"]]),  y = bin_predictions(int1,int2,pred_list[["obs"]]))






















