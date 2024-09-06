

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

### Cross validation errors ###
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

### Estimated coefficients ###
library(tidyr)
library(dplyr)
par_list = readRDS("results/estimated_model_pars.Rdata")
mat = par_list[["free_upper_tri_TRUE_all"]]
mat = mat[, (ncol(mat)-4):ncol(mat)]

M = NULL
for(i in names(par_list)){
  if(grepl("TRUE",i) & !grepl("brier",i) & !grepl("rps",i)){
  mat = par_list[[i]]
  mat = mat[, (ncol(mat)-4):ncol(mat)]
  
  if(is.null(M)){ M = mat} else {M = rbind(M, mat)}
  }
}
df <- as.data.frame(M)
colnames(df) <- c("Var1", "Var2", "Var3", "Var4", "Var5")
df_long <- pivot_longer(df, cols = everything(), names_to = "Variable", values_to = "Value")
ggplot(df_long, aes(x = Variable, y = Value)) +
  geom_boxplot() +
  theme_minimal()






####################
### THE TRIPTYCH ###
####################
#i.e. binary evaluation

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

library(ggplot2)
m = 5; int1 = m-1; int2 = m; 

### Reliability diagrams ###
source("binary_eval_funcs/reliability_diag.R")
pdf(file = "figures/reliability_diagram_empirical_corr.pdf",width = 5, height = 4) 
par(mar = c(4, 4, 1, 1))
rd1 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["empirical_dist_corr"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()
pdf(file = "figures/reliability_diagram_olr_cov.pdf",width = 5, height = 4) 
par(mar = c(4, 4, 1, 1))
rd2 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["olr_cov"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()
pdf(file = "figures/reliability_diagram_mjp_free.pdf",width = 5, height = 4) 
par(mar = c(4, 4, 1, 1))
rd3 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["free_upper_tri_TRUE_all"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()
pdf(file = "figures/reliability_diagram_ensemble.pdf",width = 5, height = 4) 
par(mar = c(4, 4, 1, 1))
rd4 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["ensemble"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()

### ROC curves ###
library("ROSE")
rc0 = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["empirical_dist"]]), n.thresholds = 200)
rc1 = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["empirical_dist_corr"]]),add.roc = T, n.thresholds = 200)
rc2 = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["olr_cov"]]), add.roc = T, n.thresholds = 200)
rc3 = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["free_upper_tri_TRUE_all"]]), add.roc = T, n.thresholds = 200)
rc4 = roc.curve(zero_bin_obs(int1,int2,pred_list[["obs"]]), bin_predictions(int1,int2,pred_list[["ensemble"]]), add.roc = T, n.thresholds = 200)
nams = c(paste("Empirical dist. = ", format(round(rc1$auc, 3), nsmall = 2)), 
         paste("OLR with covs. = ", format(round(rc2$auc, 3), nsmall = 2)),
         paste("MJP with covs. = ", format(round(rc3$auc, 3), nsmall = 2)),
         paste("Ensemble = ", format(round(rc4$auc, 3), nsmall = 2))
        )
RC1 = data.frame(rc1$false.positive.rate,rc1$true.positive.rate, rep(nams[1], length(rc1$true.positive.rate)) ); colnames(RC1) = c("false.positive.rate", "true.positive.rate", "group");
RC2 = data.frame(rc2$false.positive.rate,rc2$true.positive.rate, rep(nams[2], length(rc2$true.positive.rate)) ); colnames(RC2) = c("false.positive.rate", "true.positive.rate", "group");
RC3 = data.frame(rc3$false.positive.rate,rc3$true.positive.rate, rep(nams[3], length(rc3$true.positive.rate)) ); colnames(RC3) = c("false.positive.rate", "true.positive.rate", "group");
RC4 = data.frame(rc4$false.positive.rate,rc4$true.positive.rate, rep(nams[4], length(rc4$true.positive.rate)) ); colnames(RC4) = c("false.positive.rate", "true.positive.rate", "group");
combined_data <- rbind(RC1, RC2, RC3, RC4)
combined_data$group = factor(combined_data$group, levels = c(nams[1],nams[2], nams[3], nams[4]))
text.size = 13
pdf(file = "figures/ROCs.pdf",width = 6, height = 4.8) 
ggplot(combined_data, aes(x = false.positive.rate, y = true.positive.rate, color = group)) +
  geom_line(linewidth = 0.75) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "Group") +
  theme_minimal() +
  theme(text = element_text(size = text.size, family = "serif"), 
                        panel.background = element_rect(fill = "white", color = "black"), 
                        panel.grid.minor = element_line(color = "lightgray"),
                        legend.position = c(0.8, 0.6),  # Moves the legend closer to the bottom
                        legend.direction = "horizontal",
                        legend.background = element_rect(fill = "white", color = "black", size = 0.25, linetype = "solid"),  # Black border and white background
  ) + 
  guides(color = guide_legend(title = "", ncol = 1))
dev.off()

### Murphy diagrams ###
source("binary_eval_funcs/murphy_diag.R")
md1 = murphydiagram(f1 = bin_predictions(int1,int2,pred_list[["empirical_dist_corr"]]), f2 = bin_predictions(int1,int2,pred_list[["olr_cov"]]),  y = bin_predictions(int1,int2,pred_list[["obs"]]))
md2 = murphydiagram(f1 = bin_predictions(int1,int2,pred_list[["free_upper_tri_TRUE_all"]]), f2 = bin_predictions(int1,int2,pred_list[["ensemble"]]),  y = bin_predictions(int1,int2,pred_list[["obs"]]))
a1 = sum(md1$y[c(-1),1] * diff(md1$tsep))
a2 =sum(md1$y[c(-1),2] * diff(md1$tsep))
a3 =sum(md2$y[c(-1),1] * diff(md2$tsep))
a4 =sum(md2$y[c(-1),2] * diff(md2$tsep))
nams = c(paste("Empirical dist. = ", format(round(a1, 3), nsmall = 2)), 
         paste("OLR with covs. = ", format(round(a2, 3), nsmall = 2)),
         paste("MJP with covs. = ", format(round(a3, 3), nsmall = 2)),
         paste("Ensemble = ", format(round(a4, 3), nsmall = 2))
)
MD1 = data.frame(md1$tsep,md1$y[,1], rep(nams[1], length(md1$y[,1])) ); colnames(MD1) = c("par", "val", "group");
MD2 = data.frame(md1$tsep,md1$y[,2], rep(nams[2], length(md1$y[,2])) ); colnames(MD2) = c("par", "val", "group");
MD3 = data.frame(md2$tsep,md2$y[,1], rep(nams[3], length(md2$y[,1])) ); colnames(MD3) = c("par", "val", "group");
MD4 = data.frame(md2$tsep,md2$y[,2], rep(nams[4], length(md2$y[,2])) ); colnames(MD4) = c("par", "val", "group");
combined_data <- rbind(MD1, MD2, MD3, MD4)
combined_data$group = factor(combined_data$group, levels = c(nams[1],nams[2], nams[3], nams[4]))
text.size = 13
pdf(file = "figures/muprhy.pdf",width = 6, height = 4.8) 
ggplot(combined_data, aes(x = par, y = val, color = group)) +
  geom_line(linewidth = 0.75) +
  labs(x = expression(omega), y = expression(paste("Empirical score, ", S(omega),"")), color = "Group") +
  theme_minimal() + 
  theme(text = element_text(size = text.size, family = "serif"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        panel.grid.minor = element_line(color = "lightgray"),
        legend.position = c(0.8, 0.6),  # Moves the legend closer to the bottom
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color = "black", size = 0.25, linetype = "solid"),  # Black border and white background
  ) + 
  guides(color = guide_legend(title = "", ncol = 1))
dev.off()
