

rm(list = ls()) #clear memory
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")
pred_list = readRDS("results/model_predictions_V2.Rdata")
pred_list[["ensemble"]] = (pred_list[["olr_cov"]] + pred_list[["free_upper_tri_TRUE_softplus_softplus_all"]] )/2#+ pred_list[["empirical_dist_corr"]])/3
pred_list[["ensemble2"]] = (pred_list[["olr_cov"]] + pred_list[["free_upper_tri_TRUE_softplus_softplus_all"]] + pred_list[["empirical_dist_corr"]])/3
#logs = read.csv("results/log_score.csv")
#rps = read.csv("results/rps_score.csv")
#brier = read.csv("results/brier_score.csv")




#######################
### Forecast Scores ###
#######################
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
    log_err[count, i] = -mean(logscore_vectors(m, pred, obs))
    RPS_err[count, i] = mean(rps_vectors(m, pred, obs))
    Brier_e[count, i] = mean(BrierScore_vectors(m, pred, obs))
  }
}
rownames(log_err) = names(pred_list)
rownames(RPS_err) = names(pred_list)
rownames(Brier_e) = names(pred_list)
errs = cbind(rowMeans(RPS_err),rowMeans(log_err),rowMeans(Brier_e))

#Errors in link function combinations
link_combinations <- gsub(".*?(exp|softplus|square)_(exp|softplus|square)_.*", "\\1_\\2", rownames(errs[1:54,]))
rps_values <- errs[1:54, 1]
average_rps <- tapply(rps_values, link_combinations, mean)
print(average_rps)


namfoo = c("softplus_softplus_all", "uniform", "dist_corr", "olr", "opr", "ocllr", "ensemble")
idx = rowSums(sapply(namfoo, FUN = grepl, x = rownames(errs))) > 0
foo = errs[idx, ]
nams =
  c("uniform", "empirical_dist_corr", 
    "olr", "olr_cov", "opr", "opr_cov", "ocllr", "ocllr_cov", 
    "gerlang_FALSE_softplus_softplus_all", "gerlang_TRUE_softplus_softplus_all", "gerlang_relax_FALSE_softplus_softplus_all", "gerlang_relax_TRUE_softplus_softplus_all", "free_upper_tri_FALSE_softplus_softplus_all", "free_upper_tri_TRUE_softplus_softplus_all",
    "ensemble", "ensemble2")
nams_idx = rep(NA, length(nams));for(i in 1:length(nams)){nams_idx[i] = which(nams[i] == rownames(foo))}
library(kableExtra)
naive = expand.grid(
  Covariates = c("No"), #covariates
  Model = c("Naïve predictors"), #näive,olr,mjp 
  param = c("Uniform distribution", "Empirical distribution")
)
regression = expand.grid(
  Covariates = c("No","Yes"), #covariates
  Model = c("Ordered multinomial regression"), #näive,olr,mjp 
  param = c("Logistic link", "Probit link", "Cloglog link")
)
mjp = expand.grid(
  Covariates = c("No","Yes"), #covariates
  Model = c("Markov jump process"), #näive,olr,mjp 
  param = c("Generalized Erlang, $\\bm{A}'$", "Parameterized upper, $\\bm{A}''$", "Free upper, $\\bm{A}'''$")
)
ensemble = expand.grid(
  Covariates = c("Yes"), #covariates
  Model = c("Ensemble"), #näive,olr,mjp 
  param = c("MJP + OLR", "Naïve + MJP + OLR")
)
collapse_rows_dt = rbind(naive, regression, mjp,ensemble)
collapse_rows_dt <- collapse_rows_dt[c("Model", "param", "Covariates")]

par_list  = pred_list
npars = rep(0,length(nams)); count = 0;
for(i in nams[1:(length(nams)-2)]){
  count = count + 1
  if(count > 3){npars[count] = (NCOL(par_list[[i]])); }
  
}
npars[length(nams)-1] =  NCOL(par_list[["free_upper_tri_TRUE_softplus_softplus_all"]]) + NCOL(par_list[["olr_cov"]])
npars[length(nams)] =  NCOL(par_list[["free_upper_tri_TRUE_softplus_softplus_all"]]) + NCOL(par_list[["olr_cov"]])
collapse_rows_dt$npars = npars
collapse_rows_dt$RPS = foo[nams_idx, 1]
collapse_rows_dt$LogS = foo[nams_idx, 2]
collapse_rows_dt$Brier = foo[nams_idx, 3]
colnames(collapse_rows_dt) = c("Method", "Model", "Covariates", "\\# Parameters", "RPS", "Log S.", "Brier S.")
row_group_label_fonts <- list(
  list(bold = T, italic = F),
  list(bold = F, italic = F)
)
kableExtra::kbl(collapse_rows_dt,booktabs = T, align = c("l","l","c","c","c","c","c"), linesep = '', format = "latex",escape = FALSE, digits = 3) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts)


##############################
### Estimated coefficients ###
##############################
library(kableExtra)
par_list = readRDS("results/estimated_model_pars.Rdata")
nams = c("olr_cov", "opr_cov", "ocllr_cov", "gerlang_TRUE_all", "gerlang_relax_TRUE_all", "free_upper_tri_TRUE_all")
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
count = 0; k = 10;  
for(i in nams){
  mat = par_list[[i]]
  mat = mat[, (ncol(mat)-4):ncol(mat)]
  if(count == 0){ estim = apply(mat, MARGIN = 2, FUN = quantile, probs = c(0, 0.5, 1)); NAMS = rep(i, 3);} else { estim = rbind(estim, apply(mat, MARGIN = 2, FUN = quantile, probs = c(0.05, 0.5, 0.95))); NAMS = c(NAMS, rep(i,3))}
  count = count + 1
}
df = as.data.frame(round(estim,digits=2)); 
df$X = NAMS
df$quantiles = rep(c("0%", "50%", "100%"), length(nams))
get_quantile_strings = function(df, colname){
  quants = paste(df[df$quantiles == "0%",colname], df[df$quantiles == "100%",colname], sep = ", ")
  sep_brackets1 = rep("[", length(quants)); sep_brackets2 = rep("]", length(quants)); 
  return(#paste(df[df$quantiles == "50%",colname], 
               paste(paste(sep_brackets1, quants, sep = ""), sep_brackets2, sep = "")
               #)
  )
}
foo = sapply(exo.cols,FUN = get_quantile_strings, df = df)
row.names(foo) = nams
colnames(foo) = c("Tonnage", "Line speed", "Rail profile", "Steel hardness", "Curvature")
foo = as.data.frame((foo))
foo$Model = NA
foo[grepl("cov", rownames(foo) ), "Model"] = "Ordered multinomial regression"
foo[!grepl("cov", rownames(foo) ), "Model"] = "Markov jump process"
foo$param = NA
foo[grepl("cov", rownames(foo) ), "param"] = c("Logistic link", "Probit link", "Cloglog link")
foo[!grepl("cov", rownames(foo) ), "param"] = c("Generalized Erlang, $\\bm{A}'$", "Parameterized upper, $\\bm{A}''$", "Free upper, $\\bm{A}'''$")

foo = foo[, c("Model", "param","Tonnage", "Line speed", "Rail profile", "Steel hardness", "Curvature")]
rownames(foo) = 1:NROW(foo)
colnames(foo) = c("param", "Model","Tonnage", "Line speed", "Rail profile", "Steel hardness", "Curvature")
row_group_label_fonts <- list(
  list(bold = T, italic = F),
  list(bold = F, italic = F)
)
kableExtra::kbl(foo,booktabs = T, align = c("l","l","c","c","c","c","c"), linesep = '', format = "latex",escape = FALSE) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(1:2, latex_hline = 'major',row_group_label_position = 'stack',row_group_label_fonts = row_group_label_fonts)


###############################################
### plot of transition probability matrices ###
###############################################
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")
library(Matrix)
library(ggplot2)
library(ggpubr)
library(MASS)
library(reshape2)

states = c("3", "2B", "2A", "1", "0")
par_list = readRDS("results/estimated_model_pars.Rdata")
d = read.csv("defect_data.csv")
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
z = as.matrix(d[,exo.cols])
idx = 1 #sample(1:NROW(d), size = 1)
text.size <- 11
ndays = 365*8

#### gerlang ###
nam = "gerlang_TRUE_all"
m=5; 
npars = m-1
lambda_base = exp(par_list[[nam]][1,1:npars])
covs = par_list[[nam]][1,(npars+1):(npars + 5)]
lambda = lambda_base * exp(sum(covs * z[idx,]))
sol1 = matrix(NA, nrow = ndays, ncol = m)
count = 0
for(t in (1:ndays)/365){
  count = count + 1
  A1 = as.matrix(Matrix::expm(make_A1(m, lambda)*t))
  sol1[count, ] = c(1,0,0,0,0) %*% A1
}
dpp <- as.data.frame(sol1)
colnames(dpp) <- states
dpp$time <- 1:ndays
dpp_long <- tidyr::gather(dpp, key = "Column", value = "Probability", -time)
dpp_long$Column <- factor(dpp_long$Column, levels = c("3", "2B", "2A", "1", "0")) 
p1 <- ggplot(data = dpp_long, aes(x = time/365, y = Probability, color = Column)) +
  geom_line(size = 0.75) +
  theme(
    text = element_text(size = text.size, family = "serif"),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_line(color = "lightgray"),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal", # Set legend direction to horizontal
    legend.background = element_rect(fill = "transparent", color = NA), # Set transparent background
    legend.key = element_rect(fill = "transparent", color = NA) # Set transparent background for legend key
  ) + xlab("Time [years]") + ylab("Probability") + labs(color = "Classes")
colnames(A1) = 1:5
rownames(A1) = 5:1
longData<-melt(A1)
longData<-longData[longData$value!=0,]
p11 = ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill = value)) + 
  geom_text(aes(label = sprintf("%.5f", value)), color = "white", size = 3, family = "serif") +  
  scale_fill_gradient(low = "grey60", high = "black", guide = "none") + 
  scale_y_discrete(limits = c("0", "1", "2A", "2B", "3")) +
  scale_x_discrete(limits = c("3", "2B", "2A", "1", "0")) +
  theme(
    axis.text.x = element_text(size = 9, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 11),
    text = element_text(size = text.size, family = "serif"),  # Ensure overall text is serif
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_line(color = "lightgray")
  ) + 
  xlab("To class") + 
  ylab("From class")

p11

#### reparameterized upper ###
nam = "gerlang_relax_TRUE_all"
m=5; 
npars = m-1
lambda_base = exp(par_list[[nam]][1,1:npars])
covs = par_list[[nam]][1,(npars+1):(npars + 5)]
lambda = lambda_base * exp(sum(covs * z[idx,]))
sol1 = matrix(NA, nrow = ndays, ncol = m)
count = 0
for(t in (1:ndays)/365){
  count = count + 1
  A1 = as.matrix(Matrix::expm(make_A2(m, lambda)*t))
  sol1[count, ] = c(1,0,0,0,0) %*% A1
}
dpp <- as.data.frame(sol1)
colnames(dpp) <- states
dpp$time <- 1:ndays
dpp_long <- tidyr::gather(dpp, key = "Column", value = "Probability", -time)
dpp_long$Column <- factor(dpp_long$Column, levels = c("3", "2B", "2A", "1", "0")) 
p2 <- ggplot(data = dpp_long, aes(x = time/365, y = Probability, color = Column)) +
  geom_line(size = 0.75) +
  theme(
    text = element_text(size = text.size, family = "serif"),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_line(color = "lightgray"),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal", # Set legend direction to horizontal
    legend.background = element_rect(fill = "transparent", color = NA), # Set transparent background
    legend.key = element_rect(fill = "transparent", color = NA) # Set transparent background for legend key
  ) + xlab("Time [years]") + ylab("Probability") + labs(color = "Classes")
colnames(A1) = 1:5
rownames(A1) = 5:1
longData<-melt(A1)
longData<-longData[longData$value!=0,]
p22 = ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill = value)) + 
  geom_text(aes(label = sprintf("%.5f", value)), color = "white", size = 3, family = "serif") +  
  scale_fill_gradient(low = "grey60", high = "black", guide = "none") + 
  scale_y_discrete(limits = c("0", "1", "2A", "2B", "3")) +
  scale_x_discrete(limits = c("3", "2B", "2A", "1", "0")) +
  theme(
    axis.text.x = element_text(size = 9, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 11),
    text = element_text(size = text.size, family = "serif"),  # Ensure overall text is serif
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_line(color = "lightgray")
  ) + 
  xlab("To class") + 
  ylab("From class")


#### free upper ###
nam = "free_upper_tri_TRUE_all"
m=5; 
npars = m*(m-1)/2
lambda_base = exp(par_list[[nam]][1,1:npars])
covs = par_list[[nam]][1,(npars+1):(npars + 5)]
lambda = lambda_base * exp(sum(covs * z[idx,]))
sol1 = matrix(NA, nrow = ndays, ncol = m)
count = 0
for(t in (1:ndays)/365){
  count = count + 1
  A1 = as.matrix(Matrix::expm(make_A3(m, lambda)*t))
  sol1[count, ] = c(1,0,0,0,0) %*% A1
}
dpp <- as.data.frame(sol1)
colnames(dpp) <- states
dpp$time <- 1:ndays
dpp_long <- tidyr::gather(dpp, key = "Column", value = "Probability", -time)
dpp_long$Column <- factor(dpp_long$Column, levels = c("3", "2B", "2A", "1", "0")) 
p3 <- ggplot(data = dpp_long, aes(x = time/365, y = Probability, color = Column)) +
  geom_line(size = 0.75) +
  theme(
    text = element_text(size = text.size, family = "serif"),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_line(color = "lightgray"),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal", # Set legend direction to horizontal
    legend.background = element_rect(fill = "transparent", color = NA), # Set transparent background
    legend.key = element_rect(fill = "transparent", color = NA) # Set transparent background for legend key
  ) + xlab("Time [years]") + ylab("Probability") + labs(color = "Classes")
colnames(A1) = 1:5
rownames(A1) = 5:1
longData<-melt(A1)
longData<-longData[longData$value!=0,]
p33 = ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill = value)) + 
  geom_text(aes(label = sprintf("%.5f", value)), color = "white", size = 3, family = "serif") +  
  scale_fill_gradient(low = "grey60", high = "black", guide = "none") + 
  scale_y_discrete(limits = c("0", "1", "2A", "2B", "3")) +
  scale_x_discrete(limits = c("3", "2B", "2A", "1", "0")) +
  theme(
    axis.text.x = element_text(size = 9, angle = 0, vjust = 0.3),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 11),
    text = element_text(size = text.size, family = "serif"),  # Ensure overall text is serif
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_line(color = "lightgray")
  ) + 
  xlab("To class") + 
  ylab("From class")


pdf(file = "figures/trans_dist_gerlang.pdf",width = 4, height = 3) 
p1
dev.off()
pdf(file = "figures/trans_dist_param_upper.pdf",width = 4, height = 3) 
p2
dev.off()
pdf(file = "figures/trans_dist_free_upper.pdf",width = 4, height = 3) 
p3
dev.off()

pdf(file = "figures/tpm_gerlang.pdf",width = 4, height = 3) 
p11
dev.off()
pdf(file = "figures/tpm_param_upper.pdf",width = 4, height = 3) 
p22
dev.off()
pdf(file = "figures/tpm_free_upper.pdf",width = 4, height = 3) 
p33
dev.off()


### OLR logistic ###
olr_cov = MASS::polr(as.factor(s2) ~ as.factor(s1) + t + MBT.norm + speed.norm + profil.norm + steel.norm + invRad.norm, data = d, method = "logistic")
d.test = d[idx, ]
d.test$s1 = 1
sol1 = matrix(NA, nrow = ndays, ncol = m)
count = 0
for(t in (1:ndays)/365){
  count = count + 1
  d.test$t = t
  sol1[count, ] = predict(olr_cov, newdata = d.test, type = "p");
}
dpp <- as.data.frame(sol1)
colnames(dpp) <- states
dpp$time <- 1:ndays
dpp_long <- tidyr::gather(dpp, key = "Column", value = "Probability", -time)
dpp_long$Column <- factor(dpp_long$Column, levels = c("3", "2B", "2A", "1", "0")) 
p4 <- ggplot(data = dpp_long, aes(x = time/365, y = Probability, color = Column)) +
  geom_line(size = 0.75) +
  theme(
    text = element_text(size = text.size, family = "serif"),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_line(color = "lightgray"),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal", # Set legend direction to horizontal
    legend.background = element_rect(fill = "transparent", color = NA), # Set transparent background
    legend.key = element_rect(fill = "transparent", color = NA) # Set transparent background for legend key
  ) + xlab("Time [years]") + ylab("Probability") + labs(color = "Classes")

##########################################
### Multicategory reliability diagrams ###
##########################################
library(reshape2)
library(ggplot2)
library(cowplot)

forecast_category <- function(forecast, quantiles) {
  cumulative <- cumsum(forecast)
  res = sapply(quantiles, function(q) which(cumulative >= q)[1])
  return(res)
}
compute_reliability <- function(obs, z, q, qmin, qmax) {
  if(obs > z){
    p = 0
  } else if(obs < z){
    p = 1
  } else if(obs == z){
    if(qmax != qmin){
      p = (q-qmin) / (qmax - qmin) 
    } else{
      p = 0.5
    }
  }
  return(p)
}
multi.reliable = function(obs, y, quantiles){
  L = length(quantiles)
  forecast_quantile_matrix <- t(apply(y, 1, forecast_category, quantiles = quantiles))
  Cq = matrix(NA, nrow = NROW(y), ncol = L)
  for(i in 1:NROW(y)){
    Q = which(forecast_quantile_matrix[i,] == obs[i])
    if(length(Q) > 0) {
      qmin = quantiles[min(Q)]
      qmax = quantiles[max(Q)]
    } else {
      qmin = NA
      qmax = NA
    }
    for(j in 1:L){
      Cq[i,j] = compute_reliability(obs[i], forecast_quantile_matrix[i,j],quantiles[j], qmin = qmin, qmax = qmax)
    }
  }
  Cqave = colMeans(Cq)
  
  forecast_diff_matrix <- forecast_quantile_matrix - obs
  abs_error_matrix = abs(forecast_diff_matrix)
  ave_cat_err = mean(colMeans(abs_error_matrix))
  
  return(list(
    Cq = Cqave, 
    ave_cat_err = ave_cat_err
  ))
}
make_boots = function(obs, y, quantiles, iters){
  res = matrix(NA, nrow = iters, ncol = length(quantiles))
  err = matrix(NA, nrow = iters, ncol = length(quantiles))
  for(i in 1:iters){
    idx = sample(1:NROW(obs), replace = T, size = NROW(obs) )
    res[i,] = multi.reliable(obs[idx], y[idx, ], quantiles)$Cq
    forecast_quantile_matrix <- t(apply(y[idx, ], 1, forecast_category, quantiles = quantiles))
    forecast_diff_matrix <- forecast_quantile_matrix - obs[idx]
    abs_error_matrix = abs(forecast_diff_matrix)
    err[i, ] = colMeans(abs_error_matrix)
  }
  return(list(
    Cq = apply(res, 2, quantile,c(0.1,0.9)), 
    quan_cat_err = apply(err, 2, quantile,c(0.1,0.9)) 
    ) )
}


checkerboard_plot = function(obs, y, quantiles){
  forecast_quantile_matrix <- t(apply(y, 1, forecast_category, quantiles = quantiles))
  forecast_diff_matrix <- forecast_quantile_matrix - obs
  diff_df <- data.frame(forecast_diff_matrix)
  colnames(diff_df) = quantiles
  diff_df$Observation <- obs
  diff_df$ID <- 1:nrow(diff_df)  # Create an ID for each observation
  diff_df_long <- melt(diff_df, id.vars = c("ID", "Observation"), variable.name = "Quantile", value.name = "Forecast_Diff")
  text.size=10

  diff_df_long$Quantile <- as.numeric(as.character(diff_df_long$Quantile))
  
  p1 = ggplot(diff_df_long, aes(x = Quantile, y = Forecast_Diff)) +
    geom_tile(aes(fill = after_stat(count)), color = "white", stat = "bin2d", binwidth =  c(length(quantiles)/100, 1)) +
    scale_fill_gradient(low = "white", high = "black") +
    labs(x = "Forecast quantiles", y = "Category error") +
    scale_x_continuous(breaks = seq(min(diff_df_long$Quantile), max(diff_df_long$Quantile), by = 0.4)) + # Adjust 'by' as needed
    theme(
      legend.position = "none", 
      text = element_text(size = text.size, family = "serif"),  
      panel.background = element_rect(fill = NA, color = NA),  
      plot.background = element_rect(fill = NA, color = NA),  
      panel.grid.minor = element_blank(),  
      panel.grid.major = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
    )
  return(p1)
}


reliability_plot = function(Cq, CqCI, quantiles ){
  data <- data.frame(
    Time = quantiles,
    Value = r$Cq,
    Lower = rb$Cq[1,],
    Upper = rb$Cq[2,]
  )
  text.size = 11
  p2 = ggplot(data, aes(x = Time, y = Value)) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +  
    geom_line(color = "black") +                
    geom_point(color = "black", size = 0.5) + 
    #xlim(0,1) + 
    ylim(0,1) +
    #geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "black", alpha = 0.2) +  
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.05, color = "black") +
    labs(x = "Quantile of forecast distribution", 
         y = "Observations below forecast quantile") +
    theme(
      legend.position = "none", 
      text = element_text(size = text.size, family = "serif"),  
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.minor = element_blank(),  
      panel.grid.major = element_blank(),  
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_x_continuous(breaks = quantiles) + 
    annotate("text", x = 0.25, y = 1, 
             label = paste("Average category error =", format(round(r$ave_cat_err,2), nsmall = 2)), 
             size = 3.5, hjust = 0.5, color = "black", family = "serif") +  
    annotate("text", x = 0.25, y = 0.925, 
             label = paste0("CI 90% = [", format(round(rowMeans(rb$quan_cat_err)[1],2), nsmall = 2), ", ", format(round(rowMeans(rb$quan_cat_err)[2],2), nsmall = 2), "]"), 
             size = 3.5, hjust = 0.5, color = "black", family = "serif")
  return(p2)
}


pred_list = readRDS("results/model_predictions.Rdata")
pred_list[["ensemble"]] = (pred_list[["olr_cov"]] + pred_list[["free_upper_tri_TRUE_all"]] )/2
nams = c("empirical_dist_corr","olr_cov", "free_upper_tri_TRUE_all", "ensemble")
for(i in 1:length(nams)){
  obs = pred_list[["obs"]]
  OBS = apply(obs, 1, which.max)
  y = pred_list[[nams[i]]]
  quantiles = seq(0.05, 0.95, by = 0.1)
  iters = 200
  text.size = 11
  
  r = multi.reliable(OBS, y,quantiles)
  rb = make_boots(OBS, y, quantiles, iters)
  
  p1 = checkerboard_plot(OBS, y, quantiles)
  p2 = reliability_plot(r, rb, quantiles)
  p3 = ggdraw() + draw_plot(p2) + draw_plot(p1, x = 0.6, y = 0.15, width = .35, height = .5)
  
  nfil = paste0("figures/reliability_diag_", nams[i],".pdf")
  pdf(file = nfil,width = 4, height = 3) 
  print(p3)
  dev.off()
  graphics.off()
  
  print(nams[i])
}


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
pdf(file = "figures/reliability_diagram_empirical_corr.pdf",width = 4, height = 3.5) 
par(mar = c(4, 4, 1, 1))
rd1 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["empirical_dist_corr"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()
pdf(file = "figures/reliability_diagram_olr_cov.pdf",width = 4, height = 3.5) 
par(mar = c(4, 4, 1, 1))
rd2 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["olr_cov"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()
pdf(file = "figures/reliability_diagram_mjp_free.pdf",width = 4, height = 3.5) 
par(mar = c(4, 4, 1, 1))
rd3 = ReliabilityDiagram2(bin_predictions(int1,int2,pred_list[["free_upper_tri_TRUE_all"]]), zero_bin_obs(int1,int2,pred_list[["obs"]]), plot = T, plot.refin = F, attributes = T, bins = 10)
dev.off()
pdf(file = "figures/reliability_diagram_ensemble.pdf",width = 4, height = 3.5) 
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
nams = c(paste("Empi. dist. = ", format(round(rc1$auc, 3), nsmall = 2)), 
         paste("Cumu. link = ", format(round(rc2$auc, 3), nsmall = 2)),
         paste("MJP = ", format(round(rc3$auc, 3), nsmall = 2)),
         paste("Ensemble = ", format(round(rc4$auc, 3), nsmall = 2))
        )
RC1 = data.frame(rc1$false.positive.rate,rc1$true.positive.rate, rep(nams[1], length(rc1$true.positive.rate)) ); colnames(RC1) = c("false.positive.rate", "true.positive.rate", "group");
RC2 = data.frame(rc2$false.positive.rate,rc2$true.positive.rate, rep(nams[2], length(rc2$true.positive.rate)) ); colnames(RC2) = c("false.positive.rate", "true.positive.rate", "group");
RC3 = data.frame(rc3$false.positive.rate,rc3$true.positive.rate, rep(nams[3], length(rc3$true.positive.rate)) ); colnames(RC3) = c("false.positive.rate", "true.positive.rate", "group");
RC4 = data.frame(rc4$false.positive.rate,rc4$true.positive.rate, rep(nams[4], length(rc4$true.positive.rate)) ); colnames(RC4) = c("false.positive.rate", "true.positive.rate", "group");
combined_data <- rbind(RC1, RC2, RC3, RC4)
combined_data$group = factor(combined_data$group, levels = c(nams[1],nams[2], nams[3], nams[4]))
text.size = 11
pdf(file = "figures/ROCs.pdf",width = 4, height = 3) 
ggplot(combined_data, aes(x = false.positive.rate, y = true.positive.rate, color = group)) +
  geom_line(linewidth = 0.75) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "Group") +
  theme_minimal() +
  theme(text = element_text(size = text.size, family = "serif"), 
                        panel.background = element_rect(fill = "white", color = "black"), 
                        panel.grid.minor = element_line(color = "lightgray"),
                        legend.position = c(0.75, 0.25),  # Moves the legend closer to the bottom
                        legend.direction = "horizontal",
                        legend.background = element_rect(fill = "transparent", color = NA),  # Set transparent background for legend
                        legend.key = element_rect(fill = "transparent", color = NA)
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
nams = c(paste("Empi. dist. = ", format(round(a1, 3), nsmall = 2)), 
         paste("Cumu. link = ", format(round(a2, 3), nsmall = 2)),
         paste("MJP = ", format(round(a3, 3), nsmall = 2)),
         paste("Ensemble = ", format(round(a4, 3), nsmall = 2))
)
MD1 = data.frame(md1$tsep,md1$y[,1], rep(nams[1], length(md1$y[,1])) ); colnames(MD1) = c("par", "val", "group");
MD2 = data.frame(md1$tsep,md1$y[,2], rep(nams[2], length(md1$y[,2])) ); colnames(MD2) = c("par", "val", "group");
MD3 = data.frame(md2$tsep,md2$y[,1], rep(nams[3], length(md2$y[,1])) ); colnames(MD3) = c("par", "val", "group");
MD4 = data.frame(md2$tsep,md2$y[,2], rep(nams[4], length(md2$y[,2])) ); colnames(MD4) = c("par", "val", "group");
combined_data <- rbind(MD1, MD2, MD3, MD4)
combined_data$group = factor(combined_data$group, levels = c(nams[1],nams[2], nams[3], nams[4]))
text.size = 11
pdf(file = "figures/muprhy.pdf",width = 4, height = 3) 
ggplot(combined_data, aes(x = par, y = val, color = group)) +
  geom_line(linewidth = 0.75) +
  labs(x = expression(omega), y = expression(paste("Empirical score, ", S(omega),"")), color = "Group") +
  theme_minimal() + 
  theme(text = element_text(size = text.size, family = "serif"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        panel.grid.minor = element_line(color = "lightgray"),
        legend.position = c(0.75, 0.25),  # Moves the legend closer to the bottom
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent", color = NA),  # Set transparent background for legend
        legend.key = element_rect(fill = "transparent", color = NA)
  ) + 
  guides(color = guide_legend(title = "", ncol = 1))
dev.off()
