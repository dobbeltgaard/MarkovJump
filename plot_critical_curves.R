rm(list = ls()); gc() #clear memory
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states) ; track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
library(MASS); library(Rcpp); library(RcppEigen); library(TMB);
sourceCpp("FUNCS_MJP_with_eigen.cpp")
tmb_nam = "FUNCS_MJP_with_TMB"
compile(paste0(tmb_nam, ".cpp"))
dyn.load(dynlib(tmb_nam))


gen = "bidiagonal"; warping = "no_warp"; cov = 2; log_bin = T; rps_bin = T;
if(gen == "gerlang"){generator_type = 0; beta_base = rep(-1,m-1); }
if(gen == "gerlang_relax"){generator_type = 1; beta_base = rep(-1,m-1); }
if(gen == "free_upper_tri"){generator_type = 2; beta_base = rep(-1,m*(m-1)/2); }
if(gen == "bidiagonal"){generator_type = 3; beta_base = rep(-1,2*m-3); }
if(gen == "tridiagonal"){generator_type = 4; beta_base = rep(-1,3*m-6); }
if(warping == "no_warp"){ warp = F; xi = c(); }
if(warping == "warp"){warp = T;  xi = rep(-1, m-1); }

if(cov == 2 ){ beta = c(beta_base, xi, rep(0,length(exo.cols) * (m-1))); state_covs = T; } #!
if(cov == 1 ){ beta = c(beta_base, xi, rep(0,length(exo.cols))); state_covs = F; } #!
if(cov == 0){ beta = c(beta_base,xi); state_covs = F; } #!

d.train = d;
data.train <- list(s1 = d.train$s1,s2 = d.train$s2,u = d.train$t,z = as.matrix(d.train[, exo.cols]),m = m,generator_type = generator_type,cov_type = as.integer(cov),use_log_score = as.integer(log_bin),use_rps_score = as.integer(rps_bin), use_brier_score = 0)
parameters <- list(theta = beta)

l <- MakeADFun(data = data.train, parameters = parameters, DLL = tmb_nam, silent = TRUE)
foo = NULL; foo <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))

for(ii in 1:5){
  qt = apply(d.train[, exo.cols], 2, quantile, probs = c(0.1, 0.5, 0.9))
  means = apply(d.train[, exo.cols], 2, mean)

  col = ii
  exos = matrix(NA, 3, length(exo.cols))
  for(i in 1:length(exo.cols)){
    for(j in 1:3){
      if(i == col & j != 2){exos[j,i] = qt[j,i]} 
      else { exos[j,i] = means[i]}#qt[2,i]}
    }
  }
  
  year = 5; time = seq(0,year,year/(365*year))
  PREDS <- vector("list", length(time)); names(PREDS) <- as.character(time)
  
  init =rep(1,3)
  for (k in seq_along(time)) {
    j <- time[k]; u <- rep(j, length(init))
    pred1 = MJP_predict(m = m, s1 = init, u = u, pars = foo$par, z = as.matrix(exos), generator = gen, link_type_base = "exp", link_type_covs = "exp", covs_bin = isTRUE(cov>0), transient_dist_method = "pade", warping = warp, state_covs = state_covs) #!
    PREDS[[k]] <- pred1
  }
  
  risk_curves <- t(sapply(PREDS, function(P) rowSums(P[, (m-1):m, drop = FALSE])))
  risk <- data.frame(t = time,low = risk_curves[, 1],mid  = risk_curves[, 2],high  = risk_curves[, 3])
  
  base_size = 11; k = 1; rel_number = 1
  
  df <- risk
  library(ggplot2)
  library(geomtextpath)
  
  lab_mid  <- transform(df, y = mid,  label = "Median")
  lab_low  <- transform(df, y = low,  label = "Low (10%)")
  lab_high <- transform(df, y = high, label = "High (90%)")
  
  p = ggplot(df, aes(x = t)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
    #geom_line(aes(y = mid), linewidth = 0.8) +
    geom_textline(data = lab_mid,  aes(x = t, y = mid),linewidth = 0.7, hjust = 0.5, size = 4, fontface = 1, label = deparse(bquote(Q[50] ~ '  ')), parse = T, family = "serif") +
    geom_textline(data = lab_mid,  aes(x = t, y = low),linewidth = 0.1, hjust = 0.95, size = 4, fontface = 1, label = deparse(bquote(Q[10] ~ '  ')), parse = T, family = "serif") +
    geom_textline(data = lab_mid,  aes(x = t, y = high),linewidth = 0.1, hjust = 0.95, size = 4, fontface = 1, label = deparse(bquote(Q[90] ~ '  ')), parse = T, family = "serif") +
    theme(
      text             = element_text(size = base_size * k, family = "serif"),
      axis.text.x      = element_text(size = rel(rel_number), vjust = 0.3),
      axis.text.y      = element_text(size = rel(rel_number)),
      legend.title     = element_blank(),
      legend.text      = element_text(size = rel(0.75*rel_number)),
      legend.key.width = unit(base_size/11*0.75, "lines"),
      legend.key.height= unit(base_size/11*0.6, "lines"),
      panel.background = element_rect(fill = "white", color = "grey"),
      panel.grid.minor = element_line(color = "lightgray"),
      legend.position  = c(0.5, 0.9),
      legend.direction = "horizontal",
      legend.background= element_rect(fill = "transparent", color = NA),
      legend.key       = element_rect(fill = "transparent", color = NA)
    ) +
    guides(color = guide_legend(nrow = 1)) +
    xlab("Time (years)") + ylab("Probability")
  
  #pdf(file = nfil,width = 4, height = 3) 
  #ggsave(paste0("figures/crit_curves_",col,".pdf"), p, width = 15, height = 11, units = "in")
  
  nfil = paste0("figures/crit_curves_",col,".pdf")
  pdf(file = nfil,width = 3, height = 3) 
  print(p)
  dev.off()
  graphics.off()
}
