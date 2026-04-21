
rm(list = ls()); gc()

d = read.csv("defect_data.csv")
states <- c(1,2,3,4,5)
m <- length(states) 
track <- unique(d$Track)
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")


library(MASS); library(Rcpp); library(RcppEigen); library(TMB);
sourceCpp("FUNCS_MJP_with_eigen.cpp")
tmb_nam = "FUNCS_MJP_with_TMB_warped"
compile(paste0(tmb_nam, ".cpp"))
dyn.load(dynlib(tmb_nam))
warpings = c("warp") 

### INITIALIZE k-fold ###
set.seed(1); k= 5; bin_km= 0.1
d$pos_bin=floor(d$pos / bin_km); d$group_id=interaction(d$Track0, d$pos_bin, drop = TRUE); 
groups=levels(d$group_id); G=length(groups)
fold_id=sample(rep(1:k, length.out = G)); names(fold_id)=groups; fold_size=integer(k)

PREDS = list(); PARS = list(); CONV <- list(); GRADNORM <- list()
log_bin = F; rps_bin = F; brier_bin = F;
parameterizations = c("gerlang", "gerlang_relax","free_upper_tri","bidiagonal", "tridiagonal")
scores = c("all")
links = c("exp")
covsforms = c(0,1,2)


if (!dir.exists("estimates_V3")) dir.create("estimates_V3")
if (!dir.exists("predictions_V3")) dir.create("predictions_V3")
if (!dir.exists("diagnostics_V3")) dir.create("diagnostics_V3")

for(i in 1:5){
  start_time <- Sys.time()
  
  #Estimation and test set splits
  test_grps=names(fold_id)[fold_id == i]
  pred.idx=which(d$group_id %in% test_grps)
  d.test=d[pred.idx, , drop = FALSE]; d.train=d[-pred.idx, , drop = FALSE]
  fold_size[i] <- nrow(d.test)
  
  obs = make_Ptu_obs(m, d.test$s2)
  write.csv(obs, file = paste0("predictions_V3/", "obs", "_fold_", i, ".csv"), row.names = FALSE)
  
  ####################################################
  ### Estimate and predict with all MJP variations ###
  ####################################################
  count = 0
  for(gen in parameterizations){
    if(gen == "gerlang"){generator_type = 0; beta_base = rep(-1,m-1); }
    if(gen == "gerlang_relax"){generator_type = 1; beta_base = rep(-1,m-1); }
    if(gen == "free_upper_tri"){generator_type = 2; beta_base = rep(-1,m*(m-1)/2); }
    if(gen == "bidiagonal"){generator_type = 3; beta_base = rep(-1,2*m-3); }
    if(gen == "tridiagonal"){generator_type = 4; beta_base = rep(-1,3*m-6); }
    for(cov in covsforms){
      for(baselink in links){
        for(covslink in links){
          for(warping in warpings){
            if(warping == "no_warp"){ warp = F; xi = c(); }
            if(warping == "warp"){warp = T;  xi = rep(-1, m-1); }
            for(score in scores){
              if(score == "log" | score == "all"){log_bin = T}
              if(score == "rps" | score == "all"){rps_bin = T}
              #if(score == "brier" | score == "all"){brier_bin = T}
              count = count + 1
              
              nam = paste(gen,paste0("cov",cov),baselink,covslink,score,warping,sep = "_") #!
              print(nam)
              
              if(cov == 3 ){ beta = c(beta_base, xi, rep(0,length(exo.cols) * (m-1))); state_covs = T; } #!
              if(cov == 2 ){ beta = c(beta_base, xi, rep(0,length(exo.cols) * (m-1))); state_covs = T; } #!
              if(cov == 1 ){ beta = c(beta_base, xi, rep(0,length(exo.cols))); state_covs = F; } #!
              if(cov == 0){ beta = c(beta_base,xi); state_covs = F; } #!
              
              data.train <- list(s1 = d.train$s1,s2 = d.train$s2,u = d.train$t,z = as.matrix(d.train[, exo.cols]),m = m,generator_type = generator_type,cov_type = as.integer(cov),use_log_score = as.integer(log_bin),use_rps_score = as.integer(rps_bin), use_brier_score = 0)
              parameters <- list(theta = beta)
              
              #Estimate and store model pars
              #foo = NULL; foo = optim( par = beta, fn = MJP_score, m = m, s1 = d.train$s1, s2 = d.train$s2, u = d.train$t, z = as.matrix(d.train[,exo.cols]), generator = gen, link_type_base = baselink, link_type_covs = covslink, covs_bin = cov, likelihood_bin = log_bin, rps_bin = rps_bin, brier_bin = brier_bin, transient_dist_method = "eigenvalue_decomp", warping = warp, method = "BFGS", control = list(maxit = 1000)) #estimate model 
              l <- MakeADFun(data = data.train, parameters = parameters, DLL = tmb_nam, silent = TRUE)
              foo = NULL; foo <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))
              grad_opt = l$gr(foo$par); grad_norm = sqrt(sum(grad_opt^2)) / length(data.train$s1)
              if(foo$convergence != 0) print("not converged")
              rm(l); invisible(gc())
              
              #compute AIC
              data.train.AIC <- list(s1 = d.train$s1,s2 = d.train$s2,u = d.train$t,z = as.matrix(d.train[, exo.cols]),m = m,generator_type = generator_type,cov_type = as.integer(cov),use_log_score = as.integer(1),use_rps_score = as.integer(0), use_brier_score = 0)
              parameters.AIC <- list(theta = foo$par)
              lAIC <- MakeADFun(data = data.train.AIC, parameters = parameters.AIC, DLL = tmb_nam,silent = TRUE)
              nll_hat <- lAIC$fn(foo$par) 
              AIC_hat <- 2 * length(foo$par) + 2 * nll_hat
              
              write.csv(data.frame(convergence = foo$convergence,grad_norm = grad_norm, AIC = AIC_hat),
                        file = paste0("diagnostics_V3/", nam, "_fold_", i, ".csv"), row.names = FALSE)
              
              write.csv(data.frame(par = foo$par), file = paste0("estimates_V3/", nam, "_fold_", i, ".csv"), row.names = FALSE)
              if(i == 1){ PARS[[nam]] = foo$par} else { PARS[[nam]] = rbind(PARS[[nam]],foo$par)} #store model estimates
              #if(i == 1){PARS[[nam]] = foo$par; CONV[[nam]] = foo$convergence; GRADNORM[[nam]] = grad_norm} else {PARS[[nam]] = rbind(PARS[[nam]], foo$par); CONV[[nam]] = c(CONV[[nam]], foo$convergence); GRADNORM[[nam]] = c(GRADNORM[[nam]], grad_norm); }
              
              #Make and store predictions and compute error scores
              pred = MJP_predict(m = m, s1 = d.test$s1, u = d.test$t, pars = foo$par, z = as.matrix(d.test[,exo.cols]), generator = gen, link_type_base = baselink, link_type_covs = covslink, covs_bin = isTRUE(cov>0), transient_dist_method = "pade", warping = warp, state_covs = state_covs) #!
              write.csv(pred, file = paste0("predictions_V3/", nam, "_fold_", i, ".csv"), row.names = FALSE)
              if(i == 1){ PREDS[[nam]] = pred} else { PREDS[[nam]] = rbind(PREDS[[nam]],pred)} #store predictions
              
              log_bin = F; rps_bin = F; brier_bin = F; #reset score bools
            }
          }
        }
      }
    }
  }
  end_time <- Sys.time(); runtime <- end_time - start_time
  cat("Fold ", i, "runtime:", runtime, "\n") #approx 1.65 mins per fold
}
print(fold_size) #[1] 703 725 751 752 727


