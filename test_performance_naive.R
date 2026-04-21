


### TEST OBJECTIVES ###
# - additive covariates, nope
# - state dependent coefficients, yes
# - random effects, nope




## BASELINE :: FREE UPPER TRI + COVARIATES ##

rm(list = ls()); gc()
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states) ; track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

library(MASS); library(Rcpp); library(RcppEigen); library(TMB);
sourceCpp("FUNCS_MJP_with_eigen.cpp")
tmb_nam = "FUNCS_MJP_with_TMB"
compile(paste0(tmb_nam, ".cpp")); dyn.load(dynlib(tmb_nam))

data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = as.integer(T), use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = as.integer(F))
beta = runif(m*(m-1)/2 + length(exo.cols),0,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam, silent = TRUE)
foo <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))

pred1 =MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "free_upper_tri", link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "pade")
obs = make_Ptu_obs(m, d$s2)
rps_mjp1 = rps_vectors(m, pred1, obs)
log_mjp1 = logscore_vectors(m, pred1, obs)
mean(rps_mjp1) # === 0.4378974
mean(log_mjp1) # === 1.020314

# ============================================================================ #



## NEW IMPLEMENT :: FREE UPPER TRI + COVARIATES ##

rm(list = ls()); gc()
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states) ; track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

library(MASS); library(Rcpp); library(RcppEigen); library(TMB);
sourceCpp("FUNCS_MJP_with_eigen.cpp")
tmb_nam = "FUNCS_MJP_with_TMB"
compile(paste0(tmb_nam, ".cpp")); dyn.load(dynlib(tmb_nam))

data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = 2, use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = as.integer(F))
beta = runif(m*(m-1)/2 + length(exo.cols)*(m-1),0,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam, silent = TRUE)
foo <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))

pred1 =MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "free_upper_tri", link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "pade", state_covs = T ,warping = F)
obs = make_Ptu_obs(m, d$s2)
rps_mjp1 = rps_vectors(m, pred1, obs)
log_mjp1 = logscore_vectors(m, pred1, obs)
mean(rps_mjp1) # === 0.4378974
mean(log_mjp1) # === 1.020314

## NEW IMPLEMENT + WARPED :: FREE UPPER TRI + COVARIATES ##


rm(list = ls()); gc()
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states) ; track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

library(MASS); library(Rcpp); library(RcppEigen); library(TMB);
sourceCpp("FUNCS_MJP_with_eigen.cpp")
tmb_nam = "FUNCS_MJP_with_TMB_warped"
compile(paste0(tmb_nam, ".cpp")); dyn.load(dynlib(tmb_nam))

data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = as.integer(2), use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = as.integer(F))
beta = runif(m*(m-1)/2 + m-1 + length(exo.cols)*(m-1),-2,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam, silent = TRUE)
foo <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))

pred1 =MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "free_upper_tri", link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "pade", state_covs = T, warping = T)
obs = make_Ptu_obs(m, d$s2)
rps_mjp1 = rps_vectors(m, pred1, obs)
log_mjp1 = logscore_vectors(m, pred1, obs)
mean(rps_mjp1) # === 0.4204479
mean(log_mjp1) # === 0.9814224



#### APPENDED STATES

rm(list = ls()); gc()
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states) ; track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

library(MASS); library(Rcpp); library(RcppEigen); library(TMB);
sourceCpp("FUNCS_MJP_with_eigen.cpp")
tmb_nam = "FUNCS_MJP_with_TMB_appended"
compile(paste0(tmb_nam, ".cpp")); dyn.load(dynlib(tmb_nam))

k = 3
data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(0),cov_type = as.integer(T), use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = as.integer(F), k = k)
beta = runif(m-1 + length(exo.cols),0,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam, silent = TRUE)
foo <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))


sim <- l$simulate(par = foo$par)
pred1 <- sim$pred_mat

#pred1 =MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "gerlang", link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "pade", append = TRUE, k = k)
obs = make_Ptu_obs(m, d$s2)
rps_mjp1 = rps_vectors(m, pred1, obs)
log_mjp1 = logscore_vectors(m, pred1, obs)
mean(rps_mjp1) # ===  0.5137974
mean(log_mjp1) # === 1.536097


A_hat = make_A1(m,1:(m-1))
B <- expand_A1(m, A_hat, k)
range(rowSums(B))
range(diag(B))
min(B - diag(diag(B)))


# ============================================================================ #



## State dependent covs :: -||- ##

rm(list = ls()); gc()
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states) ; track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
#d$group <- as.integer(factor(d$Track0)) - 1; G <- length(levels(factor(d$Track0)))

library(MASS); library(Rcpp); library(RcppEigen); library(TMB);
sourceCpp("FUNCS_MJP_with_eigen.cpp")
tmb_nam = "FUNCS_MJP_with_TMB_additive"
dll <- dynlib(tmb_nam)
if (basename(dll) %in% names(getLoadedDLLs())) { dyn.unload(dll); gc() }
TMB::compile(paste0(tmb_nam, ".cpp"))
dyn.load(dll)

data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = as.integer(T), use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = as.integer(F))
beta = runif(m*(m-1)/2 + length(exo.cols)*(m-1),0,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam,checkParameterOrder=FALSE)
foo <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))
#l$env$last.par.best <- foo$par; pl <- l$env$parList(); b_hat <- as.numeric(pl$b_group)


pred1 = MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]),  generator = "free_upper_tri", link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "pade", state_covs = TRUE)
obs = make_Ptu_obs(m, d$s2)
rps_mjp2 = rps_vectors(m, pred1, obs)
log_mjp2 = logscore_vectors(m, pred1, obs)
mean(rps_mjp2) # === 0.4320385
mean(log_mjp2) # === 1.000163





## EXTENSION 1 :: MIXTURES ##


rm(list = ls()); gc()
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states) ; track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

set.seed(1); k= 5; bin_km= 0.1
d$pos_bin=floor(d$pos / bin_km); d$group_id=interaction(d$Track0, d$pos_bin, drop = TRUE); 
groups=levels(d$group_id); G=length(groups)
fold_id=sample(rep(1:k, length.out = G)); names(fold_id)=groups; fold_size=integer(k)
i = 1
test_grps=names(fold_id)[fold_id == i]
pred.idx=which(d$group_id %in% test_grps)
d.test=d[pred.idx, , drop = FALSE]; d.train=d[-pred.idx, , drop = FALSE]


library(MASS); library(Rcpp); library(RcppEigen); library(TMB);
sourceCpp("FUNCS_MJP_with_eigen.cpp")
tmb_nam = "FUNCS_MJP_with_TMB_mixture"
compile(paste0(tmb_nam, ".cpp")); dyn.load(dynlib(tmb_nam))

K <- 6
base_len <- m*(m-1)/2
mix_len  <- K-1
cov_len  <- length(exo.cols)
par_len  <- mix_len + K*base_len + cov_len
parameters <- list(theta = rep(0, par_len))  # or random init
data <- list(s1 = d.train$s1,s2 = d.train$s2,u = d.train$t,z = as.matrix(d.train[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = as.integer(T), use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = as.integer(F), K = K, ridge = 1e-4)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam, silent = TRUE)
foo <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))

foo

#sim <- l$simulate(par = foo$par)
#pred1 <- sim$pred_mat

pred1 = MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]),  generator = "free_upper_tri", link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "pade", state_covs = FALSE, mixture = TRUE, K = K)


#pred1 =MJP_predict(m = m, s1 = d$s1, u = d$t, pars = foo$par , z = as.matrix(d[,exo.cols]), generator = "gerlang", link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "pade", append = TRUE, k = k)
obs = make_Ptu_obs(m, d$s2)
rps_mjp1 = rps_vectors(m, pred1, obs)
log_mjp1 = logscore_vectors(m, pred1, obs)
mean(rps_mjp1) # ===  
mean(log_mjp1) # === 
























