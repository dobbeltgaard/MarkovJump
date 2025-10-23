
rm(list = ls()); gc(); 

d = read.csv("defect_data.csv")
states <- c(1,2,3,4,5)
m <- length(states) 
track <- unique(d$Track)
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

library(MASS); library(Rcpp); library(RcppEigen); library(TMB); #library(trust)
sourceCpp("FUNCS_MJP_with_eigen.cpp")

reload_tmb <- function(src, flags = "-O0 -g", clean = TRUE) {
  stopifnot(file.exists(src))
  base <- sub("[.]cpp$", "", normalizePath(src, winslash = "/"))
  so   <- TMB::dynlib(base)             # absolute path to the .so/.dll
  dlls <- getLoadedDLLs()
  if (basename(so) %in% names(dlls) || so %in% vapply(dlls, `[[`, "", "path")) {
    try(dyn.unload(so), silent = TRUE)
  }
  TMB::compile(src, flags = flags, clean = clean)
  stopifnot(file.exists(so))
  dyn.load(so)
  invisible(so)
}

d.train = d
d.test = d
obs = make_Ptu_obs(m, d.test$s2)



gen = "bidiagonal"; 
if(gen == "gerlang"){generator_type = 0; beta_base = rep(-1,m-1); }
if(gen == "gerlang_relax"){generator_type = 1; beta_base = rep(-1,m-1); }
if(gen == "free_upper_tri"){generator_type = 2; beta_base = rep(-1,m*(m-1)/2); }
if(gen == "bidiagonal"){generator_type = 3; beta_base = rep(-1,2*m-3); }
if(gen == "tridiagonal"){generator_type = 4; beta_base = rep(-1,3*m-6); }
xi = rep(-1,m-1); 

tmb_nam = "FUNCS_MJP_with_TMB_warped"
reload_tmb(paste0(tmb_nam, ".cpp"))
#compile(paste0(tmb_nam, ".cpp")); dyn.load(dynlib(tmb_nam))
data <- list(s1 = d.train$s1,s2 = d.train$s2,u = d.train$t,z = as.matrix(d.train[, exo.cols]),m = m,generator_type = generator_type,cov_type = as.integer(T),use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = 0)
beta = c(beta_base, xi, rep(0,length(exo.cols)))
parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam)
foo1 = NULL; #foo1 <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))
foo1 <- nlminb(l$par, l$fn, l$gr, l$he)
pred = MJP_predict(m = m, s1 = d.test$s1, u = d.test$t, pars = foo1$par, z = as.matrix(d.test[,exo.cols]), generator = gen, link_type_base = "exp", link_type_covs ="exp", covs_bin = T, transient_dist_method = "pade", warping = T)
err1 = rps_vectors(m, pred, obs)
mean(err1)


tmb_nam = "FUNCS_MJP_with_TMB_warped_v2"
reload_tmb(paste0(tmb_nam, ".cpp"))
#compile(paste0(tmb_nam, ".cpp"));dyn.load(dynlib(tmb_nam))
data <- list(s1 = d.train$s1,s2 = d.train$s2,u = d.train$t,z = as.matrix(d.train[, exo.cols]),m = m,generator_type = generator_type,cov_type = as.integer(T),use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = 0,
             lambda_reg=0)
beta = foo1$par #c(beta_base, xi, rep(0,length(exo.cols)))
parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam)
foo2 = NULL; foo2 <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))
foo2





################################
### SIMULATION STUDY TESTING ###
################################
library(tibble)
library(dplyr)
library(purrr)
library(Matrix)  # for expm()




simulate_MCM = function(m, A, N,Tt){
  s1 = sample(1:m, N, replace = TRUE, prob = 1/(1:m)/sum(1/(1:m)))
  s2 = rep(NA, N)
  tau = rlnorm(N, meanlog = Tt, 1)
  for(i in 1:N){
    e <- rep(0, m); e[s1[i]] <- 1
    v2 = as.vector(t(e) %*% as.matrix(expm(A*tau[i])))
    s2[i] = sample(1:m, 1, prob = v2, replace = T)
  }
  return(tibble(s1 = s1, s2 = s2, t = tau))
}


m <- 5
A <- make_A3(m, rep(0.5, m^3))
set.seed(1)

d = as.data.frame(simulate_MCM(m, A, 3600, 0.5))
split = floor(NROW(d)*0.7)
d.train = d[1:split, ]
d.test = d[(split+1):NROW(d),]




library(MASS); library(Rcpp); library(RcppEigen); library(TMB); #library(trust)
sourceCpp("FUNCS_MJP_with_eigen.cpp")

gen = "free_upper_tri"; 
if(gen == "gerlang"){generator_type = 0; beta_base = rep(-1,m-1); }
if(gen == "gerlang_relax"){generator_type = 1; beta_base = rep(-1,m-1); }
if(gen == "free_upper_tri"){generator_type = 2; beta_base = rep(-1,m*(m-1)/2); }
if(gen == "bidiagonal"){generator_type = 3; beta_base = rep(-1,2*m-3); }
if(gen == "tridiagonal"){generator_type = 4; beta_base = rep(-1,3*m-6); }
xi = rep(-1,m-1); 

tmb_nam = "FUNCS_MJP_with_TMB"
reload_tmb(paste0(tmb_nam, ".cpp"))
data <- list(s1 = d.train$s1,s2 = d.train$s2,u = d.train$t,z = as.matrix(1), m = m, generator_type = generator_type,cov_type = as.integer(F),use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = 0)
beta = c(beta_base)
parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam)
foo1 = NULL; foo1 <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))
exp(foo1$par)

obs = make_Ptu_obs(m, d.test$s2)
pred = MJP_predict(m = m, s1 = d.test$s1, u = d.test$t, pars = foo1$par, z = as.matrix(1), generator = gen, link_type_base = "exp", link_type_covs ="exp", covs_bin = F, transient_dist_method = "pade", warping = F)
err1 = rps_vectors(m, pred, obs)
err12 = logscore_vectors(m, pred, obs)
mean(err1)
-mean(err12)

tmb_nam = "FUNCS_MJP_with_TMB_warped"
reload_tmb(paste0(tmb_nam, ".cpp"))
data <- list(s1 = d.train$s1,s2 = d.train$s2,u = d.train$t,z = as.matrix(1), m = m, generator_type = generator_type,cov_type = as.integer(F),use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = 0)
beta = c(beta_base, xi)
parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam)
foo2 = NULL; foo2 <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))
exp(foo2$par)

obs = make_Ptu_obs(m, d.test$s2)
pred = MJP_predict(m = m, s1 = d.test$s1, u = d.test$t, pars = foo2$par, z = as.matrix(1), generator = gen, link_type_base = "exp", link_type_covs ="exp", covs_bin = F, transient_dist_method = "pade", warping = T)
err2 = rps_vectors(m, pred, obs)
err22 = logscore_vectors(m, pred, obs)
mean(err2)
-mean(err22)
l$fn(foo2$par)


#GERLANG simulation, GERLANG fit (correct model specification)
#without warping: RPS = 0.3294563; logS = -0.8505563;
#with warping: RPS = 0.3304099; logS = -0.8533674; 

#GERLANG simulation, TRIDIAGONAL fit (overcomplex model specification)
#without warping: RPS = 0.3295472; LogS = -0.8514408
#with warping: RPS =  0.3306621; LogS = -0.8551121; 

#TRIDIAGONAL simulation, GERLANG fit (oversimple model specification)
#without warping: RPS = 0.4061664; LogS = -0.9568657;
#with warping: RPS = 0.3896755; LogS = -0.8943784; 



