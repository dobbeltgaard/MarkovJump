
###################################
### MAKE PREDICTIONS FOR MARCEL ###
###################################
rm(list = ls()); gc()
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states) ; track <- unique(d$Track); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

set.seed(1); k= 5; bin_km= 0.1
d$pos_bin=floor(d$pos / bin_km); d$group_id=interaction(d$Track0, d$pos_bin, drop = TRUE); 
groups=levels(d$group_id); G=length(groups)
fold_id=sample(rep(1:k, length.out = G)); names(fold_id)=groups; fold_size=integer(k)
i = 1; test_grps=names(fold_id)[fold_id == i]
pred.idx=which(d$group_id %in% test_grps)
d.test=d[pred.idx, , drop = FALSE]; d.train=d[-pred.idx, , drop = FALSE]

library(MASS); library(Rcpp); library(RcppEigen); library(TMB);
sourceCpp("FUNCS_MJP_with_eigen.cpp")
tmb_nam = "FUNCS_MJP_with_TMB"
compile(paste0(tmb_nam, ".cpp")); dyn.load(dynlib(tmb_nam))

K = 1
data <- list(s1 = d.train$s1,s2 = d.train$s2,u = d.train$t,z = as.matrix(d.train[, exo.cols]),m = m,generator_type = as.integer(2),cov_type = 1, use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = as.integer(F) )#, K = K, ridge = 10e-4)
beta = runif(K-1 + m*(m-1)/2 * K + length(exo.cols),-1,0); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam, silent = TRUE)
foo <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))


year = 2
time = seq(0,year,year/(52*year))
PREDS <- vector("list", length(time))
names(PREDS) <- as.character(time)

for (k in seq_along(time)) {
  j <- time[k]
  u <- rep(j, NROW(d.test))
  pred1 = MJP_predict(m = m, s1 = d.test$s1, u = u, pars = foo$par , z = as.matrix(d.test[,exo.cols]), generator = "free_upper_tri", link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "pade", state_covs = F, warping = F, mixture = FALSE)
  PREDS[[k]] <- pred1
}


#MJP_predict(m = m, s1 = c(1), u = (0.01), pars = foo$par , z = as.matrix(d.test[,exo.cols]), generator = "free_upper_tri", link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "pade", state_covs = F, warping = F, mixture = TRUE, K = K)

PREDS_df <- do.call(rbind, lapply(seq_along(PREDS), function(k) {
  mat <- PREDS[[k]]
  df  <- as.data.frame(mat)
  names(df) <- paste0("V", seq_len(ncol(mat)))
  df$time   <- time[k]
  df$ID     <- seq_len(nrow(mat))
  df$Track0 <- d.test$Track0
  df[, c("time", "ID", "Track0", paste0("V", seq_len(ncol(mat))))]
}))

colnames(PREDS_df)[2:8] <- c("defect_ID","Line","class_3", "class_2B", "class_2A", "class_1", "class_0")
write.csv(PREDS_df, "predictions_decision_making.csv", row.names = FALSE)


## Make TPMs for Marcel
library(jsonlite)
set.seed(1997); D <- NROW(d.test); 
times <- c(0, 0.5, 1, 2, 3.5)
dtimes <- diff(times); H <- length(times)
idx <- sample(1:nrow(d.test), size = D, replace = FALSE); dat <- d.test[idx, ]
TPM <- vector("list", D)
for (s in seq_len(D)) {
  defect_list <- list()
  defect_list$defect_id <- paste0("defect_", s)
  #defect_list$times <- as.list(times)
  initial = rep(0.0,m); initial[dat$s1[s]] = 1.0
  defect_list$initial = initial
  defect_list$TPM <- list()
  for (tt in dtimes) {
    tpm <- matrix(NA_real_, m, m)
    for (i in seq_len(m)) {tpm[i,] = MJP_predict(m = m,s1 =i,u=tt ,pars = foo$par, z = as.matrix(dat[s, exo.cols, drop = FALSE]), generator = "free_upper_tri", link_type_base = "exp",link_type_covs = "exp",covs_bin = TRUE,transient_dist_method = "pade",state_covs = FALSE, warping = FALSE,mixture = FALSE)}
    defect_list$TPM[[paste0("dt_", tt)]] = tpm
  }
  TPM[[s]] <- defect_list
}

write_json(TPM, "tpm_list.json", pretty = TRUE, auto_unbox = TRUE)


TPMS <- jsonlite::fromJSON("tpm_list.json", simplifyDataFrame = FALSE)
TPMS[[1]]$TPM$dt_0.5

entropy <- function(x, eps = 1e-12, normalize = FALSE) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  s <- sum(x)
  if (s <= 0) return(NA_real_)
  p <- x / s
  p <- pmax(p, eps)
  H <- -sum(p * log(p))
  if (normalize) H <- H / log(length(p))
  H
}

entropy(t(c(1,0,0,0,0)) %*% TPMS[[3]]$TPM$dt_0.5)
entropy(t(c(1,0,0,0,0)) %*% TPMS[[5]]$TPM$dt_0.5)
entropy(t(c(1,0,0,0,0)) %*% TPMS[[6]]$TPM$dt_0.5)


### Monte Carlo paths
set.seed(1997)
D = 10; S = 1000
times = c(0,0.5,1,2); H = length(times)
idx = sample(1:NROW(d.test), size = D, replace = FALSE)
dat = d.test[idx, ]

paths <- array(NA_integer_, dim = c(S, D, H),dimnames = list(scenario = paste0("s", 1:S),defect= paste0("d", 1:D),time  = as.character(times)))
paths[, , 1] <- matrix(rep(dat$s1, each = S), nrow = S, ncol = D, byrow = FALSE)
draw_cat <- function(p) sample.int(length(p), size = 1, prob = p)

for (s in seq_len(S)) {
  cur_states <- as.integer(paths[s, , 1])
  for (h in 2:H) {
    dt <- times[h] - times[h - 1]; u=rep(dt, D)
    P <- MJP_predict(m = m,s1 = cur_states,u  = u,pars = foo$par, z = as.matrix(dat[, exo.cols, drop = FALSE]), generator = "free_upper_tri", link_type_base = "exp",link_type_covs = "exp",covs_bin = TRUE,transient_dist_method = "pade",state_covs = FALSE, warping = FALSE,mixture = FALSE)
    nxt_states <- integer(D)
    for (d in seq_len(D)) nxt_states[d] <- draw_cat(P[d, ])
    paths[s, , h] <- nxt_states
    cur_states <- nxt_states
  }
}

paths_df <- do.call(rbind, lapply(seq_len(S), function(s) {
  do.call(rbind, lapply(seq_len(D), function(d) { data.frame(scenario = s,defect= d,time= times,state = as.integer(paths[s, d, ]),Track0   = dat$Track0[d],stringsAsFactors = FALSE)}))
}))

library(jsonlite)

full_tree <- list(
  meta = list(m = m,S = S,D = D,H = H),times = times,
  defects = data.frame(defect = seq_len(D), row_id = idx,Track0 = dat$Track0,s1 = dat$s1,stringsAsFactors = FALSE),
  paths = lapply(seq_len(S), function(s) {list(scenario = s,states = lapply(seq_len(H), function(h) as.integer(paths[s, , h])))})
)

writeLines(toJSON(full_tree, auto_unbox = TRUE, digits = NA, pretty = FALSE), "scenario_full.json")





#############################################
### SIMULATED DATA AND ESTIMATION TESTING ###
#############################################
rm(list = ls()); gc()
library(MASS); library(Rcpp); library(RcppEigen); library(TMB);
sourceCpp("FUNCS_MJP_with_eigen.cpp")
tmb_nam = "FUNCS_MJP_with_TMB_warped"
compile(paste0(tmb_nam, ".cpp")); dyn.load(dynlib(tmb_nam))
set.seed(1)
d = read.csv("defect_data.csv"); states <- c(1,2,3,4,5); m <- length(states); exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

gen = "bidiagonal"; warping = "warp"; const = 5
if(gen == "gerlang"){generator_type = 0; beta_base = (m-1):1/const; }
if(gen == "gerlang_relax"){generator_type = 1; beta_base = (m-1):1/const; }
if(gen == "free_upper_tri"){generator_type = 2; beta_base = (m*(m-1)/2):1/const; }
if(gen == "bidiagonal"){generator_type = 3; beta_base = (2*m-3):1/const; }
if(gen == "tridiagonal"){generator_type = 4; beta_base = (3*m-6):1/const; }
if(warping == "no_warp"){ warp = F; xi = c(); }
if(warping == "warp"){warp = T;  xi = runif(m-1, min = 0, max = 3); }

K = 1
data <- list(s1 = d$s1,s2 = d$s2,u = d$t,z = as.matrix(d[, exo.cols]),m = m,generator_type = generator_type,cov_type = 1, use_log_score = as.integer(T),use_rps_score = as.integer(T), use_brier_score = as.integer(F) )#, K = K, ridge = 10e-4)
beta = runif(K-1 + length(beta_base) * K + length(xi) + length(exo.cols),-1,0); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam, silent = TRUE)
foo <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))


probs = table(d$s1)/sum(table(d$s1)) #m = 5; prob = 1; probs = c() #for(i in 1:m){prob = prob/2; probs = c(probs, prob)}
nsim=10000#NROW(d)
s1 = sample(1:m, replace = TRUE, size = nsim, prob = probs)
u <- sample(d$t, size = nsim, replace = TRUE) #u = runif(nsim, min = 0, max = 2)
zsam = apply(X = as.matrix(d[, exo.cols]), MARGIN = 2, FUN = sample, replace = TRUE, size = nsim)

par = foo$par#c(-beta_base,-xi)

pred = MJP_predict(m = m, s1 = s1, u = u, pars = par , z = zsam, generator = gen, link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "pade", state_covs = F, warping = warp)
s2 <- apply(pred, 1, function(p) {sample.int(m, size = 1, prob = p)})

data <- list(s1 = s1,s2 = s2,u = u,z = zsam,m = m,generator_type = as.integer(generator_type),cov_type = 1, use_log_score = as.integer(T),use_rps_score = as.integer(F), use_brier_score = as.integer(F))
beta = runif(length(par),0,1); parameters <- list(theta = beta)
l <- MakeADFun(data = data, parameters = parameters, DLL = tmb_nam, silent = TRUE)
foo_recover <- nlminb(l$par, l$fn, l$gr, control = list(eval.max = 2000, iter.max = 2000))

write.csv(x = data.frame(y1=s1,y2=s2,tau=u,z1=zsam[,1],z2=zsam[,2],z3=zsam[,3],z4=zsam[,4],z5=zsam[,5]),file = "synthetic_data.csv", row.names = F)

#COMPARE ESTIMATES WITH KNOWN VALUES
lambda_len <-length(beta_base)
lambda_hat <- exp(foo_recover$par[1:lambda_len])
lambda_true <- exp(foo$par[1:lambda_len])
cbind(lambda_true = lambda_true, lambda_hat = lambda_hat, diff = lambda_hat-lambda_true)

idx <- (lambda_len + 1):(lambda_len + m - 1)
normalize_gm <- function(x) x / exp(mean(log(x)))
xi_hat <- normalize_gm(exp(foo_recover$par[idx]))
xi_true <- normalize_gm(exp(par[idx]))
cbind(xi_true = xi_true, xi_hat = xi_hat, diff = xi_hat-xi_true)

beta_len <- length(exo.cols)
idx_beta <- (lambda_len + (m - 1) + 1):(lambda_len + (m - 1) + beta_len)
beta_hat  <- foo_recover$par[idx_beta]
beta_true <- par[idx_beta]
cbind(beta_true = beta_true,beta_hat  = beta_hat,diff = beta_hat - beta_true)

sum(log(xi_hat))   
sum(log(xi_true))
mean(log(xi_hat)); mean(log(xi_true))


