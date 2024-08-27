



rm(list = ls()) #clear memory
d = read.csv("defect_data.csv")
states <- c(1,2,3,4,5) #define states
m <- length(states) #number of states
track <- unique(d$Track) #define investigated tracks
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
z = as.matrix(d[, exo.cols])

library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")
library("bench")
source("speed_testing/score_functions_implementations_R.R")



beta0 <- c(c(0.1,0.2,0.3,0.4),
           #rep(0.25, (m-1)),
           #c(0.1,0.2,0.3,0.4, 0.2,0.2,0.3,0.4,0.1,9),
           rep(0.1,length(exo.cols)))  





time_implementation_transient_dist <- 
  bench::mark(
  log.lik.cov(m = m, s1 = d$s1, s2 = d$s2, u = d$t, z = z, beta = beta0, states = states),
  rps.score.cov(m = m, s1 = d$s1, s2 = d$s2, u = d$t, z = z, beta = beta0, states = states),
  brier.score.cov(m = m, s1 = d$s1, s2 = d$s2, u = d$t, z = z, beta = beta0, states = states),
  
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization"),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = F, rps_bin = T, brier_bin = F, transient_dist_method = "uniformization"),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = F, rps_bin = F, brier_bin = T, transient_dist_method = "uniformization"),
  
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "pade"),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = F, rps_bin = T, brier_bin = F, transient_dist_method = "pade"),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = F, rps_bin = F, brier_bin = T, transient_dist_method = "pade"),
  
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "eigenvalue_decomp"),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = F, rps_bin = T, brier_bin = F, transient_dist_method = "eigenvalue_decomp"),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = F, rps_bin = F, brier_bin = T, transient_dist_method = "eigenvalue_decomp"),
  
  , relative = F, check = F)

time_implementation_transient_dist
library(kableExtra)
options(knitr.table.format = "latex")
mat = t(matrix(time_implementation_transient_dist$median, nrow = 3))
rownames(mat) = c("r_pade", "cpp_uni", "cpp_pade", "cpp_eigen")
colnames(mat) = c("logS", "RPS", "Brier")

kable(
  mat,
  booktabs = T, digits = 2, row.names = T
)




time_uniformization_eps <- 
  bench::mark(
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-52)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-48)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-44)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-40)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-36)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-32)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-28)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-24)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-20)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-16)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-12)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-8)),
  MJP_score(m = m, s1 = d$s1, s2 = d$s2, u = d$t,pars = beta0, z = z, generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F, transient_dist_method = "uniformization", eps = 2^(-4)),
  , relative = T, check = F)




library(kableExtra)



options(knitr.table.format = "latex")
kable(t(matrix(c(rep(NA, length(time_uniformization_eps$median)),time_uniformization_eps$median), ncol = 2)), 
      booktabs = T, digits = 2, row.names = T)


domain = seq(52,4,by=-4)
missing_prop_mass = rep(0, length(domain))
count = 0
for(i in domain){count = count + 1;  missing_prop_mass[count] = 2^(-i)}

plot(missing_prop_mass, time_uniformization_eps$median)
