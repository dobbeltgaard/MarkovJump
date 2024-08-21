
rm(list = ls()) #clear memory

#setwd("C:/Users/atbd/OneDrive - COWI/Desktop/Bane NOR transmission rates estimation") #set proper working directory
#setwd("C:/Users/askbi/OneDrive - COWI/Desktop/Bane NOR transmission rates estimation")
#source("funcs.R") #define functions from another R file

library(RcppEigen); library(Rcpp)

source("funcs_diagonalization.R")
source("funcs_discrete_loglik.R")
source("funcs_forecasting.R")
source("funcs_helping.R")
source("funcs_mcem.R")


load("defects_covs_base.RData"); D <- foo; rm(foo) #load dat

#########################
### Data preparation ####
#########################
states <- c(1,2,3,4,5) #define states
#states <- c(2,3,4,5)
m <- length(states) #number of states
track <- unique(D$Track) #define investigated tracks
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
#idx <- D$Track0 %in% track #choose all tracks 
#idx <- Dtot$Track0 == "nord" #select single tracks

idx.remove <- 
  (D$`s-` == 4 & D$`s` == 4 & D$u > 300) | 
  (D$`s-` == 4 & D$`s` == 5 & D$u > 300) | 
  (D$`s-` == 5 & D$`s` == 5 & D$u > 300) 


d <- D[!idx.remove, c("pos", "s-", "s", "u", exo.cols, "Track0")] #define data set of interest



#########################
### DIRECT MAXIMATION ###
#########################

beta0 <- c(rep(0.25, (m-1)), rep(0.1,length(exo.cols)))  
res3 <- mle.cov(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, 
                z = as.matrix(d[,exo.cols]), beta = beta0, 
                pvalues.bin =  T, method = "BFGS")


res3
base.foo <- 1/exp(res3$lambda[1:(m-1)])*365
exo.foo <- apply(sweep((d[,exo.cols]), 2, res3$lambda[m:(m-1+length(exo.cols))], "*"), 
                 MARGIN = 2, quantile, probs = 0.5)

base.foo
base.foo/exp(sum(exo.foo))
sum(base.foo/exp(sum(exo.foo)))/365



#RPS estimation
source("funcs_rps.R")

exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")
beta0 <- c(rep(0.25, (m-1)), rep(0.1,length(exo.cols)))  

obs = make_Ptu_obs(m, d.test$`s`)
res4 <- rps.estim.cov(m = m, s1 = d$`s-`, s2 = d.train$s, u = d.train$u/365, z = as.matrix(d.train[,exo.cols]), beta = beta0, method = "BFGS")
res4pred = jump_prediction_cov(m, s1 = d.test$`s-`, u = d.test$u/365, z = as.matrix(d.test[,exo.cols]), res4$lambda)
mean(rps_vectors(m, res4pred, obs))



res4
base.foo <- 1/exp(res4$lambda[1:(m-1)])*365
exo.foo <- apply(sweep((d[,exo.cols]), 2, res3$lambda[m:(m-1+length(exo.cols))], "*"), 
                 MARGIN = 2, quantile, probs = 0.5)

base.foo
base.foo/exp(sum(exo.foo))
sum(base.foo/exp(sum(exo.foo)))/365





z = as.matrix(d[,exo.cols]);
beta0 <- c(-0.5,1.02,-0.3,0.2,0.5,0.6,0.7,0.8,0.9)
beta0 <- x$par

sourceCpp("funcs_discrete_loglik.cpp")

discrete_loglik_eigen_grad_cpp(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, pars = beta0)
pracma::grad(f = discrete_loglik_eigen_cpp, x0 = beta0, m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z)


log.lik.cov2(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, beta = beta0, states = states)
discrete_loglik_cpp(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, pars = beta0)
discrete_loglik_eigen_cpp(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, pars = beta0)


pracma::grad(f = log.lik.cov2, x0 = beta0, m = m, s1 = d$s1, s2 = d$s2, u = d$u/365, z = z, states = states)
pracma::grad(f = discrete_loglik_cpp, x0 = beta0, m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z)
pracma::grad(f = discrete_loglik_eigen_cpp, x0 = beta0, m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z)

eigenspace.U.grad(5,c(-1,3,4,2.5), 4) 
eigenspace_U_grad(5, c(-1,3,4,-2.5), 3)


s1 = d$`s-`; s2 = d$`s`; u = d$u/365; z = as.matrix(d[,exo.cols]);



x = optim( par = beta0, fn = discrete_loglik_cpp, gr = discrete_loglik_eigen_grad_cpp, m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, method = "BFGS") 
x = optim( par = beta0, fn = discrete_loglik_cpp, gr = NULL, m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, method = "BFGS") 
2 * (x$value + length(x$par))


beta0 <- c(-1.5,1.02,-0.3,0.4, 0.5,0.6,0.7,0.8,0.9)


x1 = optim( par = beta0, fn = log.lik.cov3, gr = NULL, m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, method = "BFGS") 
2 * (x$value + length(x$par))

x2 = optim( par = beta0, fn = log.lik.cov3, gr = log.lik.grad.cov3, m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, method = "BFGS") 
2 * (x$value + length(x$par))



log.lik.cov3(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, beta = beta0, states = states)




library("bench")
bench::mark(
  log.lik.cov(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, beta = beta0, states = states),
  log.lik.cov2(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, beta = beta0, states = states),
  discrete_loglik_cpp(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, pars = beta0),
  discrete_loglik_eigen_cpp(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = z, pars = beta0), 
  relative = T, check = F)


#########################
## Monte Carlo EM algo ##
#########################
beta0 <- c(rep(0.25, (m-1)), rep(0.1,length(exo.cols)))  
res4 <- MCEM(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, z = as.matrix(d[,exo.cols]), 
             beta0 = beta0, states = states, iters = 10, step.size = NULL, n.samples = NULL )



#saveRDS(res4, file = "MCEM_estimates.rds")
res4 <- readRDS("MCEM_estimates.rds")


#########################
### K-fold Cross Val. ###
#########################

k <- 10
rand.idx <- sample(x = 1:NROW(d),size = NROW(d),replace = F)
beta0 <- c(rep(1.25, (m-1)), rep(0.1,length(exo.cols)))  
pars <- matrix(NA, ncol = length(beta0), nrow = k)
start <- 1
for(i in 1:k){
  if(i == k){
    end <- NROW(d)
  } else {
    end <- start + floor(NROW(d)/k) 
  }
  pred.idx <- rand.idx[start:end]
  
  res3 <- mle.cov(m = m, s1 = d$`s-`[-pred.idx], s2 = d$s[-pred.idx], u = d$u[-pred.idx]/365, 
                  z = as.matrix(d[-pred.idx,exo.cols]), beta = beta0, states = states, 
                  A_param = NULL, pvalues.bin =  T, method = "BFGS", weight = NULL)
  
  pars[i,] <- res3$lambda 
  
  start <- end + 1
  print(i)
}



library(ggplot2)
library(tidyr)
require(gridExtra)
ESTIM <- as.data.frame(pars)
ESTIM[,1:4] <- exp(ESTIM[,1:4])/365
ESTIM$iteration <- 1:nrow(ESTIM)
data_frame_long1 <- gather(ESTIM[,c(1:4,10)], key = "trajectory", value = "value", -iteration)
data_frame_long2 <- gather(ESTIM[,c(5:9,10)], key = "trajectory", value = "value", -iteration)
text.size <- 16
p1 <- ggplot(data_frame_long1, aes(x = iteration, y = value, color = trajectory)) +
  geom_line(linewidth = 0.75) +
  labs(title = "",
       x = "Iteration",
       y = "Value",
       color = "Trajectory") + 
  theme(text = element_text(size = text.size, family = "serif"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        panel.grid.minor = element_line(color = "lightgray"),
        legend.key=element_blank(),legend.position = "bottom", 
        legend.text = element_text(size = text.size)) + 
  xlab("k-fold iterations") + ylab("MLE") + 
  scale_color_discrete(name = "Base transition rates", 
                       labels = c(expression(lambda["3,2B"]^{(0)}), 
                                  expression(lambda["2B,2A"]^{(0)}), 
                                  expression(lambda["2A,1"]^{(0)}), 
                                  expression(lambda["1,0"]^{(0)}))) 

p2 <- ggplot(data_frame_long2, aes(x = iteration, y = value, color = trajectory)) +
  geom_line(linewidth = 0.75) +
  labs(title = "",
       x = "Iteration",
       y = "Value",
       color = "Trajectory") + 
  theme(text = element_text(size = text.size, family = "serif"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        panel.grid.minor = element_line(color = "lightgray"),
        legend.key=element_blank(),legend.position = "bottom", 
        legend.text = element_text(size = text.size)) + 
  xlab("k-fold iterations") + ylab("MLE") + 
  scale_color_discrete(name = "Covariates", labels = c(expression(beta["Tonnage"]), 
                                                       expression(beta["Speed"]),
                                                       expression(beta["Profile"]),
                                                       expression(beta["Grade"]), 
                                                       expression(beta["Curve"])))

grid.arrange(p1, p2, nrow=1)


####################################
### Simulated model verification ###
####################################
d.sim = d
exo.cols <-  c("MBT.norm","invRad.norm")
beta0 <- c(0.9,1.0,1.1,1.2, 0.1, 0.5)  
for(i in 1:NROW(d)){
  foo <- sim.Markov.cov(m = m, states = states, beta = beta0, 
                        A_param = NULL, s1 = d.sim$`s-`[i], u = d$u[i]/365, z = as.matrix(d.sim[i,exo.cols]), 
                        n.samples = 1, step.size = 1/365)
  d.sim$`s`[i] <- foo[NROW(foo)] 
  print(i)
}

res1 <- mle.cov(m = m, s1 = d.sim$`s-`, s2 = d.sim$s, u = d.sim$u/365, 
                z = as.matrix(d.sim[,exo.cols]), beta = rep(0.1,6), states = states, 
                A_param = NULL, pvalues.bin =  T, method = "BFGS", weight = NULL)
res1



#################################
## Table descriptive rail stat ##
#################################
library(xtable)
#descriptive stats
count.track <- matrix(NA, nrow = length(track)+1, ncol = 4); count  <- 0
for(i in track){
  count <- count + 1
  idx <- !idx.remove & D$Track0 == i
  foo <- D[idx, ]
  count.track[count,1] <- sum(foo$Track == i)
  count.track[count,2] <- mean(foo$MBT )
  count.track[count,3] <- mean(foo$speed )
  count.track[count,4] <- mean(abs(foo$curve.rad), na.rm = T)
}
count.track[length(track)+1, 1] <- sum(count.track[1:length(track),1])
count.track[length(track)+1, 2] <- mean(count.track[1:length(track),2])
count.track[length(track)+1, 3] <- mean(count.track[1:length(track),3])
count.track[length(track)+1, 4] <- mean(count.track[1:length(track),4])


rownames(count.track) <- c(sapply(track, get_full_name), "Total")
xtable(count.track, type = "latex", digits =matrix(rep(c(0,0,1,1,1),length(track)+1), ncol = 5, byrow =T) )



#estimation table
res2 = res4
res2$lambda = res4$beta
n.digits <- 6
estim <- matrix(NA, nrow = m-1 + length(exo.cols), ncol = 3)
estim[,1] <- round(c(exp(res2$lambda[1:(m-1)])/365, 
                     res2$lambda[m:length(res2$lambda)]), n.digits)
estim[,2] <- paste0("[",
                    round(c(exp(res2$CI.LO[1:(m-1)])/365, 
                            res2$CI.LO[m:length(res2$lambda)]),n.digits),
                    ",",
                    round(c(exp(res2$CI.UP[1:(m-1)])/365, 
                            res2$CI.UP[m:length(res2$lambda)]),n.digits),
                    "]")
estim[,3] <- c(round(1/exp(res2$lambda[1:(m-1)])*365,1),
               round(apply(sweep((d[,exo.cols]), 2, res2$lambda[m:(m-1+length(exo.cols))], "*"), 
                           MARGIN = 2, quantile, probs = 0.5),2)) 

print(xtable(estim , type = "latex", display=c("s","s", "s","s")), math.style.exponents = TRUE)


foo <- apply(d[,exo.cols]* res2$lambda[m:(m-1+length(exo.cols))], 
             MARGIN = 2, quantile, probs = 0.5)

head(foo)

res2$lambda[m:(m-1+length(exo.cols))]  head(d[,exo.cols])

sweep(head(d[,exo.cols]), 2, res2$lambda[m:(m-1+length(exo.cols))], "*")
apply(
  sweep((d[,exo.cols]), 2, res3$lambda[m:(m-1+length(exo.cols))], "*"), 
  MARGIN = 2, quantile, probs = 0.5)


##########################################################
### TABLE FOR MEDIAN TIME TO TRANSITION FOR ALL TRACKS ###
##########################################################
library(stringr)
library(xtable)
#library(xlsx)
#library(openxlsx)
#install.packages("openxlsx")


base.foo1 <- exp(res3$lambda[1:(m-1)])
base.foo2 <- exp(res4$beta[1:(m-1)])

qprob = .5
#base.foo <- exp(res3$lambda[1:(m-1)])
Tracks <- unique(d$Track0)
ave.trans <- matrix(NA, ncol = 2*m, nrow = length(Tracks)+1); count <- 0
for(track in Tracks){ #loop through all tracks
  count <- count + 1
  idx <- d$Track0 == track #find track observations
  exo.foo1 <- apply(sweep((d[idx,exo.cols]), 2, res3$lambda[m:(m-1+length(exo.cols))], "*"), 
                    MARGIN = 2, quantile, probs = qprob)
  exo.foo2 <- apply(sweep((d[idx,exo.cols]), 2, res4$beta[m:(m-1+length(exo.cols))], "*"), 
                    MARGIN = 2, quantile, probs = qprob)
  
  tr1 <- base.foo1*exp(sum(exo.foo1)) / 365
  tr2 <- base.foo2*exp(sum(exo.foo2)) / 365
  
  ave.trans[count, ] <- c(round(1/tr2, 2), sum(round(1/tr2, 2)), 
                                     round(1/tr1, 2), sum(round(1/tr1, 2)))
}
exo.foo1 <- apply(sweep((d[,exo.cols]), 2, res3$lambda[m:(m-1+length(exo.cols))], "*"), 
                  MARGIN = 2, quantile, probs = qprob)
exo.foo2 <- apply(sweep((d[,exo.cols]), 2, res4$beta[m:(m-1+length(exo.cols))], "*"), 
                  MARGIN = 2, quantile, probs = qprob)

tr1 <- base.foo1*exp(sum(exo.foo1)) / 365
tr2 <- base.foo2*exp(sum(exo.foo2)) / 365
ave.trans[count+1, ] <- c(round(1/tr2, 2), sum(round(1/tr2, 2)), 
                        round(1/tr1, 2), sum(round(1/tr1, 2)))

rownames(ave.trans) <- c(sapply(Tracks, get_full_name), "Total")
colnames(ave.trans) <- rep(c("3 -> 2B", "2B -> 2A", "2A -> 1", "1 -> 0", "Sum"),2)


#write.xlsx(as.data.frame(ave.trans), "transition_rates_all_tracks.xlsx", rowNames = T )
ave.trans
foo <- round(ave.trans,1)
print(xtable( ave.trans, type = "latex",digits = 1), include.rownames=T)

####################################
#### MAKE MCEM convergence plot ####
####################################
library(ggplot2)
library(tidyr)
require(gridExtra)
ESTIM <- as.data.frame(res4$beta_iter)
ESTIM[,1:4] <- exp(ESTIM[,1:4])/365
ESTIM$iteration <- 1:nrow(ESTIM)
data_frame_long1 <- gather(ESTIM[,c(1:4,10)], key = "trajectory", value = "value", -iteration)
data_frame_long2 <- gather(ESTIM[,c(5:9,10)], key = "trajectory", value = "value", -iteration)
text.size <- 16
p1 <- ggplot(data_frame_long1, aes(x = iteration, y = value, color = trajectory)) +
  geom_line(linewidth = 0.75) +
  labs(title = "",
       x = "Iteration",
       y = "Value",
       color = "Trajectory") + 
  theme(text = element_text(size = text.size, family = "serif"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        panel.grid.minor = element_line(color = "lightgray"),
        legend.key=element_blank(),legend.position = "bottom", 
        legend.text = element_text(size = text.size)) + 
  xlab("MCEM iterations") + ylab("MLE") + 
  scale_color_discrete(name = "Base transition rates", 
                       labels = c(expression(lambda["3,2B"]^{(0)}), 
                                  expression(lambda["2B,2A"]^{(0)}), 
                                  expression(lambda["2A,1"]^{(0)}), 
                                  expression(lambda["1,0"]^{(0)}))) 

p2 <- ggplot(data_frame_long2, aes(x = iteration, y = value, color = trajectory)) +
  geom_line(linewidth = 0.75) +
  labs(title = "",
       x = "Iteration",
       y = "Value",
       color = "Trajectory") + 
  theme(text = element_text(size = text.size, family = "serif"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        panel.grid.minor = element_line(color = "lightgray"),
        legend.key=element_blank(),legend.position = "bottom", 
        legend.text = element_text(size = text.size)) + 
  xlab("MCEM iterations") + ylab("MLE") + 
  scale_color_discrete(name = "Covariates", labels = c(expression(beta["Tonnage"]), 
                                                       expression(beta["Speed"]),
                                                       expression(beta["Profile"]),
                                                       expression(beta["Grade"]), 
                                                       expression(beta["Curve"])))

grid.arrange(p1, p2, nrow=1)





#########################
### FORECASTING PLOTS ###
#########################
idx <- 2052
x <- sim.Markov.cov(m = m, states = (states), beta = res3$lambda, 
                    A_param = NULL, s1 = 2, u = 1500/365, z = as.matrix(d[idx,exo.cols]), 
                    n.samples = 100, step.size = 1/365)

#d[which(d$Track0 == "ofot"),]


marker_coordinates <- data.frame(
  x = c(0),  # Replace with the actual x coordinate
  y = c(2)   # Replace with the actual y coordinate
)
df <- as.data.frame(x)
df$time <- cumsum(1:NROW(df) / 365)
library(tidyr)
df_long <- gather(df[1:1500, ], key = "variable", value = "value", -time)
library(ggplot2)
text.size <- 16
# Plot
p1 <-ggplot(df_long, aes(x = time, y = value, color = variable)) +
  geom_line(size = 0.75) + 
  geom_point(data = marker_coordinates, aes(x = x, y = y), color = "black", size = 4) +  # Add point marker
  theme(text = element_text(size = text.size, family = "serif"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        panel.grid.minor = element_line(color = "lightgray"),
        legend.key=element_blank(),legend.position = "none"
  ) + 
  scale_y_discrete(limits = convert.to.num.inv(states)) + 
  xlab("Time [days]") + ylab("Defect class") + 
  annotate("text", x=min(df_long$time), y=marker_coordinates$y, size=text.size/3, label = "Observation", 
           family="serif", vjust = 1.1, hjust = -.10)



p2 <- forecast.state.probs2(m = m, states = convert.to.num.inv(states), z = as.matrix(d[idx,exo.cols]), 
                            class = "2B", pars = res3$lambda, span = 1500/365, position = round(d$pos[idx]), 
                            track = d$Track0[idx])

grid.arrange(p2[[2]], p1, ncol = 2)


d[idx, ]

D[D$pos == 17.11,]

################
### FIGURE 1 ###
################
load("defects_covs_base.RData"); D <- foo; rm(foo) #load dat

#########################
### Data preparation ####
#########################
states <- c(1,2,3,4,5) #define states
#states <- c(2,3,4,5)
m <- length(states) #number of states
track <- unique(D$Track) #define investigated tracks
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", 
              "invRad.norm")#, , "dist.to.station.trans") #covariate space
#idx <- D$Track0 %in% track #choose all tracks 
#idx <- Dtot$Track0 == "nord" #select single tracks

idx.remove <- 
  (D$`s-` == 4 & D$`s` == 4 & D$u > 300) | 
  (D$`s-` == 4 & D$`s` == 5 & D$u > 300) | 
  (D$`s-` == 5 & D$`s` == 5 & D$u > 300) 


d <- D[!idx.remove, c("pos", "s-", "s", "u", exo.cols, "Track0")] #define data set of interest



#histograms of u for different defect classes
library(ggplot2)
require(gridExtra)
convert.to.num.inv <- function(dat){
  x <- dat; 
  x[dat == 1] <- "3"
  x[dat == 2] <- "2B"
  x[dat == 3] <- "2A"
  x[dat == 4] <- "1"
  x[dat == 5] <- "0"
  return(x)
}
Dtot <- d
Dtot$`s-` <- convert.to.num.inv(Dtot$`s-`)
Dtot$`s` <- convert.to.num.inv(Dtot$`s`)
states <- c("3", "2B", "2A", "1", "0") #define states
plist <- list(); text.size <- 14; text.size2 = 5; 
icount <- 0; jcount <- 0; pcount <- 0; sum = 0;
for(i in states){
  icount <- icount + 1
  for(j in states){
    jcount <- jcount + 1
    if(jcount >= icount){
      pcount <- pcount + 1
      idx <- Dtot$`s-` == i & Dtot$s == j
      lab <- bquote(.(i) ~ symbol('\256') ~ .(j))
      dtemp = Dtot[idx, ]
      lab2 <- paste("#Obs. = ", NROW(dtemp))
      sum = sum + NROW(dtemp);
      p <- ggplot(dtemp, aes(x=u)) + 
        geom_histogram(aes(y =  after_stat(count / sum(count)) ), binwidth = 100,color="black", fill="grey") + 
        theme(text = element_text(size = text.size, family = "serif"),
              panel.background = element_rect(fill = "white", color = "white"),
              panel.grid.minor = element_line(color = "lightgray"),
              plot.background = element_rect(color = "black"), 
              #axis.title.x=element_blank(),
              #axis.title.y=element_blank()
              #plot.margin = margin(5, 5, 5, 5, "mm")
              ) + 
        scale_y_continuous(labels = scales::percent) +
        scale_x_continuous(breaks = 
                             c(300,600,900,1200), limits = c(0,1200) ) +
        labs(x = "Trans. times [days]", y = "Freqency") + 
        annotate("text", x=Inf, y=Inf,size=text.size2, label = lab, family="serif", vjust = 1, hjust = 1) +
        annotate("text", x=Inf, y=Inf,size=text.size2, label = lab2, family="serif", vjust = 2.5, hjust = 1) 
      
        # if(icount < 6){
        #   scale_x_continuous(breaks = 
        #                        #round(c(max(Dtot$u[idx])/4, max(Dtot$u[idx])/4*2, max(Dtot$u[idx])/4*3 )/25,0)*25
        #                        c(300,600,900,1200)
        #   )  
        # } else {
        # scale_x_continuous(breaks = 
        #                      #round(c(max(Dtot$u[idx])/4, max(Dtot$u[idx])/4*2, max(Dtot$u[idx])/4*3 )/25,0)*25
        #                       c(90,180,270)
        #                      )}
      plist[[pcount]] <- p
    } else {
      pcount <- pcount + 1
      plist[[pcount]] <- ggplot(data = data.frame(), aes()) + theme_void() 
    }
    
  }
  jcount <- 0
}
p1 <- grid.arrange(grobs = plist, ncol = m, nrow = m,widths = rep(1, m), heights = rep(1, m))
sum

ggsave("hists_u.pdf", p1, width = 15, height = 10, units = "in")



dtemp = Dtot[idx, ]
p <- ggplot(dtemp, aes(x=u)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))),color="black", fill="grey") + 
  theme(text = element_text(size = text.size, family = "serif"),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_line(color = "lightgray"),
        plot.background = element_rect(color = "black"),
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks =  c(300,600,900,1200), limits = c(0,1200) ) +
  labs(x = "Trans. times [days]", y = "Freqency") + 
  annotate("text", x=Inf, y=Inf,size=text.size2, label = lab, family="serif", vjust = 1, hjust = 1) 
p






#############################
### Plot of survival func ###
#############################
library(expm)

beta0 <- c(rep(0.25, (m-1)), rep(0.1,length(exo.cols)))  
res3 <- mle.cov(m = m, s1 = d$`s-`, s2 = d$s, u = d$u/365, 
                z = as.matrix(d[,exo.cols]), beta = beta0, 
                pvalues.bin =  T, method = "BFGS")


# res3
# base.foo <- 1/exp(res3$lambda[1:(m-1)])*365
# exo.foo <- apply(sweep((d[,exo.cols]), 2, res3$lambda[m:(m-1+length(exo.cols))], "*"), 
#                  MARGIN = 2, quantile, probs = 0.5)
# 
# 
# base.foo
# base.foo/exp(sum(exo.foo))
# sum(base.foo/exp(sum(exo.foo)))/365

lambda = exp(res3$lambda[1:(m-1)])
exo = res3$lambda[m:(m-1+length(exo.cols))]

#begin loop
results <- data.frame()
for(i in 1:NROW(d)){
  impact = exp(sum(d[i, exo.cols] * exo))
  A = construct.A3(m,lambda) * impact
  Q = A[1:(m-1), 1:(m-1)]
  
  time <- seq(0, 40, length.out = 81)
  pi <- c(1, 0, 0, 0)
  
  survival_function <- function(t) {
    exp_Qt <- expm::expm(Q * t)
    1 - pi %*% exp_Qt %*% rep(1, nrow(Q))
  }
  
  survival_values <- sapply(time, survival_function)
  temp_df <- data.frame(time = time, survival = survival_values, Q_id = as.factor(i))
  results <- rbind(results, temp_df)
  
}

#choose an appriopiate subset of survival functions evenly distributed over degradation speed
foo.res = results[which(results$time == 10),]
foo.res <- (foo.res[order(foo.res$survival),])
foo.idx = floor(sort(c(6,10,20,40,seq(1,NROW(foo.res), length.out = 100))))
chosenQ = foo.res[foo.idx,"Q_id"] 
idx = results$Q_id %in% chosenQ


library(dplyr)

# Compute average survival function
avg_survival <- results %>%
  group_by(time) %>%
  summarize(avg_survival = mean(survival, na.rm = TRUE))


#create shading data
curve1_id <- chosenQ[1]
curve2_id <- chosenQ[length(chosenQ)]
curve1_data <- results[results$Q_id == curve1_id, ]
curve2_data <- results[results$Q_id == curve2_id, ]
shading_data <- merge(curve1_data, curve2_data, by = "time", suffixes = c("_1", "_2"))
library(dplyr)
shading_data <- shading_data %>%
  rename(ymin = survival_2, ymax = survival_1)


text.size = 12
p <- ggplot() +
  geom_line(data = results[idx,], aes(x = time, y = survival, color = Q_id), size = 0.001, alpha=0.8) +
  geom_ribbon(data = shading_data, aes(x = time, ymin = ymin, ymax = ymax), fill = "gray", alpha = 0.5) +
  geom_line(data = avg_survival, aes(x = time, y = avg_survival), color = "#F8766D", size = 1.5, linetype = "dashed") +
  labs(x = "Time [years]",
       y = "Probability") +
  guides(color = "none") +  # Remove the legend
  theme_minimal() +
  scale_color_manual(values = rep("black", length(unique(results$Q_id[idx]))))+
  theme(text = element_text(size = text.size, family = "serif"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        panel.grid.minor = element_line(color = "lightgray"),
        legend.key=element_blank(),legend.position = "none"
  )

p

ggsave("survival_funcs.pdf", p, width = 18, height = 10, units = "cm")


############################
### Data scenarios notes ###
############################

##############
# GPS_repair #
##############
#N = 257, 
#large p-values, too uncertain data foundation
#Conclusion: Not usable for this type of analysis

##############
### repair ###
##############
#N = 1735, 
#state 3 not included
#significant: MBT (+), profile (-), grade (+), and curvature (+) 
#1/exp(res3$lambda[1:(m-1)])*365 = 758.4775  740.0473 4454.6413


##############
#### base ####
##############

#CASE: all data
#N = 4272
#base.foo = 703.5859  394.6391 1278.1599 3030.1093
#base.foo/exp(sum(exo.foo)) = 659.9665  370.1731 1198.9193 2842.2551
#covs significant and interpretable

#CASE: remove 1-1 and 0-0
#N = 3530
#base.foo = 715.1533  402.8296 1308.6908 1136.9907
#base.foo/exp(sum(exo.foo)) = 663.9797  374.0047 1215.0459 1055.6321
#covs significant and interpretable


#CASE: remove 1-1 and 0-0 and 1-0
#N = 3446
#base.foo = 700.1811  392.0427 1268.6674 1790.1393
#base.foo/exp(sum(exo.foo)) = 669.8716  375.0719 1213.7491 1712.6474
#covs significant and interpretable


#CASE: remove 1-1 and 0-0 and 1-0 over 650 days
#N = 3930
#base.foo = 712.4657  400.2408 1295.6633 1864.6132
#base.foo/exp(sum(exo.foo)) = 661.8524  371.8078 1203.6200 1732.1519
#covs significant and interpretable


#CASE: remove 1-1 and 0-0 and 1-0 over 300 days
#N = 3658
#base.foo = 708.5189  398.0408 1288.7816 1481.4379
#base.foo/exp(sum(exo.foo)) = 665.2714  373.7447 1210.1152 1391.0119
#base.foo/exp(sum(exo.foo)) v2 = 418.6170 235.1760 761.4558 875.2837

#covs significant and interpretable






