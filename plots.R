


rm(list = ls()) #clear memory

### Read data and define state space and covariates ###
d = read.csv("defect_data.csv")
states <- c(1,2,3,4,5) #define states
m <- length(states) #number of states
track <- unique(d$Track) #define investigated tracks
exo.cols <- c("MBT.norm","speed.norm","profil.norm", "steel.norm", "invRad.norm")

### Read libaries ###
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")




#############################
### Plot of survival func ###
#############################
library(expm)


z = as.matrix(d[, exo.cols])
beta0 <- c(c(0.1,0.2,0.3,0.4),
           rep(0.1,length(exo.cols)))

res3 =
  optim( par = beta0, fn = MJP_score, m = m, s1 = d$s1, s2 = d$s2, u = d$t, z = z,
         generator="gerlang", covs_bin = T, likelihood_bin = T, rps_bin = F, brier_bin = F,
         transient_dist_method = "eigenvalue_decomp",
         method = "BFGS", control = list(maxit = 1000))

lambda = exp(res3$par[1:(m-1)])
exo = res3$par[m:(m-1+length(exo.cols))]

#begin loop
results <- data.frame()
for(i in 1:NROW(d)){
  impact = exp(sum(d[i, exo.cols] * exo))
  A = make_A1(m,lambda) * impact
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
library(ggplot2)
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
       y = "Failure probability") +
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




####################################################
### histograms of u for different defect classes ###
####################################################
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
Dtot$`s-` <- convert.to.num.inv(Dtot$s1)
Dtot$`s` <- convert.to.num.inv(Dtot$s2)
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


