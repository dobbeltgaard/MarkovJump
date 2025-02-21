rm(list = ls()) #clear memory

### Read data and define state space and covariates ###
d = read.csv("defect_data.csv")


change_in_growth_rate = function(data, beta, val1, val2){
  
  var1_hat = (val1 - min(data)) / (max(data) - min(data))
  var2_hat = (val2 - min(data)) / (max(data) - min(data))
  
  t1 = exp(var1_hat*beta)^-1
  t2 = exp(var2_hat*beta)^-1
  
  return((t1-t2)/t2*100)
}

change_in_growth_rate2 = function(data, beta, val1, val2){
  
  #var1_hat = (val1 - min(data)) / (max(data) - min(data))
  #var2_hat = (val2 - min(data)) / (max(data) - min(data))
  
  t1 = exp(val1*beta)^-1
  t2 = exp(val2*beta)^-1
  
  return((t1-t2)/t2*100)
}


#Tonnage analysis
as.numeric(quantile(d$MBT))
as.numeric(quantile(d$MBT, probs = c(0.2,0.5,0.8)))
change_in_growth_rate(data = d$MBT, beta = 4.7*10^-1,
                      as.numeric(quantile(d$MBT,probs = 0.2)) , 
                      as.numeric(quantile(d$MBT,probs = 0.8)))
change_in_growth_rate(data = d$MBT, beta = 4.4*10^-1,2,8)
change_in_growth_rate(data = d$MBT, beta = 4.7*10^-1,2,8)

                      
#Speed analysis
as.numeric(quantile(d$speed, probs = c(0.2,0.5,0.8)))
change_in_growth_rate(data = d$speed, beta = 3.2*10^-1,
                      as.numeric(quantile(d$speed,probs = 0.2)) , 
                      as.numeric(quantile(d$speed,probs = 0.8)))
change_in_growth_rate(data = d$speed, beta = 3.2*10^-1,70,110)
change_in_growth_rate(data = d$speed, beta = 3.4*10^-1,70,110)


#Curvature analysis
d[50:100, c("curve.rad","invRad.norm")]
change_in_growth_rate2(data = d$invRad.norm, beta = 4.6*10^-1,0,0.15) #straight to 1200m
change_in_growth_rate2(data = d$invRad.norm, beta = 4.6*10^-1,0.15, 0.63) #

#profile analysis
d[50:100, c("profil","profil.norm")]
change_in_growth_rate2(data = d$profil.norm, beta = -7.2*10^-1,0.65,1) #straight to 1200m
change_in_growth_rate2(data = d$profil.norm, beta = -6.9*10^-1,0.65,1) #straight to 1200m

names(d)





