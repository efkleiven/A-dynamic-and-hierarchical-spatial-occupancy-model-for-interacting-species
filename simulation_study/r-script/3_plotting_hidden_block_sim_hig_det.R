########################################################################################################################################
##                             Plotting the model output of  the Spatial Dynamic two-species occupancy hidden block model             ##
##                                  by EFK                                                                                            ##
########################################################################################################################################

rm(list=ls())

# Call jags(and other packages)
library(jagsUI)

# set working directory
setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models")    # for norpec-server
#setwd("H:/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models")  # for laptop

#load model
setwd("./hidden_block_sim/model_output")
load("mod_hidden_sim_hig_det.rda")

mod_hidden_sim_hig_det[[2]]$summary
#######################
# Model Convergance  ##
#######################

#rhat
mod_hidden_sim_hig_det[[1]]$Rhat[1:19]
hist(mod_hidden_sim_hig_det[[1]]$Rhat$z)
hist(mod_hidden_sim_hig_det[[1]]$Rhat$x)

#neff
mod_hidden_sim_hig_det[[1]]$n.eff[1:19]
hist(mod_hidden_sim_hig_det[[1]]$n.eff$z)
hist(mod_hidden_sim_hig_det[[1]]$n.eff$x)

# number of parameters with rhat>1.1
length(mod_hidden_sim_hig_det[[1]]$summary[,8][mod_hidden_sim_hig_det[[1]]$summary[,8]>1.1])
max(mod_hidden_sim_hig_det[[1]]$summary[,8], na.rm = T)

# traceplots
traceplot(mod_hidden_sim_mid_det_1)

traceplot(mod_hidden_sim_hig_det[[1]], parameters = "z[4,1,1]")
traceplot(mod_hidden_ko_s1_70_nt20, parameters = "GamA")
traceplot(mod_hidden_ko_s1_70_nt20, parameters = "gamB")
traceplot(mod_hidden_ko_s1_70_nt20, parameters = "GamB")

plot(mod_hidden_ko_s1_70$sims.list$gamB[1:(length(mod_hidden_ko_s1_70$sims.list$gamB)/4)], type="l")
plot(mod_hidden_ko_s1_70$sims.list$gamB[(length(mod_hidden_ko_s1_70$sims.list$gamB)/4):(2*length(mod_hidden_ko_s1_70$sims.list$gamB)/4)], type="l")
plot(mod_hidden_ko_s1_70$sims.list$gamB[(2*length(mod_hidden_ko_s1_70$sims.list$gamB)/4):(3*length(mod_hidden_ko_s1_70$sims.list$gamB)/4)], type="l")
plot(mod_hidden_ko_s1_70$sims.list$gamB[(3*length(mod_hidden_ko_s1_70$sims.list$gamB)/4):(length(mod_hidden_ko_s1_70$sims.list$gamB))], type="l")

###############################
# Model summary
###############################

# make a list with estimated parameter values from all 50 simulations 
dat<-c()

for (i in 1:length(mod_hidden_sim_hig_det)){
  dat$GamA[i] <- mod_hidden_sim_hig_det[[i]]$mean$GamA
  dat$GamB[i] <- mod_hidden_sim_hig_det[[i]]$mean$GamB
  dat$GamAB[i] <- mod_hidden_sim_hig_det[[i]]$mean$GamAB
  dat$GamBA[i] <- mod_hidden_sim_hig_det[[i]]$mean$GamBA
  dat$EpsA[i] <- mod_hidden_sim_hig_det[[i]]$mean$EpsA
  dat$EpsB[i] <- mod_hidden_sim_hig_det[[i]]$mean$EpsB
  dat$EpsAB[i] <- mod_hidden_sim_hig_det[[i]]$mean$EpsAB
  dat$EpsBA[i] <- mod_hidden_sim_hig_det[[i]]$mean$EpsBA
  dat$gamA[i] <- mod_hidden_sim_hig_det[[i]]$mean$gamA
  dat$gamB[i] <- mod_hidden_sim_hig_det[[i]]$mean$gamB
  dat$gamAB[i] <- mod_hidden_sim_hig_det[[i]]$mean$gamAB
  dat$gamBA[i] <- mod_hidden_sim_hig_det[[i]]$mean$gamBA
  dat$epsA[i] <- mod_hidden_sim_hig_det[[i]]$mean$epsA
  dat$epsB[i] <- mod_hidden_sim_hig_det[[i]]$mean$epsB
  dat$epsAB[i] <- mod_hidden_sim_hig_det[[i]]$mean$epsAB
  dat$epsBA[i] <- mod_hidden_sim_hig_det[[i]]$mean$epsBA
  dat$pA[i] <- mod_hidden_sim_hig_det[[i]]$mean$pA
  dat$pB[i] <- mod_hidden_sim_hig_det[[i]]$mean$pB
  }

# the same for initial state (psi)
psi<- array(NA, dim=c(length(mod_hidden_sim_hig_det),4,3))

for(i in 1: length(mod_hidden_sim_hig_det)){
  for(j in 1:3){
  psi[i,,j] <- mod_hidden_sim_hig_det[[i]]$mean$psi[,j]
}}

# plot estimated block parameters

low <- sapply(dat[1:8],min)
high <- sapply(dat[1:8],max)
mu <- sapply(dat[1:8],mean)
true <- c(0.5,0.1,0.05,0.4,0.05,0.6,0.5,0.2)
x <- c(1,2,3,4,5,6,7,8)

png("modelpreformance_sim_hig_det_block.png", width=920, height=540)

plot(NULL, type="n", xlab="", ylab="", xlim=c(0, 8), ylim=c(0, 1), axes = F, main="Hidden block SimMidDet Block level parameters")
axis(side=1, c(1,2,3,4,5,6,7,8), c("GamA","GamB","GamAB","GamBA","EpsA","EpsB","EpsAB","EpsBA"), cex.axis=1.5)
axis(side=2, c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), cex.axis=2)
points(mu, pch=16, col=1, cex=2)
arrows(x0=x, y0=low, x1=x, y1=high, angle=90, code=3, length=0.2, lwd=2)
points(true, pch="-", col=2,cex=4)

#legend("topleft",c("Model estimates", "True parameter value"), pch=16:17, col=c("black","red"), cex=1.2)
dev.off()

# plot estimated site parameters with true simulated values
low <- sapply(dat[9:18],min)
high <- sapply(dat[9:18],max)
mu <- sapply(dat[9:18],mean)
true <- c(0.5,0.3,0.1,0.7,0.3,0.8,0.9,0.1,0.9,0.8)
x <- c(1,2,3,4,5,6,7,8,9,10)



png("modelpreformance_sim_hig_det_site.png", width=1060, height=540)

plot(NULL, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 1), axes = F, main="Hiddenblock SimMidDet Site level parameters")
axis(side=1, c(1,2,3,4,5,6,7,8,9,10), c("gamA","gamB","gamAB","gamBA","epsA","epsB","epsAB","epsBA","pA","pB"), cex.axis=1.5)
axis(side=2, c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), cex.axis=2)
points(mu, pch=16, col=1, cex=2)
arrows(x0=x, y0=low, x1=x, y1=high, angle=90, code=3, length=0.2, lwd=2)
points(true, pch="-", col=2, cex=4)

#legend("topleft",c("Model estimates", "True parameter value"), pch=16:17, col=c("black","red"), cex=1.2)
dev.off()

# plot estimated initial season occupancy probability and detection probabilities with true simulated values
setwd("../plot")

dim(psi)
psi_block <- apply(psi,c(1,3),mean)
psi_block

low <- apply(psi_block,2,min)
high <- apply(psi_block,2,max)
mu <- apply(psi_block,2,mean)
x <- c(1,2,3)
sim_data <- c(0.3,0.2,0.1)

png("modelpreformance_sim_hig_det_psi.png", width=920, height=540)

plot(NULL, type="n", xlab="", ylab="", xlim=c(0, 3), ylim=c(0, 1), axes = F, main="Initial occupnacy parameters")
axis(side=1, c(1,2,3), c("psiA","psiB","psiAB"), cex.axis=1.5)
axis(side=2, c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), cex.axis=2)
points(mu, pch=16, col=1, cex=2)
arrows(x0=x, y0=low, x1=x, y1=high, angle=90, code=3, length=0.2, lwd=2)
points(sim_data, pch=17, col="red", cex=2)
legend("topleft",c("Model estimates", "True parameter value"), pch=16:17, col=c("black","red"), cex=1.2)
dev.off()



