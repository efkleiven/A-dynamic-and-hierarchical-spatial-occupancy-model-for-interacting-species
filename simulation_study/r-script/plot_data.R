
rm(list=ls())

setwd("H:/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models")   
setwd("./hidden_block_sim/data")
dir()

load("simdata_50set_50seas_8b_y_hig_det.rda")
load("simdata_50set_50seas_8b_y_hig_occ.rda")
load("simdata_50set_50seas_8b_y_low_det.rda") 
load("simdata_50set_50seas_8b_y_low_occ.rda") 
load("simdata_50set_50seas_8b_y_mid_det.rda")
load("simdata_50set_50seas_8b_y_mid_occ.rda")

load("simdata_50set_50seas_8b_z_hig_det.rda")
load("simdata_50set_50seas_8b_z_hig_occ.rda")
load("simdata_50set_50seas_8b_z_low_det.rda") 
load("simdata_50set_50seas_8b_z_low_occ.rda") 
load("simdata_50set_50seas_8b_z_mid_det.rda")
load("simdata_50set_50seas_8b_z_mid_occ.rda")

library(abind)

# merge blocks
#for obs
y_hd <- abind(y_hig_det[1,,1,,],y_hig_det[1,,2,,],y_hig_det[1,,3,,],y_hig_det[1,,4,,],along=1)
y_ho <- abind(y_hig_occ[1,,1,,],y_hig_occ[1,,2,,],y_hig_occ[1,,3,,],y_hig_occ[1,,4,,],along=1)
y_md <- abind(y_mid_det[1,,1,,],y_mid_det[1,,2,,],y_mid_det[1,,3,,],y_mid_det[1,,4,,],along=1)
y_mo <- abind(y_mid_occ[1,,1,,],y_mid_occ[1,,2,,],y_mid_occ[1,,3,,],y_mid_occ[1,,4,,],along=1)
y_ld <- abind(y_low_det[1,,1,,],y_low_det[1,,2,,],y_low_det[1,,3,,],y_low_det[1,,4,,],along=1)
y_lo <- abind(y_low_occ[1,,1,,],y_low_occ[1,,2,,],y_low_occ[1,,3,,],y_low_occ[1,,4,,],along=1)

#for true state
z_hd <- abind(z_hig_det[1,,1,],z_hig_det[1,,2,],z_hig_det[1,,3,],z_hig_det[1,,4,],along=1)
z_ho <- abind(z_hig_occ[1,,1,],z_hig_occ[1,,2,],z_hig_occ[1,,3,],z_hig_occ[1,,4,],along=1)
z_md <- abind(z_mid_det[1,,1,],z_mid_det[1,,2,],z_mid_det[1,,3,],z_mid_det[1,,4,],along=1)
z_mo <- abind(z_mid_occ[1,,1,],z_mid_occ[1,,2,],z_mid_occ[1,,3,],z_mid_occ[1,,4,],along=1)
z_ld <- abind(z_low_det[1,,1,],z_low_det[1,,2,],z_low_det[1,,3,],z_low_det[1,,4,],along=1)
z_lo <- abind(z_low_occ[1,,1,],z_low_occ[1,,2,],z_low_occ[1,,3,],z_low_occ[1,,4,],along=1)


## reduce to N or not N
# for obs
N_obs_hd <- array(NA, dim=c(48,50))
N_obs_ho <- array(NA, dim=c(48,50))
N_obs_md <- array(NA, dim=c(48,50))
N_obs_mo <- array(NA, dim=c(48,50))
N_obs_ld <- array(NA, dim=c(48,50))
N_obs_lo <- array(NA, dim=c(48,50))


for(t in 1:50){
  for(i in 1:48){
    N_obs_hd[i,t] <- ifelse(is.element(2, y_hd[i,,t])==T | is.element(4, y_hd[i,,t])==T, 1,0)
    N_obs_ho[i,t] <- ifelse(is.element(2, y_ho[i,,t])==T | is.element(4, y_ho[i,,t])==T, 1,0) 
    N_obs_md[i,t] <- ifelse(is.element(2, y_md[i,,t])==T | is.element(4, y_md[i,,t])==T, 1,0) 
    N_obs_mo[i,t] <- ifelse(is.element(2, y_mo[i,,t])==T | is.element(4, y_mo[i,,t])==T, 1,0) 
    N_obs_ld[i,t] <- ifelse(is.element(2, y_ld[i,,t])==T | is.element(4, y_ld[i,,t])==T, 1,0) 
    N_obs_lo[i,t] <- ifelse(is.element(2, y_lo[i,,t])==T | is.element(4, y_lo[i,,t])==T, 1,0) 
  }}

# for true
N_true_hd <- array(NA, dim=c(48,50))
N_true_ho <- array(NA, dim=c(48,50))
N_true_md <- array(NA, dim=c(48,50))
N_true_mo <- array(NA, dim=c(48,50))
N_true_ld <- array(NA, dim=c(48,50))
N_true_lo <- array(NA, dim=c(48,50))


for(t in 1:50){
  for(i in 1:48){
    N_true_hd[i,t] <- ifelse(is.element(2, z_hd[i,t])==T | is.element(4, z_hd[i,t])==T, 1,0)
    N_true_ho[i,t] <- ifelse(is.element(2, z_ho[i,t])==T | is.element(4, z_ho[i,t])==T, 1,0) 
    N_true_md[i,t] <- ifelse(is.element(2, z_md[i,t])==T | is.element(4, z_md[i,t])==T, 1,0) 
    N_true_mo[i,t] <- ifelse(is.element(2, z_mo[i,t])==T | is.element(4, z_mo[i,t])==T, 1,0) 
    N_true_ld[i,t] <- ifelse(is.element(2, z_ld[i,t])==T | is.element(4, z_ld[i,t])==T, 1,0) 
    N_true_lo[i,t] <- ifelse(is.element(2, z_lo[i,t])==T | is.element(4, z_lo[i,t])==T, 1,0) 
  }}

## reduce states to P or not P
#for obs
P_obs_hd <- array(NA, dim=c(48,50)) # make empty arrays with correct dimentions
P_obs_ho <- array(NA, dim=c(48,50))
P_obs_md <- array(NA, dim=c(48,50))
P_obs_mo <- array(NA, dim=c(48,50))
P_obs_ld <- array(NA, dim=c(48,50))
P_obs_lo <- array(NA, dim=c(48,50))

for(t in 1:50){
  for(i in 1:48){
    P_obs_hd[i,t] <- ifelse(is.element(3, y_hd[i,,t])==T | is.element(4, y_hd[i,,t])==T, 1,0)
    P_obs_ho[i,t] <- ifelse(is.element(3, y_ho[i,,t])==T | is.element(4, y_ho[i,,t])==T, 1,0)
    P_obs_md[i,t] <- ifelse(is.element(3, y_md[i,,t])==T | is.element(4, y_md[i,,t])==T, 1,0)
    P_obs_mo[i,t] <- ifelse(is.element(3, y_mo[i,,t])==T | is.element(4, y_mo[i,,t])==T, 1,0)
    P_obs_ld[i,t] <- ifelse(is.element(3, y_ld[i,,t])==T | is.element(4, y_ld[i,,t])==T, 1,0)
    P_obs_lo[i,t] <- ifelse(is.element(3, y_lo[i,,t])==T | is.element(4, y_lo[i,,t])==T, 1,0)
  }}

#for true
P_true_hd <- array(NA, dim=c(48,50)) # make empty arrays with correct dimentions
P_true_ho <- array(NA, dim=c(48,50))
P_true_md <- array(NA, dim=c(48,50))
P_true_mo <- array(NA, dim=c(48,50))
P_true_ld <- array(NA, dim=c(48,50))
P_true_lo <- array(NA, dim=c(48,50))

for(t in 1:50){
  for(i in 1:48){
    P_true_hd[i,t] <- ifelse(is.element(3, z_hd[i,t])==T | is.element(4, z_hd[i,t])==T, 1,0)
    P_true_ho[i,t] <- ifelse(is.element(3, z_ho[i,t])==T | is.element(4, z_ho[i,t])==T, 1,0)
    P_true_md[i,t] <- ifelse(is.element(3, z_md[i,t])==T | is.element(4, z_md[i,t])==T, 1,0)
    P_true_mo[i,t] <- ifelse(is.element(3, z_mo[i,t])==T | is.element(4, z_mo[i,t])==T, 1,0)
    P_true_ld[i,t] <- ifelse(is.element(3, z_ld[i,t])==T | is.element(4, z_ld[i,t])==T, 1,0)
    P_true_lo[i,t] <- ifelse(is.element(3, z_lo[i,t])==T | is.element(4, z_lo[i,t])==T, 1,0)
  }}

## make mean occupancy
#for obs
mean.n_obs_hd <- apply(N_obs_hd,c(2),mean) 
mean.p_obs_hd <- apply(P_obs_hd,c(2),mean)

mean.n_obs_ho <- apply(N_obs_ho,c(2),mean) 
mean.p_obs_ho <- apply(P_obs_ho,c(2),mean)

mean.n_obs_md <- apply(N_obs_md,c(2),mean) 
mean.p_obs_md <- apply(P_obs_md,c(2),mean)

mean.n_obs_mo <- apply(N_obs_mo,c(2),mean) 
mean.p_obs_mo <- apply(P_obs_mo,c(2),mean)

mean.n_obs_ld <- apply(N_obs_ld,c(2),mean) 
mean.p_obs_ld <- apply(P_obs_ld,c(2),mean)

mean.n_obs_lo <- apply(N_obs_lo,c(2),mean) 
mean.p_obs_lo <- apply(P_obs_lo,c(2),mean)

#for true
mean.n_true_hd <- apply(N_true_hd,c(2),mean) 
mean.p_true_hd <- apply(P_true_hd,c(2),mean)

mean.n_true_ho <- apply(N_true_ho,c(2),mean) 
mean.p_true_ho <- apply(P_true_ho,c(2),mean)

mean.n_true_md <- apply(N_true_md,c(2),mean) 
mean.p_true_md <- apply(P_true_md,c(2),mean)

mean.n_true_mo <- apply(N_true_mo,c(2),mean) 
mean.p_true_mo <- apply(P_true_mo,c(2),mean)

mean.n_true_ld <- apply(N_true_ld,c(2),mean) 
mean.p_true_ld <- apply(P_true_ld,c(2),mean)

mean.n_true_lo <- apply(N_true_lo,c(2),mean) 
mean.p_true_lo <- apply(P_true_lo,c(2),mean)

# plot observation per week
setwd("../plot")

png("obsTrends_all_sim.png", width=1640, height=480)
par(mfrow=c(2,3))
par(mar = c(5,7,4,2) + 0.1)

plot(mean.n_obs_hd, type="l", col=1, ylim=c(0,1), xlab="Seasons", ylab="observed occupancy", main="Trend in observations for HD", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_obs_hd, type="l", col=2, lwd=2)
lines(mean.n_true_hd, type="l", lty=2, col=3, lwd=1)
lines(mean.p_true_hd, type="l", lty=2, col=4, lwd=1)
legend("topright", legend=c("Prey", "Predator"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

plot(mean.n_obs_md, type="l", col=1, ylim=c(0,1), xlab="Seasons", ylab="observed occupancy", main="Trend in observations for MD", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_obs_md, type="l", col=2, lwd=2)
lines(mean.n_true_md, type="l", lty=2, col=3, lwd=1)
lines(mean.p_true_md, type="l", lty=2, col=4, lwd=1)
legend("topright", legend=c("Prey", "Predator"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

plot(mean.n_obs_ld, type="l", col=1, ylim=c(0,1), xlab="Seasons", ylab="observed occupancy", main="Trend in observations for LD", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_obs_ld, type="l", col=2, lwd=2)
lines(mean.n_true_ld, type="l", lty=2, col=3, lwd=1)
lines(mean.p_true_ld, type="l", lty=2, col=4, lwd=1)
legend("topright", legend=c("Prey", "Predator"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

plot(mean.n_obs_ho, type="l", col=1, ylim=c(0,1), xlab="Seasons", ylab="observed occupancy", main="Trend in observations for HO", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_obs_ho, type="l", col=2, lwd=2)
lines(mean.n_true_ho, type="l", col=3, lwd=2, lty=2)
lines(mean.p_true_ho, type="l", col=4, lwd=2, lty=2)
legend("topright", legend=c("Prey", "Predator"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

plot(mean.n_obs_mo, type="l", col="blue", ylim=c(0,1), xlab="Seasons", ylab="observed occupancy", main="Trend in observations for MO", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_obs_mo, type="l", col="red", lwd=2 )
lines(mean.n_true_mo, type="l", col="red", lwd=2, lty=2 )
lines(mean.p_true_mo, type="l", col="red", lwd=2, lty=2 )
legend("topright", legend=c("Prey", "Predator"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

plot(mean.n_obs_lo, type="l", col="blue", ylim=c(0,1), xlab="Seasons", ylab="observed occupancy", main="Trend in observations for LO", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_obs_lo, type="l", col="red", lwd=2)
lines(mean.n_true_lo, type="l", col="red", lwd=2, lty=2 )
lines(mean.p_true_lo, type="l", col="red", lwd=2, lty=2 )
legend("topright", legend=c("Prey", "Predator"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)


dev.off()