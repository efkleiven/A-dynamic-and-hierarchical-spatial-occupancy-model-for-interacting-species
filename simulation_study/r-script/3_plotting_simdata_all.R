setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models/hidden_block_sim/data")

# import packages
library(abind)

load("simdata_50set_50seas_y_hig_det.rda")
load("simdata_50set_50seas_y_hig_occ.rda")
load("simdata_50set_50seas_y_mid_det.rda")
load("simdata_50set_50seas_y_mid_occ.rda")
load("simdata_50set_50seas_y_low_det.rda")
load("simdata_50set_50seas_y_low_occ.rda")

# merge blocks
y_hd <- abind(y_hig_det[1,,1,,],y_hig_det[1,,2,,],y_hig_det[1,,3,,],y_hig_det[1,,4,,],along=1)
y_ho <- abind(y_hig_occ[1,,1,,],y_hig_occ[1,,2,,],y_hig_occ[1,,3,,],y_hig_occ[1,,4,,],along=1)
y_md <- abind(y_mid_det[1,,1,,],y_mid_det[1,,2,,],y_mid_det[1,,3,,],y_mid_det[1,,4,,],along=1)
y_mo <- abind(y_mid_occ[1,,1,,],y_mid_occ[1,,2,,],y_mid_occ[1,,3,,],y_mid_occ[1,,4,,],along=1)
y_ld <- abind(y_low_det[1,,1,,],y_low_det[1,,2,,],y_low_det[1,,3,,],y_low_det[1,,4,,],along=1)
y_lo <- abind(y_low_occ[1,,1,,],y_low_occ[1,,2,,],y_low_occ[1,,3,,],y_low_occ[1,,4,,],along=1)

## reduce to N or not N
N_hd <- array(NA, dim=c(48,50))
N_ho <- array(NA, dim=c(48,50))
N_md <- array(NA, dim=c(48,50))
N_mo <- array(NA, dim=c(48,50))
N_ld <- array(NA, dim=c(48,50))
N_lo <- array(NA, dim=c(48,50))

for(t in 1:50){
  for(i in 1:48){
    N_hd[i,t] <- ifelse(is.element(2, y_hd[i,,t])==TRUE | is.element(4, y_hd[i,,t])==TRUE, 1,0)
    N_ho[i,t] <- ifelse(is.element(2, y_ho[i,,t])==TRUE | is.element(4, y_ho[i,,t])==TRUE, 1,0) 
    N_md[i,t] <- ifelse(is.element(2, y_md[i,,t])==TRUE | is.element(4, y_md[i,,t])==TRUE, 1,0) 
    N_mo[i,t] <- ifelse(is.element(2, y_mo[i,,t])==TRUE | is.element(4, y_mo[i,,t])==TRUE, 1,0) 
    N_ld[i,t] <- ifelse(is.element(2, y_ld[i,,t])==TRUE | is.element(4, y_ld[i,,t])==TRUE, 1,0) 
    N_lo[i,t] <- ifelse(is.element(2, y_lo[i,,t])==TRUE | is.element(4, y_lo[i,,t])==TRUE, 1,0) 
  }}

## reduce states to P or not P
P_hd <- array(NA, dim=c(48,50)) # make empty arrays with correct dimentions
P_ho <- array(NA, dim=c(48,50))
P_md <- array(NA, dim=c(48,50))
P_mo <- array(NA, dim=c(48,50))
P_ld <- array(NA, dim=c(48,50))
P_lo <- array(NA, dim=c(48,50))

for(t in 1:50){
  for(i in 1:48){
    P_hd[i,t] <- ifelse(is.element(3, y_hd[i,,t])==TRUE | is.element(4, y_hd[i,,t])==TRUE, 1,0)
    P_ho[i,t] <- ifelse(is.element(3, y_ho[i,,t])==TRUE | is.element(4, y_ho[i,,t])==TRUE, 1,0)
    P_md[i,t] <- ifelse(is.element(3, y_md[i,,t])==TRUE | is.element(4, y_md[i,,t])==TRUE, 1,0)
    P_mo[i,t] <- ifelse(is.element(3, y_mo[i,,t])==TRUE | is.element(4, y_mo[i,,t])==TRUE, 1,0)
    P_ld[i,t] <- ifelse(is.element(3, y_ld[i,,t])==TRUE | is.element(4, y_ld[i,,t])==TRUE, 1,0)
    P_lo[i,t] <- ifelse(is.element(3, y_lo[i,,t])==TRUE | is.element(4, y_lo[i,,t])==TRUE, 1,0)
  }}

## make mean occupancy
mean.n_hd <- apply(N_hd,c(2),mean) 
mean.p_hd <- apply(P_hd,c(2),mean)

mean.n_ho <- apply(N_ho,c(2),mean) 
mean.p_ho <- apply(P_ho,c(2),mean)

mean.n_md <- apply(N_md,c(2),mean) 
mean.p_md <- apply(P_md,c(2),mean)

mean.n_mo <- apply(N_mo,c(2),mean) 
mean.p_mo <- apply(P_mo,c(2),mean)

mean.n_ld <- apply(N_ld,c(2),mean) 
mean.p_ld <- apply(P_ld,c(2),mean)

mean.n_lo <- apply(N_lo,c(2),mean) 
mean.p_lo <- apply(P_lo,c(2),mean)

# plot observation per week
setwd("../plot")

png("ObsTrends_all_sim.png", width=1640, height=480)
par(mfrow=c(2,3))
par(mar = c(5,7,4,2) + 0.1)

plot(mean.n_hd, type="l", col="blue", ylim=c(0,1), xlab="Primary occasion", ylab="Observed occupancy", main="HD", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_hd, type="l", col="red", lwd=2)
legend("topright", legend=c("Species A", "Species B"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

plot(mean.n_md, type="l", col="blue", ylim=c(0,1), xlab="Primary occasion", ylab="Observed occupancy", main="MD", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_md, type="l", col="red", lwd=2)
legend("topright", legend=c("Species A", "Species B"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

plot(mean.n_ld, type="l", col="blue", ylim=c(0,1), xlab="Primary occasion", ylab="Observed occupancy", main="LD", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_ld, type="l", col="red", lwd=2)
legend("topright", legend=c("Species A", "Species B"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

plot(mean.n_ho, type="l", col="blue", ylim=c(0,1), xlab="Primary occasion", ylab="Observed occupancy", main="HO", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_ho, type="l", col="red", lwd=2)
legend("topright", legend=c("Species A", "Species B"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

plot(mean.n_mo, type="l", col="blue", ylim=c(0,1), xlab="Primary occasion", ylab="Observed occupancy", main="MO", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_mo, type="l", col="red", lwd=2)
legend("topright", legend=c("Species A", "Species B"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

plot(mean.n_lo, type="l", col="blue", ylim=c(0,1), xlab="Primary occasion", ylab="Observed occupancy", main="LO", cex.axis=2,
     cex.main=2, cex.lab=2, lwd=2)
lines(mean.p_lo, type="l", col="red", lwd=2)
legend("topright", legend=c("Species A", "Species B"), lty=1,col=c("blue","red"), cex=1.5, lwd=2)

dev.off()