####################################################################################
## Simulation study                                                               ##
## A dynamic and hiearchical spatial occupancy model for interacting species      ##
## In this script we simulate data under the mid detection scenario               ##
##                                                                                ##
## by Eivind Flittie Kleiven and Frederic Barraquand                              ##
##                                                                                ##  
####################################################################################

# load libraries
library(extraDistr)
library("abind")

M <- 12    # Number of sites
B <- 4     # Number of blocks
J <- 7     # num secondary sample periods
T <- 50    # num primary sample periods
ndat <- 50 # num simulated data sets

btpm <- array(NA, dim=c(4,4)) # block transition probability matrix
stpm <- array(NA, dim=c(ndat, B,T,4,4)) # transition probability matrix
rdm <- array(NA, dim=c(4,4))  # detection probability matrix

psi <- array(NA,dim =c(4,T))               # site occupancy probability
z_mid_det <- array(dim = c(ndat, M, B, T)) # latent site state
x_mid_det <- array(dim = c(ndat, B, T))    # latent block state

y_mid_det<- array(NA, dim = c(ndat, M, B, J, T)) # Detection histories

psi[,1] <- c(0.5,0.25,0.15,0.1) # Initial site occupancy probability

# detection probabilities
pA <- 0.5
pB <- 0.5

########################
# site level paramters
########################

# colonization probability
gamA <- 0.5 
gamB <- 0.3

# conditional colonization probability
gamAB <- 0.1
gamBA <- 0.7

# extinction probability
epsA  <- 0.3
epsB  <- 0.8

# conditional extinction probability
epsAB <- 0.9
epsBA <- 0.1

########################
# block level paramters
########################

# colonization probability
GamA <- 0.5 
GamB <- 0.1

# conditional colonization probability
GamAB <- 0.05
GamBA <- 0.4

# extinction probability
EpsA  <- 0.05
EpsB  <- 0.6

# conditional extinction probability
EpsAB <- 0.5
EpsBA <- 0.2

# First year
for(d in 1:ndat){
  for(b in 1:B){
z_mid_det[d,, b,1] <- rcat(M, psi[,1]) # Initial site state

# Initial block state, which is just a function of the states of the sites within each block
x_mid_det[d,b,1] <- ifelse(sum(z_mid_det[d,,b,1]==1) == M, 1, 
                 ifelse(sum(z_mid_det[d,,b,1]==2) + sum(z_mid_det[d,,b,1]==1) == M, 2,
                        ifelse(sum(z_mid_det[d,,b,1]==3) + sum(z_mid_det[d,,b,1]==1) == M, 3, 4) ) )
}}

######
# Latent state for dynamic part of model
# btpm = block transition probability matrix. All columns sum to 1.
######

# U to ...
btpm[1, 1] <- (1-GamA) * (1-GamB)    #--|U
btpm[2, 1] <- GamA * (1-GamB)        #--|A
btpm[3, 1] <- (1-GamA) * GamB        #--|B
btpm[4, 1] <- GamA * GamB            #--|AB

# A to ...
btpm[1, 2] <- EpsA * (1-GamBA)       #--|U
btpm[2, 2] <- (1-EpsA) * (1-GamBA)   #--|A
btpm[3, 2] <- EpsA * GamBA           #--|B
btpm[4, 2] <- (1-EpsA) * GamBA       #--|AB

# B to ...
btpm[1, 3] <- (1-GamAB) * EpsB       #--|U
btpm[2, 3] <- GamAB * EpsB           #--|A
btpm[3, 3] <- (1-GamAB) * (1-EpsB)   #--|B
btpm[4, 3] <- GamAB * (1-EpsB)       #--|AB

# AB to ..
btpm[1, 4] <- EpsAB * EpsBA          #--|U
btpm[2, 4] <- (1-EpsAB) * EpsBA      #--|A
btpm[3, 4] <- EpsAB * (1-EpsBA)      #--|B
btpm[4, 4] <- (1-EpsAB) * (1-EpsBA)  #--|AB

for(d in 1:ndat){
  for(t in 1:(T-1)){
    for(b in 1:B){
# generate latent block state
x_mid_det[d,b,t+1] <- rcat(1, btpm[ ,x_mid_det[d,b,t]])

######################################################################
## stpm = site transition probability matrix. All columns sum to 1. ##
######################################################################

  # U to ...
  stpm[d, , t, 1, 1]  <- (1-gamA*((x_mid_det[d, b, t+1] == 2)+(x_mid_det[d, b, t+1] == 4))) * (1-gamB*((x_mid_det[d, b, t+1] == 3)+(x_mid_det[d, b, t+1] == 4)))  #--|U
  stpm[d, b, t, 2, 1] <- gamA *((x_mid_det[d, b, t+1] == 2)+(x_mid_det[d, b, t+1] == 4)) * (1-gamB*((x_mid_det[d, b, t+1] == 3)+(x_mid_det[d, b, t+1] == 4)))     #--|A
  stpm[d, b, t, 3, 1] <- (1-gamA*((x_mid_det[d, b, t+1] == 2)+(x_mid_det[d, b, t+1] == 4)) ) * gamB *((x_mid_det[d, b, t+1] == 3)+(x_mid_det[d, b, t+1] == 4))    #--|B
  stpm[d, b, t, 4, 1] <- gamA * gamB *(x_mid_det[d, b, t+1] == 4)                                                                                                 #--|AB
                       
                       
   # A to ...
  stpm[d, b, t, 1, 2] <- epsA * (1-gamBA*((x_mid_det[d, b, t+1] == 3)+(x_mid_det[d, b, t+1] == 4)))       #--|U
  stpm[d, b, t, 2, 2] <- (1-epsA) * (1-gamBA*((x_mid_det[d, b, t+1] == 3)+(x_mid_det[d, b, t+1] == 4)))   #--|A
  stpm[d, b, t, 3, 2] <- epsA * gamBA *((x_mid_det[d, b, t+1] == 3)+(x_mid_det[d, b, t+1] == 4))          #--|B
  stpm[d, b, t, 4, 2] <- (1-epsA) * gamBA  *((x_mid_det[d, b, t+1] == 3)+(x_mid_det[d, b, t+1] == 4))     #--|AB
                       
  # B to ...
  stpm[d, b, t, 1, 3] <- (1-gamAB *((x_mid_det[d, b, t+1] == 2)+(x_mid_det[d, b, t+1] == 4))) * epsB       #--|U
  stpm[d, b, t, 2, 3] <- gamAB *((x_mid_det[d, b, t+1] == 2)+(x_mid_det[d, b, t+1] == 4)) * epsB           #--|A
  stpm[d, b, t, 3, 3] <- (1-gamAB *((x_mid_det[d, b, t+1] == 2)+(x_mid_det[d, b, t+1] == 4))) * (1-epsB)   #--|B
  stpm[d, b, t, 4, 3] <- gamAB *((x_mid_det[d, b, t+1] == 2)+(x_mid_det[d, b, t+1] == 4)) * (1-epsB)       #--|AB
                       
  # AB to ..
  stpm[d, b, t, 1, 4] <- epsAB * epsBA          #--|U
  stpm[d, b, t, 2, 4] <- (1-epsAB) * epsBA      #--|A
  stpm[d, b, t, 3, 4] <- epsAB * (1-epsBA)      #--|B
  stpm[d, b, t, 4, 4] <- (1-epsAB) * (1-epsBA)  #--|AB

# Generate latent states of occurrence
# Later years

for(i in 1:M){ # Loop over sites
    z_mid_det[d, i,b,t+1] <- rcat(1, stpm[d, b, t, ,z_mid_det[d,i,b,t]])
    }}}}


# Generate detection/non-detection data

######
# detection matrix (OS = observed state, TS = true state)
# rdm = rho detection matrix. Each row sums to 1.
# OS along rows, TS along columns
######
# TS = U
rdm[1, 1] <- 1    #------------|OS = U
rdm[2, 1] <- 0    #------------|OS = A
rdm[3, 1] <- 0    #------------|OS = B
rdm[4, 1] <- 0    #------------|OS = AB

# TS = A
rdm[1, 2] <- 1-pA #------------|OS = U
rdm[2, 2] <- pA   #------------|OS = A
rdm[3, 2] <- 0    #------------|OS = B
rdm[4, 2] <- 0    #------------|OS = AB

# TS = B
rdm[1, 3] <- 1-pB #------------|OS = U
rdm[2, 3] <- 0    #------------|OS = A
rdm[3, 3] <- pB   #------------|OS = B
rdm[4, 3] <- 0    #------------|OS = AB

# TS = AB
rdm[1, 4] <- (1-pA) * (1-pB) #-|OS = U
rdm[2, 4] <- pA * (1-pB)     #-|OS = A
rdm[3, 4] <- (1-pA) * pB     #-|OS = B
rdm[4, 4] <- pA * pB         #-|OS = AB

for(d in 1:ndat){
  for(b in 1:B){
    for(j in 1:M){
      for(k in 1:T){
        y_mid_det[d,j, b, 1:7,k] <- rcat(7, rdm[,z_mid_det[d,j,b,k]])
        }
      }
    }
  }

# save data
setwd("H:/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models//hidden_block_sim/data")

save(y_mid_det, file="simdata_50set_50seas_y_mid_det.rda")
save(z_mid_det, file="simdata_50set_50seas_z_mid_det.rda")
save(x_mid_det, file="simdata_50set_50seas_x_mid_det.rda")


###############################
## plotting

# remove block structure to make plot on site level
z <- abind(z_mid_det[1,,1,],z_mid_det[1,,2,],z_mid_det[1,,3,],z_mid_det[1,,4,],along=1)

## reduce to N or not N
N <- array(NA, dim=c(48,50))
for(t in 1:50){
    for(i in 1:48){
      N[i,t] <- ifelse(is.element(2, z[i,t])==TRUE | is.element(4, z[i,t])==TRUE, 1,0) 
    }}

## reduce to P or not P
P <- array(NA, dim=c(48,50))
for(t in 1:50){
  for(i in 1:48){
    P[i,t] <- ifelse(is.element(3, z[i,t])==TRUE | is.element(4, z[i,t])==TRUE, 1,0) 
  }}

## proportion of sites occupied ##
mean.z <- apply(N,c(2),mean) 
mean.p <- apply(P,c(2),mean) 

# plot trend in occupancy
setwd("../plot")
png("OccTrends_mid_det.png", width=960, height=480)
plot(mean.z, type="l", col="blue", ylim=c(0,0.7), xlab="seasons", main="Trend in occupancy")
lines(mean.p, type="l", col="red")
legend("topright", legend=c("Prey", "Predator"), lty=1,col=c("blue","red"))
dev.off()

# plot trend in detections

# remove block structure to make plot on site level
y <- abind(y[,1,,],y[,2,,],y[,3,,],y[,4,,],along=1)

## reduce to N or not N
N <- array(NA, dim=c(48,50))

for(t in 1:50){
  for(i in 1:48){
    N[i,t] <- ifelse(is.element(2, y[i,,t])==TRUE | is.element(4, y[i,,t])==TRUE, 1,0) 
  }}

P <- array(NA, dim=c(48,50))
for(t in 1:50){
  for(i in 1:48){
    P[i,t] <- ifelse(is.element(3, y[i,,t])==TRUE | is.element(4, y[i,,t])==TRUE, 1,0) 
  }}

## make mean occupancy
mean.n <- apply(N,c(2),mean) 
mean.p <- apply(P,c(2),mean) 

# plot observation per week
setwd("../plot")
png("ObsTrends_mid_det.png", width=960, height=480)
plot(mean.n, type="l", col="blue", ylim=c(0,0.7), xlab="seasons", main="Trend in occupancy")
lines(mean.p, type="l", col="red")
legend("topright", legend=c("Prey", "Predator"), lty=1,col=c("blue","red"))
dev.off()
