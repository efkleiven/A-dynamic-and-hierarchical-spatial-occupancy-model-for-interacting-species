####################################################################################
## Simulation study                                                               ##
## A dynamic and hiearchical spatial occupancy model for interacting species      ##
## In this script we simulate data under the high detection senario               ##
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

btpm <- array(NA, dim=c(4,4))           # block transition probability matrix
stpm <- array(NA, dim=c(ndat, B,T,4,4)) # transition probability matrix
rdm <- array(NA, dim=c(4,4))            # detection probability matrix

psi <- array(NA,dim =c(4,T))                # site occupancy probability
z_hig_det <- array(dim = c(ndat, M, B, T))  # latent site state
x_hig_det <- array(dim = c(ndat, B, T))     # latent block state

y_hig_det<- array(NA, dim = c(ndat, M, B, J, T)) # Detection histories

psi[,1] <- c(0.5,0.25,0.15,0.1) # Initial site occupancy probability

# detection probabilities
pA <- 0.9
pB <- 0.8

########################
# site level parameters
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
z_hig_det[d,, b,1] <- rcat(M, psi[,1]) # Initial site state

# Initial block state, which is just a function of the states of the sites within each block
x_hig_det[d,b,1] <- ifelse(sum(z_hig_det[d,,b,1]==1) == M, 1, 
                      ifelse(sum(z_hig_det[d,,b,1]==2) + sum(z_hig_det[d,,b,1]==1) == M, 2,
                        ifelse(sum(z_hig_det[d,,b,1]==3) + sum(z_hig_det[d,,b,1]==1) == M, 3, 4) ) )
    } #end block loop
  }   #end simulation loop 

######
# Latent state for dynamic part of model
# btpm = Block level transition probability matrix.
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
x_hig_det[d,b,t+1] <- rcat(1, btpm[ ,x_hig_det[d,b,t]])

######################################################################
## stpm = site transition probability matrix.                       ##
######################################################################

  # U to ...
  stpm[d, b, t, 1, 1] <- (1-gamA*((x_hig_det[d, b, t+1] == 2)+(x_hig_det[d, b, t+1] == 4))) * (1-gamB*((x_hig_det[d, b, t+1] == 3)+(x_hig_det[d, b, t+1] == 4)))   #--|U
  stpm[d, b, t, 2, 1] <- gamA *((x_hig_det[d, b, t+1] == 2)+(x_hig_det[d, b, t+1] == 4)) * (1-gamB*((x_hig_det[d, b, t+1] == 3)+(x_hig_det[d, b, t+1] == 4)))      #--|A
  stpm[d, b, t, 3, 1] <- (1-gamA*((x_hig_det[d, b, t+1] == 2)+(x_hig_det[d, b, t+1] == 4)) ) * gamB *((x_hig_det[d, b, t+1] == 3)+(x_hig_det[d, b, t+1] == 4))     #--|B
  stpm[d, b, t, 4, 1] <- gamA * gamB *(x_hig_det[d, b, t+1] == 4)                                                                                                  #--|AB
                       
                       
   # A to ...
  stpm[d, b, t, 1, 2] <- epsA * (1-gamBA*((x_hig_det[d, b, t+1] == 3)+(x_hig_det[d, b, t+1] == 4)))       #--|U
  stpm[d, b, t, 2, 2] <- (1-epsA) * (1-gamBA*((x_hig_det[d, b, t+1] == 3)+(x_hig_det[d, b, t+1] == 4)))   #--|A
  stpm[d, b, t, 3, 2] <- epsA * gamBA *((x_hig_det[d, b, t+1] == 3)+(x_hig_det[d, b, t+1] == 4))          #--|B
  stpm[d, b, t, 4, 2] <- (1-epsA) * gamBA  *((x_hig_det[d, b, t+1] == 3)+(x_hig_det[d, b, t+1] == 4))     #--|AB
                       
  # B to ...
  stpm[d, b, t, 1, 3] <- (1-gamAB *((x_hig_det[d, b, t+1] == 2)+(x_hig_det[d, b, t+1] == 4))) * epsB       #--|U
  stpm[d, b, t, 2, 3] <- gamAB *((x_hig_det[d, b, t+1] == 2)+(x_hig_det[d, b, t+1] == 4)) * epsB           #--|A
  stpm[d, b, t, 3, 3] <- (1-gamAB *((x_hig_det[d, b, t+1] == 2)+(x_hig_det[d, b, t+1] == 4))) * (1-epsB)   #--|B
  stpm[d, b, t, 4, 3] <- gamAB *((x_hig_det[d, b, t+1] == 2)+(x_hig_det[d, b, t+1] == 4)) * (1-epsB)       #--|AB
                       
  # AB to ..
  stpm[d, b, t, 1, 4] <- epsAB * epsBA          #--|U
  stpm[d, b, t, 2, 4] <- (1-epsAB) * epsBA      #--|A
  stpm[d, b, t, 3, 4] <- epsAB * (1-epsBA)      #--|B
  stpm[d, b, t, 4, 4] <- (1-epsAB) * (1-epsBA)  #--|AB

# Generate latent states of occurrence
# Later years

for(i in 1:M){ # Loop over sites
    z_hig_det[d, i,b,t+1] <- rcat(1, stpm[d, b, t, ,z_hig_det[d,i,b,t]])
 }         # end site loop
    }      # end block loop
      }    # end time loop
        }  # end simulation loop

# Generate detection/non-detection data

############################################################
# detection matrix (OS = observed state, TS = true state)  #
# rdm = rho detection matrix. Each row sums to 1.          #
############################################################

# TS = U
rdm[1, 1] <- 1    #----------| OS = U
rdm[2, 1] <- 0    #----------| OS = A
rdm[3, 1] <- 0    #----------| OS = B
rdm[4, 1] <- 0    #----------| OS = AB

# TS = A
rdm[1, 2] <- 1-pA #----------| OS = U
rdm[2, 2] <- pA   #----------| OS = A
rdm[3, 2] <- 0    #----------| OS = B
rdm[4, 2] <- 0    #----------| OS = AB

# TS = B
rdm[1, 3] <- 1-pB #----------| OS = U
rdm[2, 3] <- 0    #----------| OS = A
rdm[3, 3] <- pB   #----------| OS = B
rdm[4, 3] <- 0    #----------| OS = AB

# TS = AB
rdm[1, 4] <- (1-pA) * (1-pB) # OS = U
rdm[2, 4] <- pA * (1-pB)     # OS = A
rdm[3, 4] <- (1-pA) * pB     # OS = B
rdm[4, 4] <- pA * pB         # OS = AB

for(d in 1:ndat){
  for(b in 1:B){
    for(j in 1:M){
      for(k in 1:T){
        y_hig_det[d,j, b, 1:7,k] <- rcat(7, rdm[,z_hig_det[d,j,b,k]])
        }
      }
    }
  }

# save data
setwd("..") # set directory where the simulations should be saved

save(y_hig_det, file="simdata_50set_50seas_y_hig_det.rda")
save(z_hig_det, file="simdata_50set_50seas_z_hig_det.rda")
save(x_hig_det, file="simdata_50set_50seas_x_hig_det.rda")

##############################################
## plot occupancy and observations per week ##
##############################################

# plot occupancy per week

# merge blocks
z <- abind(z_hig_det[1,,1,],z_hig_det[1,,2,],z_hig_det[1,,3,],z_hig_det[1,,4,],along=1)

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

## calculate proportion of sites occupied ##
mean.z <- apply(N,c(2),mean) 
mean.p <- apply(P,c(2),mean) 

# plot trend in occupancy

setwd("../plot") # set directory where plots should be stored

png("OccTrends_mid_det.png", width=960, height=480)  #  save plot
plot(mean.z, type="l", col="blue", ylim=c(0,0.7), xlab="seasons", main="Trend in occupancy")
lines(mean.p, type="l", col="red")
legend("topright", legend=c("Prey", "Predator"), lty=1,col=c("blue","red"))
dev.off()

##############################
# plot trend in observations #
##############################

# merge blocks
y <- abind(y_hig_det[1,,1,,],y_hig_det[1,,2,,],y_hig_det[1,,3,,],y_hig_det[1,,4,,],along=1)

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

## calculate mean obervations per week
mean.n <- apply(N,c(2),mean) 
mean.p <- apply(P,c(2),mean) 

# plot observation per week
setwd("../plot")
png("ObsTrends_mid_det.png", width=960, height=480)
plot(mean.n, type="l", col="blue", ylim=c(0,0.7), xlab="seasons", main="Trend in detection")
lines(mean.p, type="l", col="red")
legend("topright", legend=c("Prey", "Predator"), lty=1,col=c("blue","red"))
dev.off()
