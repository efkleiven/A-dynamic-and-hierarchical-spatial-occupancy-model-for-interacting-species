####################################################################################
## Simulation study                                                               ##
## A dynamic occupancy model for interacting species with two spatial scales      ##
## In this script we simulate data under the mid detection scenario               ##
##                                                                                ##
## by Eivind Flittie Kleiven and Frederic Barraquand                              ##
##                                                                                ##  
####################################################################################

# load libraries
library(extraDistr)
library("abind")

M <- 12    # Number of sites
B <- 8     # Number of blocks
J <- 7     # num secondary sample periods
T <- 50    # num primary sample periods
ndat <- 50 # num simulated data sets

btpm <- array(NA, dim=c(4,4)) # block transition probability matrix
stpm <- array(NA, dim=c(4,4,4)) # transition probability matrix
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
z_mid_det[d,, b,1] <- rcat(M, psi[,1]) # Initial site state

# Initial block state, which is just a function of the states of the sites within each block
x_mid_det[d,b,1] <- ifelse(sum(z_mid_det[d,,b,1]==1) == M, 1, 
                 ifelse(sum(z_mid_det[d,,b,1]==2) + sum(z_mid_det[d,,b,1]==1) == M, 2,
                        ifelse(sum(z_mid_det[d,,b,1]==3) + sum(z_mid_det[d,,b,1]==1) == M, 3, 4) ) )
}}

######
# Latent state for dynamic part of model
# btpm = block transition probability matrix. 
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


####################################################
## stpm = site transition probability matrix.     ##
## These are dependent on the block level state   ##
## All columns sum to 1.                          ##
####################################################

# blocks state (x) = U

# U to ...
stpm[ 1, 1, 1] <- 1          #--|U
stpm[ 2, 1, 1] <- 0          #--|A
stpm[ 3, 1, 1] <- 0          #--|B
stpm[ 4, 1, 1] <- 0          #--|AB

# A to ...
stpm[ 1, 2, 1] <- 1          #--|U
stpm[ 2, 2, 1] <- 0          #--|A
stpm[ 3, 2, 1] <- 0          #--|B
stpm[ 4, 2, 1] <- 0          #--|AB

# B to ...
stpm[1, 3, 1] <- 1           #--|U
stpm[2, 3, 1] <- 0           #--|A
stpm[3, 3, 1] <- 0           #--|B
stpm[4, 3, 1] <- 0           #--|AB

# AB to ..
stpm[1, 4, 1] <- 1           #--|U
stpm[2, 4, 1] <- 0           #--|A
stpm[3, 4, 1] <- 0           #--|B
stpm[4, 4, 1] <- 0           #--|AB

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# blocks state (x) = A

# U to ...
stpm[1, 1, 2] <- (1-gamA)    #--|U
stpm[2, 1, 2] <- gamA        #--|A
stpm[3, 1, 2] <- 0           #--|B
stpm[4, 1, 2] <- 0           #--|AB

# A to ...
stpm[1, 2, 2] <- epsA        #--|U
stpm[2, 2, 2] <- (1-epsA)    #--|A
stpm[3, 2, 2] <- 0           #--|B
stpm[4, 2, 2] <- 0           #--|AB

# B to ...
stpm[1, 3, 2] <- (1-gamAB)   #--|U
stpm[2, 3, 2] <- gamAB       #--|A
stpm[3, 3, 2] <- 0           #--|B
stpm[4, 3, 2] <- 0           #--|AB

# AB to ..
stpm[1, 4, 2] <- epsAB       #--|U
stpm[2, 4, 2] <- (1-epsAB)   #--|A
stpm[3, 4, 2] <- 0           #--|B
stpm[4, 4, 2] <- 0           #--|AB

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# blocks state (x) = B

# U to ...
stpm[1, 1, 3] <-  (1-gamB)    #--|U
stpm[2, 1, 3] <- 0            #--|A
stpm[3, 1, 3] <-  gamB        #--|B
stpm[4, 1, 3] <- 0            #--|AB

# A to ...
stpm[1, 2, 3] <- (1-gamBA)    #--|U
stpm[2, 2, 3] <- 0            #--|A
stpm[3, 2, 3] <- gamBA        #--|B
stpm[4, 2, 3] <- 0            #--|AB

# B to ...
stpm[1, 3, 3] <-  epsB        #--|U
stpm[2, 3, 3] <- 0            #--|A
stpm[3, 3, 3] <- (1-epsB)     #--|B
stpm[4, 3, 3] <- 0            #--|AB

# AB to ..
stpm[1, 4, 3] <-  epsBA       #--|U
stpm[2, 4, 3] <- 0            #--|A
stpm[3, 4, 3] <-  (1-epsBA)   #--|B
stpm[4, 4, 3] <- 0            #--|AB


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# blocks state (x) = AB

# U to ...
stpm[1, 1, 4] <- (1-gamA) * (1-gamB)    #--|U
stpm[2, 1, 4] <- gamA  * (1-gamB)       #--|A
stpm[3, 1, 4] <- (1-gamA) * gamB        #--|B
stpm[4, 1, 4] <- gamA * gamB            #--|AB

# A to ...
stpm[1, 2, 4] <- epsA * (1-gamBA)       #--|U
stpm[2, 2, 4] <- (1-epsA) * (1-gamBA)   #--|A
stpm[3, 2, 4] <- epsA * gamBA           #--|B
stpm[4, 2, 4] <- (1-epsA) * gamBA       #--|AB

# B to ...
stpm[1, 3, 4] <- (1-gamAB ) * epsB      #--|U
stpm[2, 3, 4] <- gamAB  * epsB          #--|A
stpm[3, 3, 4] <- (1-gamAB ) * (1-epsB)  #--|B
stpm[4, 3, 4] <- gamAB  * (1-epsB)      #--|AB

# AB to ..
stpm[1, 4, 4] <- epsAB * epsBA          #--|U
stpm[2, 4, 4] <- (1-epsAB) * epsBA      #--|A
stpm[3, 4, 4] <- epsAB * (1-epsBA)      #--|B
stpm[4, 4, 4] <- (1-epsAB) * (1-epsBA)  #--|AB
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate latent states of occurrence
# Later years
for(d in 1:ndat){
  for(t in 1:(T-1)){
    for(b in 1:B){
      # generate latent block state
      x_mid_det[d,b,t+1] <- rcat(1, btpm[ ,x_mid_det[d,b,t]])
      
      # generate latent site state      
      for(i in 1:M){ # Loop over sites
        z_mid_det[d, i,b,t+1] <- rcat(1, stpm[ ,z_mid_det[d,i,b,t], x_mid_det[d,b,t+1]])
        }}}}

# Generate detection/non-detection data

######
# detection matrix (OS = observed state, TS = true state)
# rdm = rho detection matrix. 
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
setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models//hidden_block_sim/data")

save(y_mid_det, file="simdata_50set_50seas_8b_y_mid_det.rda")
save(z_mid_det, file="simdata_50set_50seas_8b_z_mid_det.rda")
save(x_mid_det, file="simdata_50set_50seas_8b_x_mid_det.rda")
