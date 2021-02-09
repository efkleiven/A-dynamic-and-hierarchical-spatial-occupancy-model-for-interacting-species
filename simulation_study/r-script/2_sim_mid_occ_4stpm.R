##########################################################################################################################################
##  A dynamic and hierarchical spatial occupancy model for interacting species analyzing simulated data under a mid occupancy scenario  ##
##                                  by Eivind Flittie Kleiven and Frederic Barraquand                                                   ##
##########################################################################################################################################

# Call jags(and other packages)
rm(list=ls())

library(jagsUI)

# set working directory
setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel")
setwd("./models/hidden_block_sim")

## import data
setwd("./data")

# data simulated under this model
load("simdata_50set_50seas_y_mid_occ.rda")     # data simulated under this model
#load("SimData.RData")

# Parameters monitored
params <- c("gamA","gamB","gamAB","gamBA","epsA","epsB","epsAB","epsBA","psi",
            "GamA","GamB","GamAB","GamBA","EpsA","EpsB","EpsAB","EpsBA", "pA","pB",
            "z","x")

# MCMC settings
ni <- 5000  ;   nt <- 10   ;   nb <- 0   ;   nc <- 4    ;   na <- 1000

###################################
# make loop over 50 data sets
###################################
mod_4stpm_sim_mid_occ <- list() # make empty list to store models in
setwd("../") # set wd before loop

for(q in 1:50){
yb <- y_mid_occ[q, , , , ] 
yb <- aperm(yb, c(1,2,4,3))

# give data
data <-list(nseason = dim(yb)[3], nblock = dim(yb)[2], nsite = dim(yb)[1], nsurvey = dim(yb)[4], 
            nout=4, y = yb) 
nseason = dim(yb)[3]; nblock = dim(yb)[2]; nsite = dim(yb)[1]; nsurvey = dim(yb)[4]

# Initial values 
sp_inits <- apply(yb,c(1,2,3),max)

# give initial values
inits=function(){list( 
  z = sp_inits, pA=runif(1), pB=runif(1), gamB=runif(1), gamAB=runif(1), gamBA=runif(1),
  epsA=runif(1), epsB=runif(1), epsAB=runif(1), epsBA=runif(1), GamA=runif(1), 
  GamB=runif(1), GamAB=runif(1), GamBA=runif(1),EpsA=runif(1), 
  EpsB=runif(1), EpsAB=runif(1), EpsBA=runif(1)
)}

# loop to make cases where both state 2 and 3 is observed within the same sampling occasion have initial value 4
dataL <- array(0,dim=c(nsite,nblock,nseason))
for(j in 1:nsite){
  for(b in 1:nblock){
    for(i in 1:nseason){
      if (is.element(0, match(c(2,3),yb[j,b,i,], nomatch = FALSE, incomparables = FALSE))) 
        dataL[j,b,i] <- "FALSE"
      else
        dataL[j,b,i] <- "TRUE"
    }}}

for(j in 1:nsite){
  for(b in 1:nblock){
    for(i in 1:nseason){
      if(dataL[j,b,i]==TRUE){
        sp_inits[j,b,i] <- 4}
    }}}


# run model in jags

mod_4stpm_sim_mid_occ[[q]] <- jags(data, inits=inits, params, "mod.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)
}

# Save model
setwd("./model_output")
save(mod_4stpm_sim_mid_occ, file="mod_4stpm_sim_mid_occ.rda")

