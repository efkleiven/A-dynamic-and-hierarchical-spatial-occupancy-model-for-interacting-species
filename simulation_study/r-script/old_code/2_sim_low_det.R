##########################################################################################################################################
## A dynamic and hierarchical spatial occupancy model for interacting species analyzing simulated data under a low detection scenario   ##
##                                  by Eivind Flittie Kleiven and Frederic Barraquand                                                   ##
##########################################################################################################################################

# clear work space
rm(list=ls())

# Call jags(and other packages)
library(jagsUI)

# set working directory
setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel")
setwd("./models/hidden_block_sim")

# Specify model in BUGS language
sink("mac_dmsom_sdet.txt")
cat("
    model{
    
    ###########################################
    # State 1= Unoccupied(U), State 2= vole(A), State 3 = stoat(B), State 4 = vole & stoat(AB)  
    ########################################### 
   
    ##############
    ##  Priors  ##    
    ##############
    
    # detection parameters
    pA    ~ dunif(0,1)
    pB    ~ dunif(0,1)

    # site parameters
    gamA  ~ dunif(0,1)
    gamB  ~ dunif(0,1)
    gamAB ~ dunif(0,1)
    gamBA ~ dunif(0,1)
    epsA  ~ dunif(0,1)
    epsB  ~ dunif(0,1)
    epsAB ~ dunif(0,1)
    epsBA ~ dunif(0,1)
    
    # block parameters
    GamA  ~ dunif(0,1)
    GamB  ~ dunif(0,1)
    GamAB ~ dunif(0,1)
    GamBA ~ dunif(0,1)
    EpsA  ~ dunif(0,1)
    EpsB  ~ dunif(0,1)
    EpsAB ~ dunif(0,1)
    EpsBA ~ dunif(0,1)

  # initial state parameters
    for(b in 1:nblock){
    for(i in 1:3){
    psi[b,i] ~ dunif(0,0.5) # site
    }
  
    # First season probabilities for each state

    #site    
    for(j in 1:nsite){
    fsm[j, b, 1] <- 1-psi[b,1]-psi[b,2]-psi[b,3]  #-----------|U
    fsm[j, b, 2] <- psi[b,1]                      #-----------|A
    fsm[j, b, 3] <- psi[b,2]                      #-----------|B
    fsm[j, b, 4] <- psi[b,3]                      #-----------|AB

    
    # first season latent state
    # for sites   
      z[j, b, 1] ~ dcat( fsm[j, b, ( 1:nout )])    # site 
    } # end site loop
    
    # for block, which is just a function of the states of the sites within each block
   
     x[b,1] <- ifelse(sum(z[,b,1]==1) == nsite, 1,
                ifelse(sum(z[,b,1]==2) + sum(z[,b,1]==1) == nsite, 2,
                  ifelse(sum(z[,b,1]==3) + sum(z[,b,1]==1) == nsite, 3, 4) ) )
    
    } # end block loop

    
    #############################################
    # block level transition probability matrix #
    #############################################  
    
    ##########################################
    # btpm = block transition probability matrix. All columns sum to 1.
    #############################################
    
    # U to ...
    btpm[ 1, 1] <- (1-GamA) * (1-GamB)    #--|U
    btpm[ 2, 1] <- GamA * (1-GamB)        #--|A
    btpm[ 3, 1] <- (1-GamA) * GamB        #--|B
    btpm[ 4, 1] <- GamA * GamB            #--|AB
    
    # A to ...
    btpm[ 1, 2] <- EpsA * (1-GamBA)       #--|U
    btpm[ 2, 2] <- (1-EpsA) * (1-GamBA)   #--|A
    btpm[ 3, 2] <- EpsA * GamBA           #--|B
    btpm[ 4, 2] <- (1-EpsA) * GamBA       #--|AB
    
    # B to ...
    btpm[ 1, 3] <- (1-GamAB) * EpsB       #--|U
    btpm[ 2, 3] <- GamAB * EpsB           #--|A
    btpm[ 3, 3] <- (1-GamAB) * (1-EpsB)   #--|B
    btpm[ 4, 3] <- GamAB * (1-EpsB)       #--|AB
    
    # AB to ..
    btpm[ 1, 4] <- EpsAB * EpsBA           #--|U
    btpm[ 2, 4] <- (1-EpsAB) * EpsBA      #--|A
    btpm[ 3, 4] <- EpsAB * (1-EpsBA)      #--|B
    btpm[ 4, 4] <- (1-EpsAB) * (1-EpsBA)  #--|AB

    # latent block state for the rest of the seasons
    for(t in 1:(nseason-1)){
      for(b in 1:nblock){
        x[b, t+1] ~ dcat(btpm[(1:nout), x[b, t]])
      
    ######################################################################
    ## stpm = site transition probability matrix. All columns sum to 1. ##
    ######################################################################
    
    # U to ...
    stpm[b, t, 1, 1] <- (1-gamA*((x[b, t+1] == 2)+(x[b, t+1] == 4))) * (1-gamB*((x[b, t+1] == 3)+(x[b, t+1] == 4)))    #--|U
    stpm[b, t, 2, 1] <- gamA *((x[b, t+1] == 2)+(x[b, t+1] == 4)) * (1-gamB*((x[b, t+1] == 3)+(x[b, t+1] == 4)))       #--|A
    stpm[b, t, 3, 1] <- (1-gamA*((x[b, t+1] == 2)+(x[b, t+1] == 4))) * gamB *((x[b, t+1] == 3)+(x[b, t+1] == 4))       #--|B
    stpm[b, t, 4, 1] <- gamA * gamB *(x[b, t+1] == 4)        #--|AB
    
    # A to ...
    stpm[b, t, 1, 2] <- epsA * (1-gamBA*((x[b, t+1] == 3)+(x[b, t+1] == 4)))       #--|U
    stpm[b, t, 2, 2] <- (1-epsA) * (1-gamBA*((x[b, t+1] == 3)+(x[b, t+1] == 4)))   #--|A
    stpm[b, t, 3, 2] <- epsA * gamBA *((x[b, t+1] == 3)+(x[b, t+1] == 4))          #--|B
    stpm[b, t, 4, 2] <- (1-epsA) * gamBA  *((x[b, t+1] == 3)+(x[b, t+1] == 4))     #--|AB
    
    # B to ...
    stpm[b, t, 1, 3] <- (1-gamAB *((x[b, t+1] == 2)+(x[b, t+1] == 4))) * epsB       #--|U
    stpm[b, t, 2, 3] <- gamAB *((x[b, t+1] == 2)+(x[b, t+1] == 4)) * epsB           #--|A
    stpm[b, t, 3, 3] <- (1-gamAB *((x[b, t+1] == 2)+(x[b, t+1] == 4))) * (1-epsB)   #--|B
    stpm[b, t, 4, 3] <- gamAB *((x[b, t+1] == 2)+(x[b, t+1] == 4)) * (1-epsB)       #--|AB
    
    # AB to ..
    stpm[b, t, 1, 4] <- epsAB * epsBA          #--|U
    stpm[b, t, 2, 4] <- (1-epsAB) * epsBA      #--|A
    stpm[b, t, 3, 4] <- epsAB * (1-epsBA)      #--|B
    stpm[b, t, 4, 4] <- (1-epsAB) * (1-epsBA)  #--|AB

      # latent site state for the rest of the seasons
      for(j in 1:nsite){
        z[j, b, t+1] ~ dcat( stpm[b, t, ( 1:nout ) , z[ j, b, t]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
       
          for(day in 1:nsurvey) {      
            y[j, b, t, day] ~ dcat( dpm[ ( 1:nout ) , z[j, b, t]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
          } #end survey loop
        } # end site loop
      } # end block loop
    } # end time loop
    
    #############################################################
    ## detection matrix (OS = observed state, TS = true state) ##
    #############################################################
    
    # TS = U
    dpm[ 1, 1] <- 1                #--|OS = U
    dpm[ 2, 1] <- 0                #--|OS = A
    dpm[ 3, 1] <- 0                #--|OS = B
    dpm[ 4, 1] <- 0                #--|OS = AB
    
    # TS = A
    dpm[ 1, 2] <- 1-pA             #--|OS = U
    dpm[ 2, 2] <- pA               #--|OS = A
    dpm[ 3, 2] <- 0                #--|OS = B
    dpm[ 4, 2] <- 0                #--|OS = AB
    
    # TS = B
    dpm[ 1, 3] <- 1-pB             #--|OS = U
    dpm[ 2, 3] <- 0                #--|OS = A
    dpm[ 3, 3] <- pB               #--|OS = B
    dpm[ 4, 3] <- 0                #--|OS = AB
    
    # TS = AB
    dpm[ 1, 4] <- (1-pA) * (1-pB)  #--|OS = U
    dpm[ 2, 4] <- pA * (1-pB)      #--|OS = A
    dpm[ 3, 4] <- (1-pA) * pB      #--|OS = B
    dpm[ 4, 4] <- pA * pB          #--|OS = AB

    }# end
    ",fill = TRUE)
sink()

## import data
setwd("./data")

# data simulated under this model
#load("simdata_50set_50seas_y_low_det.rda")     # data simulated under this model
load("SimData.RData")
dim(y_low_det)

mod_hidden_sim_low_det <- list() # make empty list to store models in
setwd("../") # set wd before loop

# Parameters monitored
params <- c("gamA","gamB","gamAB","gamBA","epsA","epsB","epsAB","epsBA","psi",
            "GamA","GamB","GamAB","GamBA","EpsA","EpsB","EpsAB","EpsBA", "Psi", "pA","pB",
            "z","x")

# MCMC settings
ni <- 5000  ;   nt <- 10   ;   nb <- 0   ;   nc <- 4    ;   na <- 1000

#############################
# loop over 50 data sets

for(q in 1:50){
  yb <- y_low_det[q, , , , ] #
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
  
  mod_hidden_sim_low_det[[q]] <- jags(data, inits=inits, params, "mac_dmsom_sdet.txt", n.chains = nc,
                                      n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)
}

# Save model
setwd("./model_output")
save(mod_hidden_sim_low_det, file="mod_hidden_sim_low_det.rda")

