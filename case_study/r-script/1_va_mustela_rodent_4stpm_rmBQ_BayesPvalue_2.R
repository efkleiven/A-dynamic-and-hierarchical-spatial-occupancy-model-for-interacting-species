###################################################################################################################
##  A dynamic occupancy model for interacting species with two spatial scales model analyzing data from         ##  
##  a long-term monitoring program of small mammals on the arctic tundra                                         ##
##                 by Eivind Flittie Kleiven and Frederic Barraquand                                             ##
##    Last updated 27.1.21                                                                                       ##
###################################################################################################################

# a part of the code is modified from https://github.com/mikemeredith/AHM_code/blob/master/AHM2_ch04/AHM2_04.08.R

# Call jags(and other packages)
rm(list=ls())

library(jagsUI)

# set working directory
#setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel")
#setwd("./models/hidden_block_ko_fromautoclass/ko_vj")
setwd("~/UiT/GitProjects/A-dynamic-and-hierarchical-spatial-occupancy-model-for-interacting-species/case_study")

#
# Specify model in JAGS language
sink("mod_spatial_dynocc_bayp.txt")
cat("
    model{
    
    ###########################################
    ### Spatial Dynamic Co-Occurrence Model ###
    ###########################################
    # State 1= Unoccupied(U), State 2= rodent(A), State 3 = mustelid(B), State 4 = rodent & mustelid(AB)  
    ########################################### 
   
    ##############
    ##  Priors  ##    
    ##############
    
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
    
    # interscept det prob

    for(i in 1:2){    
      alphaA0[i] ~ dnorm(0,1)
      alphaB0[i] ~ dnorm(0,1)
    }

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

    ###############################################
    # btpm = block transition probability matrix. #
    # All columns sum to 1.                       #
    ###############################################
    
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
    btpm[ 1, 4] <- EpsAB * EpsBA          #--|U
    btpm[ 2, 4] <- (1-EpsAB) * EpsBA      #--|A
    btpm[ 3, 4] <- EpsAB * (1-EpsBA)      #--|B
    btpm[ 4, 4] <- (1-EpsAB) * (1-EpsBA)  #--|AB


      
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
  
    # latent block state for the rest of the seasons
    for(t in 1:(nseason-1)){
      for(b in 1:nblock){
        x[b, t+1] ~ dcat(btpm[(1:nout), x[b, t]])

      # latent site state for the rest of the seasons
      for(j in 1:nsite){
        z[j, b, t+1] ~ dcat( stpm[( 1:nout ) , z[ j, b, t], x[b,t+1]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
      } # end site loop
      } # end block loop
    } #close time loop
    
          for(i in 1:nvisits) {      
            y[i] ~ dcat( dpm[time_id, ( 1:nout ) , z[site_id[i], block_id[i], time_id[i]]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
          } #end nvisit loop


    
    #############################################################
    ## detection matrix (OS = observed state, TS = true state) ##
    #############################################################
    
  for(t in 1:nseason){
    
    # TS = U
    dpm[t, 1, 1] <- 1                      #--|OS = U
    dpm[t, 2, 1] <- 0                      #--|OS = A
    dpm[t, 3, 1] <- 0                      #--|OS = B
    dpm[t, 4, 1] <- 0                      #--|OS = AB
    
    # TS = A
    dpm[t, 1, 2] <- 1-pA[t]                #--|OS = U
    dpm[t, 2, 2] <- pA[t]                  #--|OS = A
    dpm[t, 3, 2] <- 0                      #--|OS = B
    dpm[t, 4, 2] <- 0                      #--|OS = AB
    
    # TS = B
    dpm[t, 1, 3] <- 1-pB[t]                #--|OS = U
    dpm[t, 2, 3] <- 0                      #--|OS = A
    dpm[t, 3, 3] <- pB[t]                  #--|OS = B
    dpm[t, 4, 3] <- 0                      #--|OS = AB
    
    # TS = AB
    dpm[t, 1, 4] <- (1-pA[t]) * (1-pB[t])  #--|OS = U
    dpm[t, 2, 4] <- pA[t] * (1-pB[t])      #--|OS = A
    dpm[t, 3, 4] <- (1-pA[t]) * pB[t]      #--|OS = B
    dpm[t, 4, 4] <- pA[t] * pB[t]          #--|OS = AB
    
    ## logit links for detection probs
    logit(pA[t]) <- alphaA0[season[t]]
    logit(pB[t]) <- alphaB0[season[t]]
    } #close time loop
    
    ## Derived parameters
    
    diff_gamA <- gamA - gamAB
    diff_gamB <- gamB - gamBA
    diff_epsA <- epsA - epsAB
    diff_epsB <- epsB - epsBA
    
    diff_GamA <- GamA - GamAB
    diff_GamB <- GamB - GamBA
    diff_EpsA <- EpsA - EpsAB
    diff_EpsB <- EpsB - EpsBA
    
    ratio_gamA <- gamA / gamAB
    ratio_gamB <- gamB / gamBA
    ratio_epsA <- epsAB / epsA
    ratio_epsB <- epsB / epsBA
    
    ratio_GamA <- GamA / GamAB
    ratio_GamB <- GamBA / GamB
    ratio_EpsA <- EpsAB / EpsA
    ratio_EpsB <- EpsB / EpsBA    
    
  # (2) GoF computation part of code
  # (based on posterior predictive distributions)
  # --------------------------------------------
  # Draw a replicate data set under the fitted model
   
     for(t in 1:nseason){
      for(b in 1:nblock){
        for(j in 1:nsite){
           for(day in 1:nsurvey){

              yrep[j, b, t, day] ~ dcat(dpm[t, , z[j, b, t]] + 0.01)

           } # end survey loop
          } # end site loop
        } # end block loop
      } # end survey loop
      
  # (2a) Computations for the GoF of the open part of the model
  # (based on number of state transitions)
  # ----------------------------------------------------------
  
  for(t in 1:nseason){
    for(b in 1:nblock){
      for(j in 1:nsite){
         for(day in 1:nsurvey) {  
          y2[(b-1)*6+j,t, day] <- yb[j,b,t,day]
          yrep2[(b-1)*6+j ,t,day] <- yrep[j, b, t, day]
          } # end survey loop

          z2[(b-1)*6+j, t] <- z[j, b, t]         
        } # end site loop
      } # end block loop
    } # end time loop


# seperate species
  yrep_22 <- yrep2==2       # select data point where status 2 (only A observed)
  yrep_2 <- yrep_22*1      # make numeric

  yrep_44 <- yrep2==4       # select data point where status 4 ( A and B observed)
  yrep_4 <- yrep_44*1      # make numeric

  yrep_A <- yrep_2+yrep_4 

  yrep_33 <- yrep2==3       # select data point where status 2 (only A observed)
  yrep_3 <- yrep_33*1      # make numeric

  yrep_B <- yrep_3+yrep_4 

  y_22 <- y2==2       # select data point where status 2 (only A observed)
  y_2 <- y_22*1      # make numeric

  y_44 <- y2==4       # select data point where status 4 ( A and B observed)
  y_4 <- y_44*1      # make numeric

  y_A <- y_2+y_4 

  y_33 <- y2==3       # select data point where status 2 (only A observed)
  y_3 <- y_33*1      # make numeric

 y_B <- y_3+y_4 
  
  # Compute observed z matrix for observed and replicated data
 for(t in 1:nseason){
    for(b in 1:nblock){
      for(j in 1:nsite){
      
        zobs_A[(b-1)*6+j,t] <- max(y_A[(b-1)*6+j,t,]) 
        zobsrep_A[(b-1)*6+j,t] <- max(yrep_A[(b-1)*6+j,t,]) # For replicated data
       
        zobs_B[(b-1)*6+j,t] <- max(y_B[(b-1)*6+j,t,])       # For observed data
        zobsrep_B[(b-1)*6+j,t] <- max(yrep_B[(b-1)*6+j,t,]) # For replicated data
        } # end site loop
      } # end block loop
    } # end time loop
        
    # Identify extinctions, persistence, colonization and non-colonizations
    
  for (t in 2:nseason){
    for(b in 1:nblock){
      for(j in 1:nsite){ 
      
      # ... for observed data
      ext_A[(b-1)*6+j,t-1] <- equals(zobs_A[(b-1)*6+j,t],0) * equals(zobs_A[(b-1)*6+j,t-1],1)
      nonext_A[(b-1)*6+j,(t-1)] <- equals(zobs_A[(b-1)*6+j,t],1) * equals(zobs_A[(b-1)*6+j,t-1],1)
      colo_A[(b-1)*6+j,(t-1)] <- equals(zobs_A[(b-1)*6+j,t],1) * equals(zobs_A[(b-1)*6+j,t-1],0)
      noncolo_A[(b-1)*6+j,(t-1)] <- equals(zobs_A[(b-1)*6+j,t],0) * equals(zobs_A[(b-1)*6+j,t-1],0)
      
      ext_B[(b-1)*6+j,(t-1)] <- equals(zobs_B[(b-1)*6+j,t],0) * equals(zobs_B[(b-1)*6+j,t-1],1)
      nonext_B[(b-1)*6+j,(t-1)] <- equals(zobs_B[(b-1)*6+j,t],1) * equals(zobs_B[(b-1)*6+j,t-1],1)
      colo_B[(b-1)*6+j,(t-1)] <- equals(zobs_B[(b-1)*6+j,t],1) * equals(zobs_B[(b-1)*6+j,t-1],0)
      noncolo_B[(b-1)*6+j,(t-1)] <- equals(zobs_B[(b-1)*6+j,t],0) * equals(zobs_B[(b-1)*6+j,t-1],0)
    
      # ... for replicated data
      extrep_A[(b-1)*6+j,(t-1)] <- equals(zobsrep_A[(b-1)*6+j,t],0) * equals(zobsrep_A[(b-1)*6+j,t-1],1)
      nonextrep_A[(b-1)*6+j,(t-1)] <- equals(zobsrep_A[(b-1)*6+j,t],1) * equals(zobsrep_A[(b-1)*6+j,t-1],1)
      colorep_A[(b-1)*6+j,(t-1)] <- equals(zobsrep_A[(b-1)*6+j,t],1) * equals(zobsrep_A[(b-1)*6+j,t-1],0)
      noncolorep_A[(b-1)*6+j,(t-1)] <- equals(zobsrep_A[(b-1)*6+j,t],0)*equals(zobsrep_A[(b-1)*6+j,t-1],0)
      
      extrep_B[(b-1)*6+j,(t-1)] <- equals(zobsrep_B[(b-1)*6+j,t],0) * equals(zobsrep_B[(b-1)*6+j,t-1],1)
      nonextrep_B[(b-1)*6+j,(t-1)] <- equals(zobsrep_B[(b-1)*6+j,t],1) * equals(zobsrep_B[(b-1)*6+j,t-1],1)
      colorep_B[(b-1)*6+j,(t-1)] <- equals(zobsrep_B[(b-1)*6+j,t],1) * equals(zobsrep_B[(b-1)*6+j,t-1],0)
      noncolorep_B[(b-1)*6+j,(t-1)] <- equals(zobsrep_B[(b-1)*6+j,t],0)*equals(zobsrep_B[(b-1)*6+j,t-1],0)
        } # end site loop
      } # end block loop
    } # end time loop
    
    # Tally up number of transitions and put into a matrix for each year
    
     for(t in 1:(nseason-1)){
      # ... for observed data
      tm_A[1,1,t] <- sum(noncolo_A[,t]) # transition mat for obs. data
      tm_A[1,2,t] <- sum(colo_A[,t])
      tm_A[2,1,t] <- sum(ext_A[,t])
      tm_A[2,2,t] <- sum(nonext_A[,t])
      
      tm_B[1,1,t] <- sum(noncolo_B[,t]) # transition mat for obs. data
      tm_B[1,2,t] <- sum(colo_B[,t])
      tm_B[2,1,t] <- sum(ext_B[,t])
      tm_B[2,2,t] <- sum(nonext_B[,t])
    
      # ... for replicated data
      tmrep_A[1,1,t] <- sum(noncolorep_A[,t]) # transition mat for rep. data
      tmrep_A[1,2,t] <- sum(colorep_A[,t])
      tmrep_A[2,1,t] <- sum(extrep_A[,t])
      tmrep_A[2,2,t] <- sum(nonextrep_A[,t])
      
      tmrep_B[1,1,t] <- sum(noncolorep_B[,t]) # transition mat for rep. data
      tmrep_B[1,2,t] <- sum(colorep_B[,t])
      tmrep_B[2,1,t] <- sum(extrep_B[,t])
      tmrep_B[2,2,t] <- sum(nonextrep_B[,t])
     }
    
  # Compute expected numbers of transitions under the model
  
  
  # remove block structure
  # first we need to calculate the occupancy of each species at each time step. 
   z_22 <- z2==2
   z_2 <- z_22*1 # make numeric
  
   z_44 <- z2==4
   z_4 <- z_44*1 # make numeric
  
   z_33 <- z2==3
   z_3 <- z_33*1 # make numeric
  
   z_A <- z_2+z_4
   z_B <- z_3+z_4
  
  # calculate occupancy
  for(t in 1:nseason){
  occA[t] <- mean(z_A[,t])
  occB[t] <- mean(z_B[,t])}
  
  # Probability of each individual transition
  for(j in 1:nsite){
    for(b in 1:nblock){
     for(t in 1:(nseason-1)){
    noncolo.exp_A[(b-1)*6+j,t] <- ifelse(z_B[(b-1)*6+j,t]==1, 1-z_A[(b-1)*6+j,t] * (1-gamAB), 1-z_A[(b-1)*6+j,t] * (1-gamA))
    colo.exp_A[(b-1)*6+j,t] <- ifelse(z_B[(b-1)*6+j,t]==1, (1-z_A[(b-1)*6+j,t]) * gamAB, 1-z_A[(b-1)*6+j,t] * gamA) 
    ext.exp_A[(b-1)*6+j,t] <- ifelse(z_B[(b-1)*6+j,t]==1, z_A[(b-1)*6+j,t] * epsAB, z_A[(b-1)*6+j,t] * epsA)
    nonext.exp_A[(b-1)*6+j,t] <- ifelse(z_B[(b-1)*6+j,t]==1, z_A[(b-1)*6+j,t] * (1-epsAB), z_A[(b-1)*6+j,t] * (1-epsA))
    
    noncolo.exp_B[(b-1)*6+j,t] <- ifelse(z_A[(b-1)*6+j,t]==1, 1-z_B[(b-1)*6+j,t] * (1-gamBA), 1-z_B[(b-1)*6+j,t] * (1-gamB))
    colo.exp_B[(b-1)*6+j,t] <- ifelse(z_A[(b-1)*6+j,t]==1, (1-z_B[(b-1)*6+j,t]) * gamBA, 1-z_B[(b-1)*6+j,t] * gamB) 
    ext.exp_B[(b-1)*6+j,t] <- ifelse(z_A[(b-1)*6+j,t]==1, z_B[(b-1)*6+j,t] * epsBA, z_B[(b-1)*6+j,t] * epsB)
    nonext.exp_B[(b-1)*6+j,t] <- ifelse(z_A[(b-1)*6+j,t]==1, z_B[(b-1)*6+j,t] * (1-epsBA), z_B[(b-1)*6+j,t] * (1-epsB))
   } # end season loop
   } # end block loop 
  } # end site loop
  
  for(t in 1:(nseason-1)){
    Etm_A[1,1,t] <- sum(noncolo.exp_A[,t])
    Etm_A[1,2,t] <- sum(colo.exp_A[,t])
    Etm_A[2,1,t] <- sum(ext.exp_A[,t])
    Etm_A[2,2,t] <- sum(nonext.exp_A[,t])
    
    Etm_B[1,1,t] <- sum(noncolo.exp_B[,t])
    Etm_B[1,2,t] <- sum(colo.exp_B[,t])
    Etm_B[2,1,t] <- sum(ext.exp_B[,t])
    Etm_B[2,2,t] <- sum(nonext.exp_B[,t])
  }
  
  # Compute Chi-square discrepancy
  for(t in 1:(nseason-1)){
    # ... for observed data
    x2Open_A[1,1,t] <- pow((tm_A[1,1,t] - Etm_A[1,1,t]), 2) / (tm_A[1,1,t]+e)
    x2Open_A[1,2,t] <- pow((tm_A[1,2,t] - Etm_A[1,2,t]), 2) / (tm_A[1,2,t]+e)
    x2Open_A[2,1,t] <- pow((tm_A[2,1,t] - Etm_A[2,1,t]), 2) / (tm_A[2,1,t]+e)
    x2Open_A[2,2,t] <- pow((tm_A[2,2,t] - Etm_A[2,2,t]), 2) / (tm_A[2,2,t]+e)

    x2Open_B[1,1,t] <- pow((tm_B[1,1,t] - Etm_B[1,1,t]), 2) / (tm_B[1,1,t]+e)
    x2Open_B[1,2,t] <- pow((tm_B[1,2,t] - Etm_B[1,2,t]), 2) / (tm_B[1,2,t]+e)
    x2Open_B[2,1,t] <- pow((tm_B[2,1,t] - Etm_B[2,1,t]), 2) / (tm_B[2,1,t]+e)
    x2Open_B[2,2,t] <- pow((tm_B[2,2,t] - Etm_B[2,2,t]), 2) / (tm_B[2,2,t]+e)
    
    # ... for replicated data
   x2repOpen_A[1,1,t] <- pow((tmrep_A[1,1,t]-Etm_A[1,1,t]),2)/(tmrep_A[1,1,t]+e)
    x2repOpen_A[1,2,t] <- pow((tmrep_A[1,2,t]-Etm_A[1,2,t]),2)/(tmrep_A[1,2,t]+e)
   x2repOpen_A[2,1,t] <- pow((tmrep_A[2,1,t]-Etm_A[2,1,t]),2)/(tmrep_A[2,1,t]+e)
    x2repOpen_A[2,2,t] <- pow((tmrep_A[2,2,t]-Etm_A[2,2,t]),2)/(tmrep_A[2,2,t]+e)
    
    x2repOpen_B[1,1,t] <- pow((tmrep_B[1,1,t]-Etm_B[1,1,t]),2)/(tmrep_B[1,1,t]+e)
    x2repOpen_B[1,2,t] <- pow((tmrep_B[1,2,t]-Etm_B[1,2,t]),2)/(tmrep_B[1,2,t]+e)
    x2repOpen_B[2,1,t] <- pow((tmrep_B[2,1,t]-Etm_B[2,1,t]),2)/(tmrep_B[2,1,t]+e)
    x2repOpen_B[2,2,t] <- pow((tmrep_B[2,2,t]-Etm_B[2,2,t]),2)/(tmrep_B[2,2,t]+e)
  }
  
  
  # Add up overall test statistic and compute fit stat ratio (open part)
  
  Chi2Open_A <- sum(x2Open_A[,,])       # Chisq. statistic for observed data
  Chi2repOpen_A <- sum(x2repOpen_A[,,]) # Chisq. statistic for replicated data
  Chi2ratioOpen_A <- Chi2Open_A / Chi2repOpen_A
  
  Chi2Open_B <- sum(x2Open_B[,,])       # Chisq. statistic for observed data
  Chi2repOpen_B <- sum(x2repOpen_B[,,]) # Chisq. statistic for replicated data
  Chi2ratioOpen_B <- Chi2Open_B / Chi2repOpen_B
  
  # (2b) Computations for the GoF of the closed part of the model
  # (based on the number of times detected, i.e., detection freqiencies)
  # --------------------------------------------------------------------
  # Compute detection frequencies for observed and replicated data

  for (j in 1:nsite){
    for(b in 1:nblock){
    for (t in 1:(nseason-1)){
      # Det. frequencies for observed and replicated data
      detfreq_A[(b-1)*6+j,t] <- sum(y_A[(b-1)*6+j,t,])
      detfreqrep_A[(b-1)*6+j,t] <- sum(yrep_A[(b-1)*6+j,t,])

      detfreq_B[(b-1)*6+j,t] <- sum(y_B[(b-1)*6+j,t,])
      detfreqrep_B[(b-1)*6+j,t] <- sum(yrep_B[(b-1)*6+j,t,])
      
      # Expected detection frequencies under the model
      for (day in 1:nsurvey){
        tmp_A[(b-1)*6+j,day,t] <- z_A[(b-1)*6+j, t] * pA[t]
        tmp_B[(b-1)*6+j,day,t] <- z_B[(b-1)*6+j, t] * pB[t]
      }
      
      E_A[(b-1)*6+j,t] <- sum(tmp_A[(b-1)*6+j,,t])     # Expected number of detections for A
      E_B[(b-1)*6+j,t] <- sum(tmp_B[(b-1)*6+j,,t])     # Expected number of detections for B
      
      # Chi-square and Freeman-Tukey discrepancy measures
      # ..... for actual data set
      x2Closed_A[(b-1)*6+j,t] <- pow((detfreq_A[(b-1)*6+j,t] - E_A[(b-1)*6+j,t]),2) / (E_A[(b-1)*6+j,t] + e)
      ftClosed_A[(b-1)*6+j,t] <- pow((sqrt(detfreq_A[(b-1)*6+j,t]) - sqrt(E_A[(b-1)*6+j,t])),2)
      
      x2Closed_B[(b-1)*6+j,t] <- pow((detfreq_B[(b-1)*6+j,t] - E_B[(b-1)*6+j,t]),2) / (E_B[(b-1)*6+j,t] + e)
      ftClosed_B[(b-1)*6+j,t] <- pow((sqrt(detfreq_B[(b-1)*6+j,t]) - sqrt(E_B[(b-1)*6+j,t])),2)
      
      # ..... for replicated data set
      x2repClosed_A[(b-1)*6+j,t] <- pow((detfreqrep_A[(b-1)*6+j,t] - E_A[(b-1)*6+j,t]),2) / (E_A[(b-1)*6+j,t] + e)
      ftrepClosed_A[(b-1)*6+j,t] <- pow((sqrt(detfreqrep_A[(b-1)*6+j,t]) - sqrt(E_A[(b-1)*6+j,t])),2)
      
      x2repClosed_B[(b-1)*6+j,t] <- pow((detfreqrep_B[(b-1)*6+j,t] - E_B[(b-1)*6+j,t]),2) / (E_B[(b-1)*6+j,t] + e)
      ftrepClosed_B[(b-1)*6+j,t] <- pow((sqrt(detfreqrep_B[(b-1)*6+j,t]) - sqrt(E_B[(b-1)*6+j,t])),2)
    } # end season loop
    } # end block loop
  } # end site loop
  
  # Add up Chi-square and FT discrepancies and compute fit stat ratio (closed part)
  Chi2Closed_A      <- sum(x2Closed_A[,])
  FTClosed_A        <- sum(ftClosed_A[,])
  Chi2repClosed_A   <- sum(x2repClosed_A[,])
  FTrepClosed_A     <- sum(ftrepClosed_A[,])
  Chi2ratioClosed_A <- Chi2Closed_A / Chi2repClosed_A
  FTratioClosed_A   <- FTClosed_A / FTrepClosed_A
  
  Chi2Closed_B      <- sum(x2Closed_B[,])
  FTClosed_B        <- sum(ftClosed_B[,])
  Chi2repClosed_B   <- sum(x2repClosed_B[,])
  FTrepClosed_B     <- sum(ftrepClosed_B[,])
  Chi2ratioClosed_B <- Chi2Closed_B / Chi2repClosed_B
  FTratioClosed_B   <- FTClosed_B / FTrepClosed_B
  
    }# end
    ",fill = TRUE)
sink()

## import data
setwd("./data") # set wd to where the data is stored

#load("occm_mustela_rodent_var_snowbed.rda")    
#load("case_study_data.RData")
load("occm_mustela_rodent_var_snowbed_rmBQ.rda")

yb <-occm_ko3[,,170:173,] # change name of imported object to fit with the rest of the code
#yb <-occm_ko3[,,1:4,] # change name of imported object to fit with the rest of the code

#replace NA's with 0's 
yb[is.na(yb)] <- 1
yb

dim(yb) # check that dimensions are ok

# which data points excist
has_data <- which(
  !is.na(yb),
  arr.ind = TRUE
)

time_id <- has_data[,3]
site_id <- has_data[,1]
block_id <- has_data[,2]

# make table with unique block, site and time combinations
dat <- tibble::as.tibble(cbind(site_id, block_id, time_id))
uniq_site_block_time <- dplyr::distinct(dat)

uniq_site <- dplyr::pull(uniq_site_block_time,1)
uniq_block <- dplyr::pull(uniq_site_block_time,2) 
uniq_time <- dplyr::pull(uniq_site_block_time,3) 

n_unique <- length(uniq_site)

y_long <- yb[!is.na(yb)]

nvisits = length(y_long)

#load cov
load("season.rda") 
season <- season[170:173]
season[season==1] <- 2
season[season==0] <- 1

# give data
data <-list(nseason = dim(yb)[3], nblock = dim(yb)[2], nsite = dim(yb)[1], nsurvey = dim(yb)[4], 
            nout=4, y = y_long, season = season, e = 0.0001, nvisits=nvisits, yb=yb,
            time_id=time_id, block_id=block_id, site_id=site_id)

# naming some parameters for loops further down in this script
nseason = dim(yb)[3]; nblock = dim(yb)[2]; nsite = dim(yb)[1]; nsurvey = dim(yb)[4]

# Initial values for state
sp_inits <- apply(yb,c(1,2,3),max)

# loop to make cases where both state 2 and 3 is observed within the same primary occasion have initial value 4
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

# give initial values
inits=function(){list( 
  z = sp_inits, alphaA0=runif(2), alphaB0=runif(2),
  gamA=runif(1,0.1,0.9), gamB=runif(1,0.1,0.9), gamAB=runif(1,0.1,0.9), gamBA=runif(1,0.1,0.9),
  epsA=runif(1,0.1,0.9), epsB=runif(1,0.1,0.9), epsAB=runif(1,0.1,0.9), epsBA=runif(1,0.1,0.9), 
  GamA=runif(1,0.1,0.9), GamB=runif(1,0.1,0.9), GamAB=runif(1,0.1,0.9), GamBA=runif(1,0.1,0.9),
  EpsA=runif(1,0.1,0.9), EpsB=runif(1,0.1,0.9), EpsAB=runif(1,0.1,0.9), EpsBA=runif(1,0.1,0.9)
)}

# Parameters monitored
params <- c("gamA","gamB","gamAB","gamBA","epsA","epsB","epsAB","epsBA","psi",
            "GamA","GamB","GamAB","GamBA","EpsA","EpsB","EpsAB","EpsBA", "pA","pB","z","x",
            "alphaA0","alphaB0", 
            "diff_gamA", "diff_gamB", "diff_epsA", "diff_epsB", "diff_GamA", "diff_GamB", "diff_EpsA", "diff_EpsB",
            "ratio_gamA", "ratio_gamB", "ratio_epsA", "ratio_epsB", "ratio_GamA", "ratio_GamB", "ratio_EpsA", "ratio_EpsB",
            'Chi2Open_A', 'Chi2repOpen_A','Chi2ratioOpen_A', 'Chi2Closed_A', 'Chi2repClosed_A', 'Chi2ratioClosed_A',
            'FTClosed_A', 'FTrepClosed_A', 'FTratioClosed_A', 'tm_A', 'tmrep_A', 'Etm_A', 
            'Chi2Open_B', 'Chi2repOpen_B','Chi2ratioOpen_B', 'Chi2Closed_B', 'Chi2repClosed_B', 'Chi2ratioClosed_B',
            'FTClosed_B', 'FTrepClosed_B', 'FTratioClosed_B', 'tm_B', 'tmrep_B', 'Etm_B')

# MCMC settings
ni <- 5   ;   nt <- 1   ;   nb <- 0   ;   nc <- 1    ;   na <- 0
#ni <- 100000   ;   nt <- 20   ;   nb <- 25000   ;   nc <- 4    ;   na <- 25000

# run model in jags
setwd("../")

va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ <- jags(data, inits=inits, params, "mod_spatial_dynocc_bayp.txt", n.chains = nc,
                                                     n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)

# Save model
#setwd("./model_output")
#save(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ, file="va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rm_BQ.rda")

#~ End of script