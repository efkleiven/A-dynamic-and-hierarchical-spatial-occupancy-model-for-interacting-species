###################################################################################################################
##  A dynamic occupancy model for interacting species with two spatial scales model analyzing data from          ##  
##  a long-term monitoring program of small mammals on the arctic tundra                                         ##
##                 by Eivind Flittie Kleiven                                                                     ##
##                                                                                                               ##
###################################################################################################################

# a part of the code is modified from https://github.com/mikemeredith/AHM_code/blob/master/AHM2_ch04/AHM2_04.08.R

# Call jags(and other packages)
rm(list=ls())

library(jagsUI)

# set working directory
setwd("C:/Users/ekl013/OneDrive - UiT Office 365/GitProjects/A-dynamic-and-hierarchical-spatial-occupancy-model-for-interacting-species/case_study")

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

          for(i in 1:nvisits) {
          
            y[i] ~ dcat( dpm[site_id[i], block_id[i], time_id[i], ( 1:nout ) , z[site_id[i], block_id[i], time_id[i]]] + 0.0005)  # +0.005 to avoide giving the dcat a prob of 0 
          
           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           # Draw a replicate data set under the fitted model
   
            yrep[i] ~ dcat(dpm[site_id[i], block_id[i], time_id[i], , z[site_id[i], block_id[i], time_id[i]]] + 0.0005)
            
          #derived parameters for Goodness-of-Fit check
            y2[(block_id[i]-1)*6+site_id[i],time_id[i], day_id[i]] <- y[i]
            
            yrep2[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]] <- yrep[i]
           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          } #end nvisit loop
          
          
    # latent block state for the rest of the seasons
    for(b in 1:nblock){
      for(t in 1:(nseason-1)){
        x[b, t+1] ~ dcat(btpm[(1:nout), x[b, t]])

      # latent site state for the rest of the seasons
      for(j in 1:nsite){
        z[j, b, t+1] ~ dcat( stpm[( 1:nout ) , z[ j, b, t], x[b,t+1]] + 0.0005)  # +0.01 to avoide giving the dcat a prob of 0 
         


         } #end site loop
      } # end time loop

    for(t in 1:(nseason)){      
      for(j in 1:nsite){
    
      # remove block structure from z for use in Goodness of Fit check
        z2[(b-1)*6+j, t] <- z[j, b, t] 
    #############################################################
    ## detection matrix (OS = observed state, TS = true state) ##
    #############################################################
          
        # TS = U
        dpm[j, b, t, 1, 1] <- 1                                 #--|OS = U
        dpm[j, b, t, 2, 1] <- 0                                 #--|OS = A
        dpm[j, b, t, 3, 1] <- 0                                 #--|OS = B
        dpm[j, b, t, 4, 1] <- 0                                 #--|OS = AB
    
        # TS = A
        dpm[j, b, t, 1, 2] <- 1-pA[j, b, t]                     #--|OS = U
        dpm[j, b, t, 2, 2] <- pA[j, b, t]                       #--|OS = A
        dpm[j, b, t, 3, 2] <- 0                                 #--|OS = B
        dpm[j, b, t, 4, 2] <- 0                                 #--|OS = AB
    
        # TS = B
        dpm[j, b, t, 1, 3] <- 1-pB[j, b, t]                     #--|OS = U
        dpm[j, b, t, 2, 3] <- 0                                 #--|OS = A
        dpm[j, b, t, 3, 3] <- pB[j, b, t]                       #--|OS = B
        dpm[j, b, t, 4, 3] <- 0                                 #--|OS = AB
    
        # TS = AB
        dpm[j, b, t, 1, 4] <- (1-pA[j, b, t]) * (1-pB[j, b, t])  #--|OS = U
        dpm[j, b, t, 2, 4] <- pA[j, b, t] * (1-pB[j, b, t])      #--|OS = A
        dpm[j, b, t, 3, 4] <- (1-pA[j, b, t]) * pB[j, b, t]      #--|OS = B
        dpm[j, b, t, 4, 4] <- pA[j, b, t] * pB[j, b, t]          #--|OS = AB
    
        ## logit links for detection probs
        logit(pA[j, b, t]) <- alphaA0[season[j, b, t]]
        logit(pB[j, b, t]) <- alphaB0[season[j, b, t]]
     
         } #end site loop
      } # end time loop
    } # end block loop  

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
    
    
  # (2a) Computations for the GoF of the open part of the model
  # (based on number of state transitions)
  # ----------------------------------------------------------

# seperate species
for(i in 1:nvisits) {
y_A[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]] <- ifelse(y2[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]]==2 || y2[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]]==4, 1, 0)
y_B[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]] <- ifelse(y2[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]]==3 || y2[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]]==4, 1, 0)

yrep_A[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]] <- ifelse(yrep2[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]]==2 || yrep2[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]]==4, 1, 0)
yrep_B[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]] <- ifelse(yrep2[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]]==3 || yrep2[(block_id[i]-1)*6+site_id[i], time_id[i], day_id[i]]==4, 1, 0)
} # end loop

  # Compute observed z matrix for observed and replicated data
for(i in 1:nz) {
    zobs_A[(block_id2[i]-1) * 6 + site_id2[i], time_id2[i]] <- max(y_A[(block_id2[i]-1) * 6 + site_id2[i], time_id2[i], ]) 
    zobsrep_A[(block_id2[i]-1)*6+site_id2[i], time_id2[i]] <- max(yrep_A[(block_id2[i]-1)*6+site_id2[i], time_id2[i], ]) # For replicated data
       
    zobs_B[(block_id2[i]-1)*6+site_id2[i], time_id2[i]] <- max(y_B[(block_id2[i]-1)*6+site_id2[i], time_id2[i],])       # For observed data
    zobsrep_B[(block_id2[i]-1)*6+site_id2[i], time_id2[i]] <- max(yrep_B[(block_id2[i]-1)*6+site_id2[i], time_id2[i],]) # For replicated data
    
    z_A[(block_id2[i]-1)*6+site_id2[i], time_id2[i]] <- ifelse(z2[(block_id2[i]-1)*6+site_id2[i], time_id2[i]]==2 || z2[(block_id2[i]-1)*6+site_id2[i], time_id2[i]]==4, 1, 0)
    z_B[(block_id2[i]-1)*6+site_id2[i], time_id2[i]] <- ifelse(z2[(block_id2[i]-1)*6+site_id2[i], time_id2[i]]==3 || z2[(block_id2[i]-1)*6+site_id2[i], time_id2[i]]==4, 1, 0)

} # end loop


    # Identify extinctions, persistence, colonization and non-colonizations
for(i in 1:nz2){
      # ... for observed data
      ext_A[(block_id3[i]-1)*6+site_id3[i], time_id3[i]] <- equals(zobs_A[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])], 0) * equals(zobs_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1], 1)
     
      nonext_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobs_A[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],1) * equals(zobs_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],1)
        
      colo_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <-  equals(zobs_A[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],1) * equals(zobs_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],0)
      
      noncolo_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobs_A[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],0) * equals(zobs_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],0)
      
      ext_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobs_B[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],0) * equals(zobs_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],1)
        
      nonext_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobs_B[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],1) * equals(zobs_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],1)
      
      colo_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobs_B[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],1) * equals(zobs_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],0)
        
      noncolo_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobs_B[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],0) * equals(zobs_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],0)
    
      # ... for replicated data
      extrep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobsrep_A[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],0) * equals(zobsrep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],1)
        
      nonextrep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobsrep_A[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],1) * equals(zobsrep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],1)
        
      colorep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobsrep_A[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],1) * equals(zobsrep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],0)
        
      noncolorep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobsrep_A[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],0) * equals(zobsrep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],0)
      
      extrep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobsrep_B[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],0) * equals(zobsrep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],1)
        
      nonextrep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobsrep_B[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],1) * equals(zobsrep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],1)
        
      colorep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobsrep_B[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],1) * equals(zobsrep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],0)
        
      noncolorep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- equals(zobsrep_B[(block_id3[i]-1)*6+site_id3[i],(time_id3[i])],0) * equals(zobsrep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]-1],0)
 
      
      # Probability of each individual transition
     
      noncolo.exp_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- ifelse(z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]==1, 1-z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * (1-gamAB), 1-z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * (1-gamA))
        
      colo.exp_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <-   ifelse(z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]==1, (1-z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]) * gamAB, 1-z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * gamA)
        
      ext.exp_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- ifelse(z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]==1, z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * epsAB, z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * epsA)
        
      nonext.exp_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- ifelse(z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]==1, z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * (1-epsAB), z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * (1-epsA))
    
      noncolo.exp_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- ifelse(z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]==1, 1-z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * (1-gamBA), 1-z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * (1-gamB))
        
      colo.exp_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- ifelse(z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]==1, (1-z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]) * gamBA, 1-z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * gamB)
        
      ext.exp_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- ifelse(z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]==1, z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * epsBA, z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * epsB)
        
      nonext.exp_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- ifelse(z_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]==1, z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * (1-epsBA), z_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] * (1-epsB))


     # Det. frequencies for observed and replicated data
      detfreq_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- sum(y_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i],])
      detfreqrep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- sum(yrep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i],])

      detfreq_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- sum(y_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i],])
      detfreqrep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- sum(yrep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i],])
      
      
      # Expected detection frequencies under the model
      for (day in 1:nsurvey){
        tmp_A[(block_id3[i]-1)*6 + site_id3[i], day, time_id3[i]] <- z_A[(block_id3[i]-1)*6 + site_id3[i], time_id3[i]] * pA[site_id3[i],block_id3[i],time_id3[i]]
        tmp_B[(block_id3[i]-1)*6 + site_id3[i], day, time_id3[i]] <- z_B[(block_id3[i]-1)*6 + site_id3[i], time_id3[i]] * pB[site_id3[i],block_id3[i],time_id3[i]]
      } # end day loop
      
      E_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- sum(tmp_A[(block_id3[i]-1)*6+site_id3[i],,time_id3[i]])     # Expected number of detections for A
      E_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- sum(tmp_B[(block_id3[i]-1)*6+site_id3[i],,time_id3[i]])     # Expected number of detections for B
      
      # Chi-square and Freeman-Tukey discrepancy measures
      # ..... for actual data setblock_id3[i]
      x2Closed_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- pow((detfreq_A[(block_id3[i]-1)*6 + site_id3[i], time_id3[i]] - E_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]),2) / (E_A[(block_id3[i]-1)*6 + site_id3[i],time_id3[i]] + e)
        
      ftClosed_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- pow((sqrt(detfreq_A[(block_id3[i]-1)*6 + site_id3[i], time_id3[i]]) - sqrt(E_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]])),2)
      
      x2Closed_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- pow((detfreq_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] - E_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]),2) / (E_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] + e)
        
      ftClosed_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- pow((sqrt(detfreq_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]) - sqrt(E_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]])),2)
      
      # ..... for replicated data set
      x2repClosed_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- pow((detfreqrep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] - E_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]),2) / (E_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] + e)
        
      ftrepClosed_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- pow((sqrt(detfreqrep_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]) - sqrt(E_A[(block_id3[i]-1)*6+site_id3[i],time_id3[i]])),2)
      
      x2repClosed_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- pow((detfreqrep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] - E_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]),2) / (E_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] + e)
        
      ftrepClosed_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]] <- pow((sqrt(detfreqrep_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]]) - sqrt(E_B[(block_id3[i]-1)*6+site_id3[i],time_id3[i]])),2)
} # end loop  


  for(t in 2:nseason){
  
    # Tally up number of transitions and put into a matrix for each year
    
    # ... for observed data
    tm_A[1,1,t] <- sum(noncolo_A[ind4[1:indl[t], t], t]) # transition mat for obs. data
    tm_A[1,2,t] <- sum(colo_A[ind4[1:indl[t], t], t])
    tm_A[2,1,t] <- sum(ext_A[ind4[1:indl[t], t], t])
    tm_A[2,2,t] <- sum(nonext_A[ind4[1:indl[t], t], t])
      
     tm_B[1,1,t] <- sum(noncolo_B[ind4[1:indl[t], t], t]) # transition mat for obs. data
    tm_B[1,2,t] <- sum(colo_B[ind4[1:indl[t], t], t])
    tm_B[2,1,t] <- sum(ext_B[ind4[1:indl[t], t], t])
    tm_B[2,2,t] <- sum(nonext_B[ind4[1:indl[t], t], t])
    
    # ... for replicated data
    tmrep_A[1,1,t] <- sum(noncolorep_A[ind4[1:indl[t], t], t]) # transition mat for rep. data
    tmrep_A[1,2,t] <- sum(colorep_A[ind4[1:indl[t], t], t])
    tmrep_A[2,1,t] <- sum(extrep_A[ind4[1:indl[t], t], t])
    tmrep_A[2,2,t] <- sum(nonextrep_A[ind4[1:indl[t], t], t])
    
    tmrep_B[1,1,t] <- sum(noncolorep_B[ind4[1:indl[t], t], t]) # transition mat for rep. data
    tmrep_B[1,2,t] <- sum(colorep_B[ind4[1:indl[t], t], t])
    tmrep_B[2,1,t] <- sum(extrep_B[ind4[1:indl[t], t], t])
    tmrep_B[2,2,t] <- sum(nonextrep_B[ind4[1:indl[t], t], t])
  
    Etm_A[1,1,t] <- sum(noncolo.exp_A[ind4[1:indl[t], t], t])
    Etm_A[1,2,t] <- sum(colo.exp_A[ind4[1:indl[t], t], t])
    Etm_A[2,1,t] <- sum(ext.exp_A[ind4[1:indl[t], t], t])
    Etm_A[2,2,t] <- sum(nonext.exp_A[ind4[1:indl[t], t], t])
    
    Etm_B[1,1,t] <- sum(noncolo.exp_B[ind4[1:indl[t], t], t])
    Etm_B[1,2,t] <- sum(colo.exp_B[ind4[1:indl[t], t], t])
    Etm_B[2,1,t] <- sum(ext.exp_B[ind4[1:indl[t], t], t])
    Etm_B[2,2,t] <- sum(nonext.exp_B[ind4[1:indl[t], t], t])
   
  # Compute Chi-square discrepancy

    # ... for observed data
    x2Open_A[1,1,t] <- pow((tm_A[1,1,t] - Etm_A[1,1,t]), 2) / (Etm_A[1,1,t]+e)
    x2Open_A[1,2,t] <- pow((tm_A[1,2,t] - Etm_A[1,2,t]), 2) / (Etm_A[1,2,t]+e)
    x2Open_A[2,1,t] <- pow((tm_A[2,1,t] - Etm_A[2,1,t]), 2) / (Etm_A[2,1,t]+e)
    x2Open_A[2,2,t] <- pow((tm_A[2,2,t] - Etm_A[2,2,t]), 2) / (Etm_A[2,2,t]+e)

    x2Open_B[1,1,t] <- pow((tm_B[1,1,t] - Etm_B[1,1,t]), 2) / (Etm_B[1,1,t]+e)
    x2Open_B[1,2,t] <- pow((tm_B[1,2,t] - Etm_B[1,2,t]), 2) / (Etm_B[1,2,t]+e)
    x2Open_B[2,1,t] <- pow((tm_B[2,1,t] - Etm_B[2,1,t]), 2) / (Etm_B[2,1,t]+e)
    x2Open_B[2,2,t] <- pow((tm_B[2,2,t] - Etm_B[2,2,t]), 2) / (Etm_B[2,2,t]+e)
    
    # ... for replicated data
    x2repOpen_A[1,1,t] <- pow((tmrep_A[1,1,t]-Etm_A[1,1,t]),2)/(Etm_A[1,1,t]+e)
    x2repOpen_A[1,2,t] <- pow((tmrep_A[1,2,t]-Etm_A[1,2,t]),2)/(Etm_A[1,2,t]+e)
    x2repOpen_A[2,1,t] <- pow((tmrep_A[2,1,t]-Etm_A[2,1,t]),2)/(Etm_A[2,1,t]+e)
    x2repOpen_A[2,2,t] <- pow((tmrep_A[2,2,t]-Etm_A[2,2,t]),2)/(Etm_A[2,2,t]+e)
    
    x2repOpen_B[1,1,t] <- pow((tmrep_B[1,1,t]-Etm_B[1,1,t]),2)/(Etm_B[1,1,t]+e)
    x2repOpen_B[1,2,t] <- pow((tmrep_B[1,2,t]-Etm_B[1,2,t]),2)/(Etm_B[1,2,t]+e)
    x2repOpen_B[2,1,t] <- pow((tmrep_B[2,1,t]-Etm_B[2,1,t]),2)/(Etm_B[2,1,t]+e)
    x2repOpen_B[2,2,t] <- pow((tmrep_B[2,2,t]-Etm_B[2,2,t]),2)/(Etm_B[2,2,t]+e)
    
    Chi2Closed_A2[t]    <- sum(x2Closed_A[ind4[1:indl[t], t], t])
    FTClosed_A2[t]        <- sum(ftClosed_A[ind4[1:indl[t], t], t])
    Chi2repClosed_A2[t]   <- sum(x2repClosed_A[ind4[1:indl[t], t], t])
    FTrepClosed_A2[t]     <- sum(ftrepClosed_A[ind4[1:indl[t], t], t])
    
    Chi2Closed_B2[t]      <- sum(x2Closed_B[ind4[1:indl[t], t], t])
    FTClosed_B2[t]        <- sum(ftClosed_B[ind4[1:indl[t], t], t])
    Chi2repClosed_B2[t]   <- sum(x2repClosed_B[ind4[1:indl[t], t], t])
    FTrepClosed_B2[t]     <- sum(ftrepClosed_B[ind4[1:indl[t], t], t])
  } # end time loop

  # Add up overall test statistic and compute fit stat ratio (open part)
  Chi2Open_A <- sum(x2Open_A[,,2:nseason])       # Chisq. statistic for observed data
  Chi2repOpen_A <- sum(x2repOpen_A[,,2:nseason]) # Chisq. statistic for replicated data
  Chi2ratioOpen_A <- Chi2Open_A / Chi2repOpen_A
  
  Chi2Open_B <- sum(x2Open_B[,,2:nseason])       # Chisq. statistic for observed data
  Chi2repOpen_B <- sum(x2repOpen_B[,,2:nseason]) # Chisq. statistic for replicated data
  Chi2ratioOpen_B <- Chi2Open_B / Chi2repOpen_B
  
  # Add up Chi-square and FT discrepancies and compute fit stat ratio (closed part)
  Chi2Closed_A    <- sum(Chi2Closed_A2[2:nseason])
  FTClosed_A        <- sum(FTClosed_A2[2:nseason])
  Chi2repClosed_A   <- sum(Chi2repClosed_A2[2:nseason])
  FTrepClosed_A     <- sum(FTrepClosed_A2[2:nseason])
  Chi2ratioClosed_A <- Chi2Closed_A / Chi2repClosed_A
  FTratioClosed_A   <- FTClosed_A / FTrepClosed_A
  
  Chi2Closed_B      <- sum(Chi2Closed_B2[2:nseason])
  FTClosed_B        <- sum(FTClosed_B2[2:nseason])
  Chi2repClosed_B   <- sum(Chi2repClosed_B2[2:nseason])
  FTrepClosed_B     <- sum(FTrepClosed_B2[2:nseason])
  Chi2ratioClosed_B <- Chi2Closed_B / Chi2repClosed_B
  FTratioClosed_B   <- FTClosed_B / FTrepClosed_B

    }# end
    ",fill = TRUE)
sink()

## import data


#setwd("./data") # set wd to where the data is stored
setwd("./data/revision2")

load("occm_vole_mustelid_snowbed_2016_2021.rda")

yb <-occm_va # change name of imported object to fit with the rest of the code

summary(yb)
dim(yb) # check that dimensions are ok

# manually remove seasons with only missing data
yb <- yb[,,-c(45,105,149,150),]

#define some parameters for later use in dataformating
nseason = dim(yb)[3]
nblock = dim(yb)[2]
nsite = dim(yb)[1]
nsurvey = dim(yb)[4]

# which data points exist
has_data <- which(
  !is.na(yb),
  arr.ind = TRUE
)

time_id <- has_data[,3]
site_id <- has_data[,1]
block_id <- has_data[,2]
day_id <- has_data[,4]

# Data in long formate without missing observations
y_long <- yb[!is.na(yb)]

# length of longdata
nvisits = length(y_long)

# make indeks for complete cases (data from all site in all days)
dim(yb)

# need to subset the dataset to only include primary occassion that have data for all secondary occassions

dat <- array(NA, dim=c(nsite,nblock,nseason))
dat2 <- array(NA, dim=c(nsite,nblock,nseason))

# find complete cases
for(i in 1:nsite){
  for(b in 1:nblock){
dat[i,b,1:nseason]<- complete.cases(yb[i,b,,])
  }}

# find complete cases
for(i in 1:nsite){
  for(b in 1:nblock){
    for(t in 1:nseason){
    dat2[i,b,t]<- !any(is.na(yb[i,b,t,]))
  }}}

# make indecied of primary occ to remove 
ind <- which(dat==T, arr.ind=T)
length(ind)

time_id2 <- ind[,3]
site_id2 <- ind[,1]
block_id2 <- ind[,2]

nz <- length(time_id2) 

# find complete cases that also have complete cases in previous time step
dat3 <- array(NA, dim=c(nsite,nblock,nseason))

for(i in 1:nsite){
  for(b in 1:nblock){
    for(t in 2:nseason){
      dat3[i,b,t]<- !any(is.na(yb[i,b,t,])) & !any(is.na(yb[i,b,t-1,]))  
    }}}

ind3 <- which(dat3==T, arr.ind=T)

time_id3 <- ind3[,3]
site_id3 <- ind3[,1]
block_id3 <- ind3[,2]

nz2 <- length(time_id3) 

# make index for sum-function in JAGS

library(tidyverse)

ind3 <- as_tibble(ind3)

ind4 <- array(NA, dim=c(nsite*nblock,nseason))
indl <- array(NA, dim=c(nseason))

for(t in 2:nseason){
  dat1 <- filter(ind3, dim3==t)
  indl[t] <- length(as.numeric(unlist((dat1[,2]-1)*6 + dat1[,1])))
  ind4[1:indl[t],t] <- as.numeric(unlist((dat1[,2]-1)*6 + dat1[,1]))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Load covariates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

setwd("../")
#load("season_cov_from_temp.rda")
load("season_cov_fromTemp_2021.rda")

# restructuring sites to match the occupancy data frame
season1 <- season_cov[6,3:4,] 
season_cov[12,3:4, ] <- season1

# manually remove seasons with only missing data
season <- season_cov[,,-c(45,105,149,150)]

# make sure dimentions are the same as for observations
season <- season_cov[7:12, ,1:nseason]

# contains NA's but we do not have any data points for these 
summary(season)

# replace NA's with mean of the same site from other blocks

for(i in 1:dim(season)[1]){
  for(j in 1:dim(season)[3]){
    for(k in 1:dim(season)[2])
    if(is.na(season[i,k,j])){
    season[i,k,j] <- ceiling(mean(season[i,,j], na.rm=T))  
    }
  }}

summary(season)

# give data
data <-list(nseason = dim(yb)[3], nblock = dim(yb)[2], nsite = dim(yb)[1], nsurvey = dim(yb)[4], 
            nout=4, y = y_long, season = season, e = 0.0001, nvisits=nvisits, 
            time_id=time_id, block_id=block_id, site_id=site_id, day_id=day_id,
            nz=nz, time_id2=time_id2, block_id2=block_id2, site_id2=site_id2,
            nz2=nz2, time_id3=time_id3, block_id3=block_id3, site_id3=site_id3,
            ind4=ind4, indl=indl)

#yb=yb2,

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
            'Chi2Open_A', 'Chi2repOpen_A','Chi2ratioOpen_A', 'Chi2Closed_A', 'Chi2repClosed_A', 'Chi2ratioClosed_A', 'tm_A', 'tmrep_A', 'Etm_A', 
            'Chi2Open_B', 'Chi2repOpen_B','Chi2ratioOpen_B', 'Chi2Closed_B', 'Chi2repClosed_B', 'Chi2ratioClosed_B',  'tm_B', 'tmrep_B', 'Etm_B',
            'FTClosed_B', 'FTrepClosed_B', 'FTratioClosed_B','FTClosed_A', 'FTrepClosed_A', 'FTratioClosed_A')


# MCMC settings
ni <- 20000   ;   nt <- 10   ;   nb <- 1000  ;   nc <- 6    ;   na <- 5000

# run model in jags
setwd("../")

out.gof1 <- jags(data, inits=inits, parameters.to.save=params, "mod_spatial_dynocc_bayp.txt", n.chains = nc,
                    n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)

# Save model
setwd("./model_output/revision2")
save(out.gof1, file="va_snowbed_mustela_rodent_gof_pri1_0005_2016_2021_2.rda")
#~ End of script