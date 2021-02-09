
# set working directory
setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel")
setwd("./models/hidden_block_sim")

# Specify model in JAGS language
sink("mod.txt")
cat("
    model{
    
    ############################################################################################
    # State 1= Unoccupied(U), State 2= vole(A), State 3 = stoat(B), State 4 = vole & stoat(AB) #  
    ############################################################################################
   
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
