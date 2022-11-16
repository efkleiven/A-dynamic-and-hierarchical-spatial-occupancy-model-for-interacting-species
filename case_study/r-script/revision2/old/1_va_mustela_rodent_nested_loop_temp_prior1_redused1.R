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
#setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel")
#setwd("./models/hidden_block_ko_fromautoclass/ko_vj")

#setwd("~/UiT/GitProjects/A-dynamic-and-hierarchical-spatial-occupancy-model-for-interacting-species/case_study")
setwd("C:/Users/ekl013/OneDrive - UiT Office 365/GitProjects/A-dynamic-and-hierarchical-spatial-occupancy-model-for-interacting-species/case_study")

#setwd("H:/UiT/GitProjects/A-dynamic-and-hierarchical-spatial-occupancy-model-for-interacting-species/case_study")

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

          for(i in 1:nvisits) {
          
            y[i] ~ dcat( dpm[site_id[i], block_id[i], time_id[i], ( 1:nout ) , z[site_id[i], block_id[i], time_id[i]]] + 0.005)  # +0.005 to avoide giving the dcat a prob of 0 
          
           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           # Draw a replicate data set under the fitted model
   
            yrep[i] ~ dcat(dpm[site_id[i], block_id[i], time_id[i], , z[site_id[i], block_id[i], time_id[i]]] + 0.005)
            
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
        z[j, b, t+1] ~ dcat( stpm[( 1:nout ) , z[ j, b, t], x[b,t+1]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
         
      # remove block structure from z for use in Goodness of Fit check
        z2[(b-1)*6+j, t+1] <- z[j, b, t+1] 

         } #end time loop
      } # end site loop

    for(t in 1:(nseason)){      
      for(j in 1:nsite){
    
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
#  yrep_22 <- yrep2==2       # select data point where status 2 (only A observed)
#  yrep_2 <- yrep_22 * 1      # make numeric

#  yrep_44 <- yrep2==4       # select data point where status 4 ( A and B observed)
#  yrep_4 <- yrep_44*1      # make numeric

#  yrep_A <- yrep_2+yrep_4 

#  yrep_33 <- yrep2==3       # select data point where status 2 (only A observed)
#  yrep_3 <- yrep_33*1      # make numeric

#  yrep_B <- yrep_3+yrep_4 

#  y_22 <- y2==2       # select data point where status 2 (only A observed)
#  y_2 <- y_22*1      # make numeric

#  y_44 <- y2==4       # select data point where status 4 ( A and B observed)
#  y_4 <- y_44*1      # make numeric

#  y_A <- y_2+y_4 

#  y_33 <- y2==3       # select data point where status 2 (only A observed)
#  y_3 <- y_33*1      # make numeric

# y_B <- y_3+y_4 
 
 
    }# end
    ",fill = TRUE)
sink()

## import data


#setwd("./data") # set wd to where the data is stored
setwd("./data/revision2")

#load("occm_mustela_rodent_var_snowbed_rmBQ.rda")
#load("occm_rodent_coatdata_2021.rda")

load("occm_vole_mustelid_snowbed_2016_2021.rda")

#yb <-occm_ko3 # change name of imported object to fit with the rest of the code
yb <-occm_va[,,170:200,] # change name of imported object to fit with the rest of the code

summary(yb)
dim(yb)
yb

#nas <- which(is.na(yb)) # funker ikke!

#replace NA's with 0's 
#yb2 <- yb

#yb[is.na(yb)] <- 1

dim(yb) # check that dimensions are ok

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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Load covariates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


setwd("../")
#load("season_cov_from_temp.rda")
load("season_cov_fromTemp_2021.rda")

# restructuring sites to match the occupancy data frame
season1 <- season_cov[6,3:4,] 
season_cov[12,3:4, ] <- season1

season <- season_cov[7:12, , ]
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
            time_id=time_id, block_id=block_id, site_id=site_id, day_id=day_id)

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
            "ratio_gamA", "ratio_gamB", "ratio_epsA", "ratio_epsB", "ratio_GamA", "ratio_GamB", "ratio_EpsA", "ratio_EpsB")


# MCMC settings
ni <- 20   ;   nt <- 1   ;   nb <- 1  ;   nc <- 1    ;   na <- 0

# run model in jags
setwd("../")

#ptm2 <- proc.time()
out.gof1 <- jags(data, inits=inits, parameters.to.save=params, "mod_spatial_dynocc_bayp.txt", n.chains = nc,
                    n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)
#proc.time() - ptm2                

# Save model
#setwd("./model_output/revision2")
#save(out.gof1, file="va_snowbed_mustela_rodent_gof_nestedloop_temp_pri1_005_2016_2021.rda")
#~ End of script