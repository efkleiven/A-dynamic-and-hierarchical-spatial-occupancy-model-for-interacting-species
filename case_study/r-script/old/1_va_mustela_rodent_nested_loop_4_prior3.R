###################################################################################################################
##  A dynamic occupancy model for interacting species with two spatial scales model analyzing data from         ##  
##  a long-term monitoring program of small mammals on the arctic tundra                                         ##
##                 by Eivind Flittie Kleiven and Frederic Barraquand                                             ##
##    Last updated 27.1.21                                                                                       ##
###################################################################################################################

# a part of the code is modified from https://github.com/mikemeredith/AHM_code/blob/master/AHM2_ch04/AHM2_04.08.R

rm(list=ls())

# Call jags(and other packages)
library(jagsUI)

# set working directory
setwd("~/UiT/GitProjects/A-dynamic-and-hierarchical-spatial-occupancy-model-for-interacting-species/case_study")

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
    gamA  ~ dbeta(4,2)
    gamB  ~ dbeta(2,4)
    gamAB ~ dbeta(2,4)
    gamBA ~ dbeta(4,2)
    epsA  ~ dbeta(2,4)
    epsB  ~ dbeta(4,2)
    epsAB ~ dbeta(4,2)
    epsBA ~ dbeta(2,4)
    
    # block parameters
    GamA  ~ dbeta(2,4)
    GamB  ~ dbeta(2,4)
    GamAB ~ dbeta(2,4)
    GamBA ~ dbeta(2,4)
    EpsA  ~ dbeta(2,4)
    EpsB  ~ dbeta(4,2)
    EpsAB ~ dbeta(4,2)
    EpsBA ~ dbeta(2,4)
    
    # interscept det prob
    for(i in 1:2){    
      alphaA0[i] ~ dnorm(0.5,1)
      alphaB0[i] ~ dnorm(-0.5,1)
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
          
            y[i] ~ dcat( dpm[time_id[i], ( 1:nout ) , z[site_id[i], block_id[i], time_id[i]]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
          
           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           # Draw a replicate data set under the fitted model
   
              #yrep[i] ~ dcat( dpm[time_id[i], , z[site_id[i], block_id[i], time_id[i]]] + 0.001 )

             #y2[(block_id[i]-1)*12 +   site_id[i], time_id[i], day_id[i]] <- yb[site_id[i], block_id[i], time_id[i], day_id[i]]
            #yrep2[(block_id[i]-1)*12 + site_id[i], time_id[i], day_id[i]] <- yrep[i]
           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          } #end nvisit loop
          
  # latent site state for the rest of the seasons
  
  for(q in 1:length(uniq_site1_t0)){
    z[uniq_site1_t0[q], uniq_block1_t0[q], uniq_time1_t0[q]+1] ~ dcat( stpm[( 1:nout ) , z[ uniq_site1_t0[q], uniq_block1_t0[q], uniq_time1_t0[q]], x[uniq_block1_t0[q], uniq_time1_t0[q]+1]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
    
    # remove block structure from z
    #z2[(uniq_block1_t0[q]-1)*12+uniq_site1_t0[q], uniq_time1_t0[q]] <- z[uniq_site1_t0[q], uniq_block1_t0[q], uniq_time1_t0[q]] 
    } # end uniq1 loop
  
  for(w in 1:length(uniq_block2_t0)){
    # latent block state for the rest of the seasons 
    x[uniq_block2_t0[w], uniq_time2_t0[w]+1] ~ dcat(btpm[(1:nout), x[uniq_block2_t0[w], uniq_time2_t0[w]]])
  }
  
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
    
  
    }# end
    ",fill = TRUE)
sink()

## import data
setwd("./data") # set wd to where the data is stored

load("occm_mustela_rodent_var_snowbed_rmBQ.rda")   

yb <-occm_ko3 # change name of imported object to fit with the rest of the code

dim(yb) # check that dimensions are ok

#replace NA's with 0's 
yb[is.na(yb)] <- 5

# which data points excist
has_data <- as.data.frame(which(
  !is.na(yb),
  arr.ind = TRUE
))


# order has_data, first on season, then day, then block then site.
has_data <- has_data[order(has_data$dim3),]
has_data <- has_data[order(has_data$dim4),]
has_data <- has_data[order(has_data$dim2),]
has_data <- has_data[order(has_data$dim1),]

head(has_data)

# make find all t-1 
#has_data_t0 <- cbind(as.integer(has_data[,1]), as.integer(has_data[,2]), as.integer(has_data[,3]-1), as.integer(has_data[,4]))

#head(has_data_t0)

# remove rowes where t=0
#has_data_t0 <- has_data_t0[!has_data_t0[,3]==0,]

# match has_data with has_data_t0(t-1)
#library(prodlim)
#ro <- row.match(as.data.frame(has_data_t0), as.data.frame(has_data))
#length(ro)

# ro is a vector giving the row numbers of has_data that matched the given row of has_data_t0, if no match found it will be an NA
# now extract all rows of has_data_t0 where there was no match (meaning that they are not already in has_data)

#table(is.na(ro))

#row_add <- has_data_t0[is.na(ro),]

# add rows to has_data
row_add <- as.data.frame(row_add)
names(row_add) <- names(has_data)

has_data <- rbind(has_data, row_add)

dim(has_data)
has_data[1,]

# make long data
y_long <- yb[!is.na(yb)]

# add t-1 seaons where they are missing
y_long <- c(y_long,rep(1, times=nrow(row_add)))

# add long data to has_data before ordering
has_data2 <- cbind(has_data, y_long)
tail(has_data2)

# order data, first on season, then day, then block then site.
has_data3 <- has_data2[order(has_data2$dim3),]
has_data3 <- has_data2[order(has_data2$dim4),]
has_data3 <- has_data2[order(has_data2$dim2),]
has_data3 <- has_data2[order(has_data2$dim1),]

# remove duplicated rows
library(dplyr)
has_data4 <- has_data3 %>% distinct(dim1, dim2, dim3,dim4, .keep_all = TRUE)

# extract time, site, block and day id
time_id  <- has_data4[,3]
site_id  <- has_data4[,1]
block_id <- has_data4[,2]
day_id   <- has_data4[,4]

# make table with unique block, site and time combinations
dat <- tibble::as_tibble(cbind(site_id, block_id, time_id))

# remove t=1 as they will be defined outside of this jags loop
dat <- filter(dat, !time_id==1)

uniq_site_block_time <- dplyr::distinct(dat)

# making vectors with indexes for site, block and time, with full time span
uniq_site1_t1 <- dplyr::pull(uniq_site_block_time,1)
uniq_block1_t1 <- dplyr::pull(uniq_site_block_time,2) 
uniq_time1_t1 <- dplyr::pull(uniq_site_block_time,3) 

# remove the last primary occasion, as we would need to have t+1 in the loop in the jags model
uniq_site_block_time_t0 <- dplyr::filter(uniq_site_block_time, !time_id==max(uniq_site_block_time$time_id))

# making vectors with indexes for site, block and time, without last primary occasion
uniq_site1_t0 <- dplyr::pull(uniq_site_block_time_t0,1)
uniq_block1_t0 <- dplyr::pull(uniq_site_block_time_t0,2) 
uniq_time1_t0 <- dplyr::pull(uniq_site_block_time_t0,3) 

# making indexes for uniqe block and time combinations
dat2 <- tibble::as_tibble(cbind(block_id, time_id))

# remove t=1 as it will defined outside this jagsloop
dat2 <- filter(dat2, !time_id==1)

uniq_block_time_t1 <- dplyr::distinct(dat2)

uniq_block2_t1 <- dplyr::pull(uniq_block_time_t1,1) 
uniq_time2_t1 <- dplyr::pull(uniq_block_time_t1,2) 

# remove the last primary occasion, as we would need to have t+1 in the loop in the jags model
uniq_block_time_t0 <- dplyr::filter(uniq_block_time_t1, !time_id==max(uniq_block_time_t1$time_id))

uniq_block2_t0 <- dplyr::pull(uniq_block_time_t0,1) 
uniq_time2_t0 <- dplyr::pull(uniq_block_time_t0,2) 

# make long data
#y_long <- yb[!is.na(yb)]

# add t-1 seaons where they are missing
#y_long <- c(y_long,rep(NA, times=nrow(row_add)))

y_long <- as.vector(has_data4$y_long)

nvisits = length(has_data4)

#load cov
load("season.rda") 
season <- season
season[season==1] <- 2
season[season==0] <- 1

filter(has_data4, dim1==1 & dim2==1 & dim3==39)
filter(has_data4, dim1==1 & dim2==1 & dim3==40)
filter(has_data4, dim1==1 & dim2==1 & dim3==41)

filter(has_data4, dim1==1 & dim2==1 & dim3==174)
filter(has_data4, dim1==1 & dim2==1 & dim3==173)

dplyr::filter(uniq_site_block_time_t0, site_id==1 & block_id==1 & time_id==39)
dplyr::filter(uniq_site_block_time_t0, time_id==1)

dplyr::filter(uniq_block_time_t0, time_id==1)

filter(has_data4, dim1==1 & dim2==1 & dim3==173)

# give data
data <-list(nseason = dim(yb)[3], nblock = dim(yb)[2], nsite = dim(yb)[1], nsurvey = dim(yb)[4], 
            nout=4, y = y_long, season = season, e = 0.0001, nvisits=nvisits,
            time_id=time_id, block_id=block_id, site_id=site_id, day_id=day_id, yb=yb,
            uniq_site1_t0=uniq_site1_t0, uniq_block1_t0=uniq_block1_t0, uniq_time1_t0=uniq_time1_t0,
            uniq_site1_t1=uniq_site1_t1, uniq_block1_t1=uniq_block1_t1, uniq_time1_t1=uniq_time1_t1,
            uniq_block2_t0=uniq_block2_t0, uniq_time2_t0=uniq_time2_t0,
            uniq_block2_t1=uniq_block2_t1, uniq_time2_t1=uniq_time2_t1)

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
            'Chi2Open_B', 'Chi2repOpen_B','Chi2ratioOpen_B', 'Chi2Closed_B', 'Chi2repClosed_B', 'Chi2ratioClosed_B')

#,'FTClosed_B', 'FTrepClosed_B', 'FTratioClosed_B','FTClosed_A', 'FTrepClosed_A', 'FTratioClosed_A'

# MCMC settings
ni <- 10000   ;   nt <- 20   ;   nb <- 2500   ;   nc <- 4    ;   na <- 2500

# run model in jags
setwd("../")

#ptm2 <- proc.time()
out.gof <- jags(data, inits=inits, parameters.to.save=params, "mod_spatial_dynocc_bayp.txt", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)
#proc.time() - ptm2                

# Save model
setwd("./model_output")
save(out.gof, file="outgof_pri3_ni10k.rda")
#~ End of script