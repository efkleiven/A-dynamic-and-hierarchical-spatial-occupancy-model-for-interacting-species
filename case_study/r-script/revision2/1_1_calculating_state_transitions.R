dim(yb)

# set wd

setwd("../../data")

dir()

load("occm_vole_mustelid_snowbed_2016_2021.rda")
yb <- occm_va 

# take max from days
yb2 <- apply(yb, c(1,2,3), function(x) max(x,na.rm=T)) 

dim(yb2)

# take max from site within block
yb3 <- apply(yb2, c(2,3), function(x) max(x,na.rm=T))

yb3[,200]

# check how often one observe a transition from 
ntime <- dim(yb)[3]

UU <- array(NA, dim=c(8,ntime))
UA <- array(NA, dim=c(8,ntime))
UB <- array(NA, dim=c(8,ntime))
UAB <- array(NA, dim=c(8,ntime))

AU <- array(NA, dim=c(8,ntime))
AA <- array(NA, dim=c(8,ntime))
AB <- array(NA, dim=c(8,ntime))
AAB <- array(NA, dim=c(8,ntime))

BU <- array(NA, dim=c(8,ntime))
BA <- array(NA, dim=c(8,ntime))
BB <- array(NA, dim=c(8,ntime))
BAB <- array(NA, dim=c(8,ntime))

ABU <- array(NA, dim=c(8,ntime))
ABA <- array(NA, dim=c(8,ntime))
ABB <- array(NA, dim=c(8,ntime))
ABAB <- array(NA, dim=c(8,ntime))

for(b in 1:8){
  for(t in 2:ntime){
    UU[b,t] <- ifelse(yb3[b,t]==1 & yb3[b,t-1]==1, 1, 0)
    UA[b,t] <- ifelse(yb3[b,t]==2 & yb3[b,t-1]==1, 1, 0)
    UB[b,t] <- ifelse(yb3[b,t]==3 & yb3[b,t-1]==1, 1, 0)
    UAB[b,t] <- ifelse( yb3[b,t]==4 & yb3[b,t-1]==1, 1, 0)
  
    AU[b,t] <- ifelse(yb3[b,t]==1 & yb3[b,t-1]==2, 1, 0)
    AA[b,t] <- ifelse(yb3[b,t]==2 & yb3[b,t-1]==2, 1, 0)
    AB[b,t] <- ifelse(yb3[b,t]==3 & yb3[b,t-1]==2, 1, 0)
    AAB[b,t] <- ifelse( yb3[b,t]==4 & yb3[b,t-1]==2, 1, 0)    
  
    BU[b,t] <- ifelse(yb3[b,t]==1 & yb3[b,t-1]==3, 1, 0)
    BA[b,t] <- ifelse(yb3[b,t]==2 & yb3[b,t-1]==3, 1, 0)
    BB[b,t] <- ifelse(yb3[b,t]==3 & yb3[b,t-1]==3, 1, 0)
    BAB[b,t] <- ifelse( yb3[b,t]==4 & yb3[b,t-1]==3, 1, 0)  
    
    ABU[b,t] <- ifelse(yb3[b,t]==1 & yb3[b,t-1]==4, 1, 0)
    ABA[b,t] <- ifelse(yb3[b,t]==2 & yb3[b,t-1]==4, 1, 0)
    ABB[b,t] <- ifelse(yb3[b,t]==3 & yb3[b,t-1]==4, 1, 0)
    ABAB[b,t] <- ifelse( yb3[b,t]==4 & yb3[b,t-1]==4, 1, 0)  
    }}

table(UU)
table(UA)
table(UB)
table(UAB)

table(AU)
table(AA)
table(AB)
table(ABA)

table(BU)
table(BA)
table(BB)
table(BAB)

table(ABU)
table(ABA)
table(ABB)
table(ABAB)
#



dim(yb)


