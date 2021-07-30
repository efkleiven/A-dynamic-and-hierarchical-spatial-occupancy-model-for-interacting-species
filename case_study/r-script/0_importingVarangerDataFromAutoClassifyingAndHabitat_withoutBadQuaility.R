# read libraries
library(lubridate)
library(dplyr)
library(tidyr)

# set working directory
wd <- "Z:/bilder_smagnagerfotobokser/Data/automatic_classification/MLWIC_classification"
setwd(wd)

# look at files in the directory
dir() 

# import files from Komag based on file position in the folder
filenames <- c("classifications_komagdalen_2016_2020_03_02.csv",
               "classifications_komagdalen_2017_2020_03_02.csv",
               "classifications_komagdalen_2018_2020_03_02.csv",
               "classifications_komagdalen_2019_2020_03_02.csv",
               "classifications_vestre_jakobselv_2019_2020_03_02.csv")

ctdata <- list()
for(i in 1:length(filenames)){
ctdata[[i]] <- read.csv(filenames[i])}

# check that impoting was fine
head(ctdata[[1]])
str(ctdata)

df <- do.call(rbind, ctdata) # merge list objects to one df

# add columns with site name and date 

datesplit <- function(filename)  # retrieve date and station from file name
{
  filename2 <- as.character(filename)
  string <- strsplit(filename2,"/")        # split filename at _
  
  station <-string[[1]][10]  # recollect the filename
  
  string2 <- strsplit(filename2,"_")
  date <- string2[[1]][10] # pick the date

  # should add month , day and week!
  
  outdf <- c(station, date) # collect the thinks you want to save
  return(outdf)
  } # end function


datesite <- as.data.frame(t(sapply(df$fileName,datesplit))) # run function and store as df
names(datesite) <- c("site", "date")

datesite$date <- as.Date(datesite$date)
datesite$day <- format(datesite$date, format = "%d")
datesite$month <- format(datesite$date, format="%m")
datesite$year <- format(datesite$date, format="%Y")
datesite$julian <- julian(datesite$date, origin=min(datesite$date))
datesite$week <- week(datesite$date)
datesite$yearweek<-format(datesite$date, format="%Y-%W")

df2 <- cbind(df,datesite) # merge date and site df with initial df
head(df2) # check that its ok

#0=bad quality, 1=empty, 2=bird, 3=vole, 4=least_weasel, 5=lemming, 6=shrew, 7=stoat, 

# select images classified with more than 90% certainty
df2$answer <- ifelse(df2$confidence1>0.90,df2$guess1,0)


# set answer to NA if quality is bad. This leads to 2768 more NA's in the final multi-state occupancy matrix
table(df2$answer)

df2$answer[df2$answer==0] <- NA
na_vec <- is.na(df2$answer)

summary(df2) # check number of NA's

#remove rows with NA's
df3 <- df2[!na_vec,]

# check that correct number of rows was removed
nrow(df2)-nrow(df3)


df3 <- df3[,-c(1:2,4:8,10:13)] # remove unuseful columns
str(df3)# check that the df is fine

table(df2$answer)

# add column for vole, lemming, stoat and least weasel
df3$vole <- ifelse(df3$answer==3,1,0)
table(df3$vole)

df3$lemming <- ifelse(df3$answer==5,1,0)
table(df3$lemming)

df3$stoat <- ifelse(df3$answer==7,1,0)
table(df3$stoat)

df3$least_weasel <- ifelse(df3$answer==4,1,0)
table(df3$least_weasel)

#############################################################
# import metadata

setwd("H:/UiT/CameraTrapsForSmallMammals/Data")
metadata <- read.csv("Small mammal camera trap metadata locations.csv", header=T, sep=";")
head(metadata)
names(metadata)[3]<-"site" # change name of site column to be identical to df2

metadata <- select(metadata, c("site","habitat"))

dat <- left_join(df3, metadata, by="site") # add on info from metadata

str(dat)

# subset only snowbed habitat (where stoat is seen)
snow <- filter(dat, habitat=="snowbed")

# sum animals per day
vole_c <- aggregate(vole~date, data=snow,sum)
vole_c$date <- as.Date(vole_c$date)

stoat_c <- aggregate(stoat~date, data=snow,sum)
stoat_c$date <- as.Date(stoat_c$date)

lemming_c <- aggregate(lemming~date, data=snow,sum)
lemming_c$date <- as.Date(lemming_c$date)

least_weasel_c <- aggregate(least_weasel~date, data=snow,sum)
least_weasel_c$date <- as.Date(least_weasel_c$date)

plot(vole~date, data=vole_c, type="l", col="blue", ylim=c(0,400))
lines(lemming~date, data=lemming_c, col="red")
lines(stoat~date, data=stoat_c, col="green")
lines(least_weasel~date, data=least_weasel_c, col="orange")
legend("topleft",legend=c("vole","lemming","stoat","least weasel"),lty=1, lwd=2,
       col=c("blue","red","green","orange"), cex=0.5)

# daily occpancy
vole_site <- aggregate(vole~date+site, data=snow, max)
vole_occ <- aggregate(vole~date, data=vole_site, mean)

stoat_site <- aggregate(stoat~date+site, data=snow,max)
stoat_occ <- aggregate(stoat~date, data=stoat_site,mean)

lemming_site <- aggregate(lemming~date+site, data=snow,max)
lemming_occ <- aggregate(lemming~date, data=lemming_site,mean)

least_weasel_site <- aggregate(least_weasel~date+site, data=snow,max)
least_weasel_occ <- aggregate(least_weasel~date, data=least_weasel_site,mean)

plot(vole~date, data=vole_occ, type="l", col="blue")
lines(lemming~date, data=lemming_occ, col="red")
lines(stoat~date, data=stoat_occ, col="green")
lines(least_weasel~date, data=least_weasel_occ, col="orange")
legend("topleft",legend=c("vole","lemming","stoat","least weasel"),lty=1, lwd=2,
       col=c("blue","red","green","orange"), cex=0.5)



##############################################
# make occupancy tables for vole and stoat   #
##############################################

occ_vole <- spread(vole_site, site, vole)
occ_stoat <- spread(stoat_site, site, stoat)
occ_least_weasel <- spread(least_weasel_site, site, least_weasel)
occ_lemming <- spread(lemming_site, site, lemming)

occ_vole2 <- array(NA, dim=c(7,203,46))
occ_stoat2 <- array(NA, dim=c(7,203,46))
occ_least_weasel2 <- array(NA, dim=c(7,203,46))
occ_lemming2 <- array(NA, dim=c(7,203,46))

for(i in 1:46){
  occ_vole2[,,i] <- matrix(unlist(split(occ_vole[1:1421,i+1], ceiling(seq_along(occ_vole[1:1421,i+1])/7))), nrow = 7, byrow = F)
  occ_stoat2[,,i] <- matrix(unlist(split(occ_stoat[1:1421,i+1], ceiling(seq_along(occ_stoat[1:1421,i+1])/7))), nrow = 7, byrow = F)
  occ_least_weasel2[,,i] <- matrix(unlist(split(occ_least_weasel[1:1421,i+1], ceiling(seq_along(occ_least_weasel[1:1421,i+1])/7))), nrow = 7, byrow = F)
  occ_lemming2[,,i] <- matrix(unlist(split(occ_lemming[1:1421,i+1], ceiling(seq_along(occ_lemming[1:1421,i+1])/7))), nrow = 7, byrow = F)
  }

str(occ_vole2)

table(occ_stoat2)
table(occ_lemming2)

occ_mustela <- occ_stoat2+occ_least_weasel2 
occ_rodent <- occ_vole2+occ_lemming2 


table(occ_mustela)
table(occ_rodent)

# save tables
setwd("H:/UiT/CameraTrapsForSmallMammals/Data/summary")
#write.table(occ_vole2, file="occ_snow_vole.txt", sep="\t", quote=F, row.names=F)
#write.table(occ_stoat2, file="occ_snow_stoat.txt", sep="\t", quote=F, row.names=F)
#write.table(occ_mustela, file="occ_snow_mustela.txt", sep="\t", quote=F, row.names=F)

####################################################################################
#### make multi scale occupancy table for mustela and rodent 
occ_rodent[occ_rodent>1] <- 1
occ_mustela[occ_mustela==1] <- 2

occ_ko3 <- occ_mustela+occ_rodent 
table(occ_ko3)

# readjust states to 1:4 instead of 0:3
occ_ko3[occ_ko3==3]<-4; occ_ko3[occ_ko3==2]<-3; occ_ko3[occ_ko3==1]<-2; occ_ko3[occ_ko3==0]<-1

# add block dim
occm_ko3 <- array(NA, dim=c(6,8,203,7))

occ_ko3 <- aperm(occ_ko3, c(3,2,1))

occm_ko3[1:5,1,,] <- occ_ko3[1:5,,] 
occm_ko3[1:5,2,,] <- occ_ko3[6:10,,] 
occm_ko3[,3,,] <- occ_ko3[11:16,,] 
occm_ko3[,4,,] <- occ_ko3[17:22,,]
occm_ko3[,5,,] <- occ_ko3[23:28,,] 
occm_ko3[,6,,] <- occ_ko3[29:34,,] 
occm_ko3[,7,,] <- occ_ko3[35:40,,] 
occm_ko3[,8,,] <- occ_ko3[41:46,,]

table(occ_ko3)
table(occm_ko3)

# save 
setwd("H:/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models/hidden_block_ko_fromautoclass/ko_vj/data")
save(occm_ko3, file="occm_mustela_rodent_var_snowbed_rmBQ.rda")

######################################################################################
### make season as cov ############
##############################

#vole_occ2$Date <- format(vole_occ2$date, format="%m-%d")

#season <- ifelse(vole_occ2$date <= as.Date("2015-11-01"), 0, 
#           ifelse(vole_occ2$date <= as.Date("2016-07-01"), 1,
#             ifelse(vole_occ2$date <= as.Date("2016-11-01"), 0,
#               ifelse(vole_occ2$date <= as.Date("2017-07-01"), 1,
#                 ifelse(vole_occ2$date <= as.Date("2017-11-01"), 0,
#                   ifelse(vole_occ2$date <= as.Date("2018-07-01"), 1,
#                     ifelse(vole_occ2$date <= as.Date("2018-11-01"), 0,
#                       ifelse(vole_occ2$date <= as.Date("2019-07-01"), 1,0))))))))

#season[is.na(season)]<-1
#length(season)
#save(season, file="season.rda")

###################################3
occ_vole$date

c <- occ_vole$date[seq(1,length(occ_vole$date), 7)]

season <- ifelse(c <= as.Date("2015-11-01"), 0, 
                 ifelse(c <= as.Date("2016-07-01"), 1,
                        ifelse(c <= as.Date("2016-11-01"), 0,
                               ifelse(c <= as.Date("2017-07-01"), 1,
                                      ifelse(c <= as.Date("2017-11-01"), 0,
                                             ifelse(c <= as.Date("2018-07-01"), 1,
                                                    ifelse(c <= as.Date("2018-11-01"), 0,
                                                           ifelse(c <= as.Date("2019-07-01"), 1,0))))))))

save(season, file="season.rda")




###############################################################################################################3
##############################################
# make plot with same week def as the occupancy table  #
##############################################

occ_vole <- spread(vole_site, site, vole)
occ_stoat <- spread(stoat_site, site, stoat)
occ_least_weasel <- spread(least_weasel_site, site, least_weasel)
occ_lemming <- spread(lemming_site, site, lemming)


# add week (as every 7 days)
n=7
occ_vole$week <- c(rep(1:(nrow(occ_vole)/n), each=n),rep(204,times=5))
occ_vole$date <- NULL 

occ_stoat$week <- c(rep(1:(nrow(occ_stoat)/n), each=n),rep(204,times=5))
occ_stoat$date <- NULL 

occ_least_weasel$week <- c(rep(1:(nrow(occ_least_weasel)/n), each=n),rep(204,times=5))
occ_least_weasel$date <- NULL 

occ_lemming$week <- c(rep(1:(nrow(occ_lemming)/n), each=n),rep(204,times=5))
occ_lemming$date <- NULL 

occ_vole2 <- aggregate(x=occ_vole, by=list(occ_vole$week), max) # max over every 7 day period for every site
occ_stoat2 <- aggregate(x=occ_stoat, by=list(occ_stoat$week), max) # max over every 7 day period for every site
occ_least_weasel2 <- aggregate(x=occ_least_weasel, by=list(occ_least_weasel$week), max) # max over every 7 day period for every site
occ_lemming2 <- aggregate(x=occ_lemming, by=list(occ_lemming$week), max) # max over every 7 day period for every site

# mean over all columns
str(occ_vole)

#vole
vole<- c()
vole$week <- occ_vole2$week

occ_vole2$week <- NULL
occ_vole2$Group.1 <- NULL

vole$occ <- rowSums(occ_vole2, na.rm=T)

#stoat
stoat <- c()
stoat$week <- occ_stoat2$week

occ_stoat2$week <- NULL
occ_stoat2$Group.1 <- NULL

stoat$occ <- rowSums(occ_stoat2, na.rm=T)

# least weasel
least_weasel<- c()
least_weasel$week <- occ_least_weasel2$week

occ_least_weasel2$week <- NULL
occ_least_weasel2$Group.1 <- NULL

least_weasel$occ <- rowSums(occ_least_weasel2, na.rm=T)

# lemming
lemming<- c()
lemming$week <- occ_lemming2$week

occ_lemming2$week <- NULL
occ_lemming2$Group.1 <- NULL

lemming$occ <- rowSums(occ_lemming2, na.rm=T)

# small rodent
sr_occ <- occ_lemming2+occ_vole2
sr_occ[sr_occ>1]<-1

sr <- c()
sr$week <- c(1:204)
sr$occ <- rowSums(sr_occ, na.rm=T)
sr$mean <- rowMeans(sr_occ, na.rm=T)

# small mustelids
mu_occ <- occ_stoat2+occ_least_weasel2
mu_occ[mu_occ>1]<-1

mu <- c()
mu$week <- c(1:204)
mu$occ <- rowSums(mu_occ, na.rm=T)
mu$mean <- rowMeans(mu_occ, na.rm=T)

table(sr$occ)
table(as.matrix(mu_occ))
str(mu_occ)
###############
# plot
setwd("H:/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models/hidden_block_ko_fromautoclass/ko_vj/plot")

# make plot
png(filename="occ_trend_VA_snowbed.png", width=800, height=400)

plot(NULL, type="n", xlim=c(1,204), ylim=c(0,50),
     xlab="week", ylab="nsite", axes=F, cex.lab=1.5)
axis(side=1, seq(0,204, by=20), cex.axis=1.4)
axis(side=2, seq(0,40, by=10), cex.axis=1.4)
lines(occ ~ week, data=sr, col="black", lwd=3, lty=1)
lines(occ ~ week, data=mu, col="gray",lwd=3, lty=1)
#lines(occ ~ week, data=least_weasel, col="red", lwd=3)
#lines(mustela ~ date, data=plotmustela, col="black")
#lines(occ ~ week, data=lemming, col="brown", lwd=3)
#lines(occ ~ week, data=vole, col="brown")

legend("topleft",legend=c("Small Rodent","Small Mustelid"),lty=1, lwd=3,
       col=c("black","gray"), cex=1)

dev.off()

# make plot of mean occupancy for active cameras
png(filename="occ_trend_VA_snowbed_mean.png", width=800, height=400)

plot(NULL, type="n", xlim=c(1,204), ylim=c(0,1),
     xlab="week", ylab="occupancy", axes=F, cex.lab=1.4)
axis(side=1, seq(0,204, by=20), cex.axis=1.4)
axis(side=2, seq(0,1, by=0.2), cex.axis=1.4)
lines(mean ~ week, data=sr, col="black", lwd=3, lty=1)
lines(mean ~ week, data=mu, col="gray",lwd=3, lty=1)
#lines(occ ~ week, data=least_weasel, col="red", lwd=3)
#lines(mustela ~ date, data=plotmustela, col="black")
#lines(occ ~ week, data=lemming, col="brown", lwd=3)
#lines(occ ~ week, data=vole, col="brown")

legend("top",legend=c("Small Rodent","Small Mustelid"),lty=1, lwd=3,
       col=c("black","gray"), cex=1.2)

dev.off()
