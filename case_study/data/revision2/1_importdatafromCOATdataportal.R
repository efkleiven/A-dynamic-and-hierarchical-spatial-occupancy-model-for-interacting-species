#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     require libraries 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     Import data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/ekl013/OneDrive - UiT Office 365/GitProjects/MustelidsAndRodents-/data/coatdata")
dir()
filenames <- dir()
filenames <- filenames[-c(1:3)] # remove not relevant filenames (because I have some folder in the dir now)

ctdata <- list()

for(i in 1:length(filenames)){
  ctdata[[i]] <- read.table(filenames[i], header=T, sep=";")}

# merge list objects to one df
df <- do.call(rbind, ctdata) 
str(df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          reshape data    ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# remove some columns that make problems in the pivot_wider
df2 <- select(df, !c(v_confidence_automatic, v_comment))

# structure dates
df2$t_date <- as.Date(df2$t_date)

# remove early and late periods
# use only data until 2021 in JABES manus
df2 <-filter(df2, sc_type_of_sites_ecological =="snowbed")
df2 <-filter(df2, t_date < "2021-07-01" & t_date > "2015-08-31" )

# some of the cameras did not get checked until september, as they were covered in snow during the july fieldwork.

# use pivot_wider to go from long to wide format (spread species as columns instead of rows)
df3 <- pivot_wider(df2, names_from=v_class_id, values_from = c(v_presence_automatic, v_presence_manual), values_fn = {max})

str(df3)

# add a column for the classification of each image
# use manual classification if it exist and automatic if not

df4 <- mutate(.data = df3, 
              species= case_when(v_presence_manual_bad_quality==1 ~ "bad_quality",
                                 v_presence_manual_empty==1 ~ "empty",
                                 v_presence_manual_cricetidae==1 ~ "vole",
                                 v_presence_manual_lem_lem==1 ~ "lem",
                                 v_presence_manual_mus_erm==1 ~ "stoat",
                                 v_presence_manual_mus_niv==1 ~ "least_weasel",
                                 v_presence_manual_sor_sp==1 ~ "shrew",
                                 v_presence_manual_aves==1 ~ "bird",
                                 v_presence_automatic_bad_quality==1 ~ "bad_quality",
                                 v_presence_automatic_empty==1 ~ "empty",
                                 v_presence_automatic_cricetidae==1 ~ "vole",
                                 v_presence_automatic_lem_lem==1 ~ "lem",
                                 v_presence_automatic_mus_erm==1 ~ "stoat",
                                 v_presence_automatic_mus_niv==1 ~ "least_weasel",
                                 v_presence_automatic_sor_sp==1 ~ "shrew",
                                 v_presence_automatic_aves==1 ~ "bird",
                                 TRUE ~ "no"))

# table species to check that it look ok
table(df4$species)
# 18 images seems to have not been annoteted, these are the 18 mink images 

# remove columns that are no longer needed
names(df4)
df5 <- select(df4,c("sn_region","sn_locality","sn_section", "sc_type_of_sites_ecological","sn_site","t_date","t_time","v_image_name", "species"))

# adding separate columns  for interesting species
df6 <- mutate(.data=df5,
              bq=case_when(species=="bad_quality"~1,
                           TRUE ~ 0),
              lemming=case_when(species=="lem" ~ 1,
                                TRUE ~ 0),
              vole=case_when(species=="vole" ~ 1,
                             TRUE ~ 0),
              rodent=case_when(species=="lem" ~ 1,
                               species=="vole" ~ 1,
                               TRUE ~ 0),
              stoat=case_when(species=="stoat" ~ 1,
                              TRUE ~ 0),
              least_weasel=case_when(species=="least_weasel" ~ 1,
                                     TRUE ~ 0),
              mustelid=case_when(species=="stoat" ~ 1,
                                 species=="least_weasel" ~ 1,
                                 TRUE ~ 0),
)

# check that species columns are ok
table(df6$lemming)
table(df6$vole)
table(df6$rodent)
table(df6$stoat)
table(df6$least_weasel)
table(df6$mustelid)

table(df5$species)

# when sites are moved (e.g. because they have been flooded) they are given new names
# manually standardizing names from the same sites

length(unique(df6$sn_site))

df6$sn_site[df6$sn_site=="ko_ga_hu_1b"] <- "ko_ga_hu_1"
df6$sn_site[df6$sn_site=="ko_ga_hu_4b"] <- "ko_ga_hu_4"
df6$sn_site[df6$sn_site=="ko_ga_sn_1b"] <- "ko_ga_sn_1"
df6$sn_site[df6$sn_site=="ko_hu_hu_3b"] <- "ko_hu_hu_3"
df6$sn_site[df6$sn_site=="ko_kj_hu_1b"] <- "ko_kj_hu_1"
df6$sn_site[df6$sn_site=="ko_kj_hu_3b"] <- "ko_kj_hu_3"
df6$sn_site[df6$sn_site=="ko_ga_sn_2b"] <- "ko_ga_sn_2"
df6$sn_site[df6$sn_site=="ko_ga_sn_4b"] <- "ko_ga_sn_4"
df6$sn_site[df6$sn_site=="ko_hu_hu_2b"] <- "ko_hu_hu_2"
df6$sn_site[df6$sn_site=="ko_ry_hu_1b"] <- "ko_ry_hu_1"
df6$sn_site[df6$sn_site=="ko_kj_hu_4b"] <- "ko_kj_hu_4"
# as of 2021 we should have 92 sites (where are the images from rÃ¸yskattfjellet?)
length(unique(df6$sn_site))

# order df6
df6 <- df6[order(df6$sn_locality, df6$sn_section, df6$sn_site, df6$t_date),]

# set species to NA it images is classified as bad quality

df6$lemming <- ifelse(df6$bq==1, NA, df6$lemming)
df6$vole <- ifelse(df6$bq==1, NA, df6$vole)
df6$rodent <- ifelse(df6$bq==1, NA, df6$rodent)
df6$stoat <- ifelse(df6$bq==1, NA, df6$stoat)
df6$least_weasel <- ifelse(df6$bq==1, NA, df6$least_weasel)
df6$mustelid <- ifelse(df6$bq==1, NA, df6$mustelid)

# daily occpancy and sum(for plotting)
# every day that now includes at least one image classified as bad quality will now be set as NA
str(df6)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save df6
#setwd("C:/Users/ekl013/OneDrive - UiT Office 365/GitProjects/MustelidsAndRodents-/data/coatdata/modified")
#save(df6, file="imageannotations.rda")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vole_site <- aggregate(vole~t_date+sn_site, data=df6, max, na.action = na.pass)
lemming_site <- aggregate(lemming~t_date+sn_site, data=df6,max, na.action = na.pass)
rodent_site <- aggregate(rodent~t_date+sn_site, data=df6, max, na.action = na.pass)

stoat_site <- aggregate(stoat~t_date+sn_site, data=df6, max, na.action = na.pass)
least_weasel_site <- aggregate(least_weasel~t_date+sn_site, data=df6, max, na.action = na.pass)
mustela_site <- aggregate(mustelid~t_date+sn_site, data=df6, max, na.action = na.pass)


##############################################
# make occupancy tables for vole and stoat   #
##############################################

# make a column for each site
occ_vole <- spread(vole_site, sn_site, vole)[,-c(1)]
occ_lemming <- spread(lemming_site, sn_site, lemming)[,-c(1)]
occ_stoat <- spread(stoat_site, sn_site, stoat)[,-c(1)]
occ_least_weasel <- spread(least_weasel_site, sn_site, least_weasel)[,-c(1)]

occ_mustela <- spread(mustela_site, sn_site, mustelid)[,-c(1)]
occ_rodent <- spread(rodent_site, sn_site, rodent)[,-c(1)]

# specifying dimensions of the multi-state occupancy table
ndays <- 7
nsite <- dim(occ_vole)[2]
nprocc <- floor(dim(occ_vole)[1]/ndays)
ntime <- ndays*nprocc

# creating empty array with correct dimensions
occ_vole2 <- array(NA, dim=c(ndays,nprocc,nsite))
occ_lemming2 <- array(NA, dim=c(ndays,nprocc,nsite))

occ_stoat2 <- array(NA, dim=c(ndays,nprocc,nsite))
occ_least_weasel2 <- array(NA, dim=c(ndays,nprocc,nsite))

occ_rodent2 <- array(NA, dim=c(ndays,nprocc,nsite))
occ_mustela2 <- array(NA, dim=c(ndays,nprocc,nsite))

# split the observations into primary and secondary occasions (different dimensions in the array) 
for(i in 1:nsite){
  occ_vole2[,,i] <- matrix(unlist(split(occ_vole[1:ntime,i], ceiling(seq_along(occ_vole[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  occ_lemming2[,,i] <- matrix(unlist(split(occ_lemming[1:ntime,i], ceiling(seq_along(occ_lemming[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  
  occ_stoat2[,,i] <- matrix(unlist(split(occ_stoat[1:ntime,i], ceiling(seq_along(occ_stoat[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  occ_least_weasel2[,,i] <- matrix(unlist(split(occ_least_weasel[1:ntime,i], ceiling(seq_along(occ_least_weasel[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  
  occ_mustela2[,,i] <- matrix(unlist(split(occ_mustela[1:ntime,i], ceiling(seq_along(occ_mustela[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  occ_rodent2[,,i] <- matrix(unlist(split(occ_rodent[1:ntime,i], ceiling(seq_along(occ_rodent[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  
}

str(occ_vole2)
table(occ_vole2)
table(occ_rodent2)
table(occ_lemming2)
table(occ_mustela2)

#############################################################
## Make multi state occupancy tables with block structure ###
#############################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For voles and mustelids

# set the mustelids to 2
occ_mustela2[occ_mustela2==1] <- 2
occ_rodent <- occ_lemming2+occ_vole2
occ_rodent[occ_rodent==2]<-1

# making a multi state array for mustelids and rodents
occ_va1 <- occ_mustela2+occ_rodent
str(occ_va1)


# changing states from 0,1,2,3 to 1,2,3,4 (needed for the categorical dist in Jags)
occ_va1[occ_va1==3]<-4
occ_va1[occ_va1==2]<-3
occ_va1[occ_va1==1]<-2
occ_va1[occ_va1==0]<-1

# add block dim
nblocks <- 8

# making empty array to fill in the multi state occupancy matrix
# a problem here is that the sites in Komag have only 11 sites per block, 
# while the sites in VJ have 12 sites per block. Now I create empty sites in komag, which dosen't make sence.
# probably this will also cause problems when we want to include covariates. 

occm_va <- array(NA, dim=c(6,nblocks,nprocc,ndays))

occ_va1 <- aperm(occ_va1, c(3,2,1))
dim(occ_va1)

unique(df6$sn_site)

######
### need to be changes to only include snowbedsites
occm_va[1:5,1,,] <- occ_va1[1:5,,] 
occm_va[1:5,2,,] <- occ_va1[6:10,,] 
occm_va[1:6,3,,] <- occ_va1[11:16,,] 
occm_va[1:6,4,,] <- occ_va1[17:22,,]
occm_va[,5,,] <- occ_va1[23:28,,] 
occm_va[,6,,] <- occ_va1[29:34,,] 
occm_va[,7,,] <- occ_va1[35:40,,] 
occm_va[,8,,] <- occ_va1[41:46,,] 

table(occm_va)
dim(occm_va)

# check number of NA's
summary(occm_va) 

# save 
setwd("C:/Users/ekl013/OneDrive - UiT Office 365/GitProjects/A-dynamic-and-hierarchical-spatial-occupancy-model-for-interacting-species/case_study/data/revision2")
save(occm_va, file="occm_vole_mustelid_snowbed_2016_2021.rda")

###########################
# make season covariate ##
###########################


datvec=seq(as.Date("2015-09-01"), as.Date("2021-07-01"), by="7 day")

seas <- ifelse(datvec < as.Date("2015-11-01"), "s", 
               ifelse(datvec < as.Date("2016-07-01"), "w",
                      ifelse(datvec < as.Date("2016-11-01"), "s",
                             ifelse(datvec < as.Date("2017-07-01"), "w",
                                    ifelse(datvec < as.Date("2017-11-01"), "s",
                                           ifelse(datvec < as.Date("2018-07-01"), "w", 
                                                  ifelse(datvec < as.Date("2018-11-01"), "s",
                                                         ifelse(datvec < as.Date("2019-07-01"), "w",
                                                                ifelse(datvec < as.Date("2019-11-01"), "s",
                                                                       ifelse(datvec < as.Date("2020-07-01"), "w",
                                                                              ifelse(datvec < as.Date("2020-11-01"), "s", "w")))))))))))

setwd("C:/Users/ekl013/OneDrive - UiT Office 365/GitProjects/A-dynamic-and-hierarchical-spatial-occupancy-model-for-interacting-species/case_study/data/revision2")
save(seas, file="seas_2021.rda")
#~ End of script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

