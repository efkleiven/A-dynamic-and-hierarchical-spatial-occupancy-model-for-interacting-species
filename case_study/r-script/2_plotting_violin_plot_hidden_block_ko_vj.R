########################################################################################################################################
##                             Plotting the model output of  the Spatial Dynamic two-species occupancy model - Mac                    ##
##                                  by EFK and FB                                                                                     ##
########################################################################################################################################

# empty environment
rm(list=ls())

# Call jags(and other packages)
library(jagsUI)
library(ggplot2)

# set working directory
setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models")    # for norpec-server
#setwd("H:/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models")  # for laptop

#load model
setwd("./hidden_block_ko_fromautoclass/ko_vj/model_output")
dir()

load("ko_s140-203_250k.rda")
load("va_s140_203_mustela_ni250k.rda")
load("ko_sdet_s161_203_ni250k.rda")
load("ko_mustela_sdet_s161_203_ni250k.rda")
load("ko_mustela_sdet_s155_202_ni250k.rda")
load("ko_mustela_sdet_s1_203_ni250k.rda")
load("va_snowbed_mustela_rodent_sdet_s1_203_ni250k.rda")
load("va_mustela_rodent_seasall_s3_200_ni250k.rda")

####################################################################################################################################
### violin plot with ggplot2

########################
##  ko_sdet_s161_203  ##
########################

dat <- data.frame(sims=unlist(ko_sdet_s161_203_ni250k$sims.list[c(1:8,10:17)]),
                  par=c(rep("\u03B3[A]", times=40000),rep("\u03B3[B]", times=40000),rep("\u03B3[AB]", times=40000),rep("\u03B3[BA]", times=40000),
                        rep("\u03B5[A]", times=40000),rep("\u03B5[B]", times=40000),rep("\u03B5[AB]", times=40000),rep("\u03B5[BA]", times=40000),
                        rep("\u0393[A]", times=40000),rep("\u0393[B]", times=40000),rep("\u0393[AB]", times=40000),rep("\u0393[BA]", times=40000),
                        rep("\u0395[A]", times=40000),rep("\u0395[B]", times=40000),rep("\u0395[AB]", times=40000),rep("\u0395[BA]", times=40000)))

# make a seperate array with estimates from all model runs for all gammas, epsilons and ps
pos <-c("\u03B3[A]", "\u03B3[B]", "\u03B3[AB]", "\u03B3[BA]",
        "\u03B5[A]", "\u03B5[B]", "\u03B5[AB]", "\u03B5[BA]",
        "\u0393[A]", "\u0393[B]", "\u0393[AB]", "\u0393[BA]",
        "\u0395[A]", "\u0395[B]", "\u0395[AB]", "\u0395[BA]")  # spesify what order to plot the simulations

dat1 <- dat[1:320000,]
pos1 <- pos[1:8]

ggplot(data=dat1, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos1)+
  labs(y="", x="")+
  ggtitle("va_sdet_s161-203")

setwd("../plot")

ggsave(
  "modperf_va_sdet_site_s161_203.png"
)


dat2 <- dat[320001:640000,]
pos2 <- pos[9:16]

ggplot(data=dat2, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos2)+
  labs(y="", x="")+
  ggtitle("va_sdet_block_s161-203")

ggsave(
  "modperf_va_sdet_block_s161_203.png"
)
##############################################################################

##################################
##   ko_mustela_sdet_s161_203   ##
##################################

dat <- data.frame(sims=unlist(ko_mustela_sdet_s161_203_ni250k$sims.list[c(1:8,10:17)]),
                  par=c(rep("\u03B3[A]", times=40000),rep("\u03B3[B]", times=40000),rep("\u03B3[AB]", times=40000),rep("\u03B3[BA]", times=40000),
                        rep("\u03B5[A]", times=40000),rep("\u03B5[B]", times=40000),rep("\u03B5[AB]", times=40000),rep("\u03B5[BA]", times=40000),
                        rep("\u0393[A]", times=40000),rep("\u0393[B]", times=40000),rep("\u0393[AB]", times=40000),rep("\u0393[BA]", times=40000),
                        rep("\u0395[A]", times=40000),rep("\u0395[B]", times=40000),rep("\u0395[AB]", times=40000),rep("\u0395[BA]", times=40000))
)

# make a seperate array with estimates from all model runs for all gammas, epsilons and ps
pos <-c("\u03B3[A]", "\u03B3[B]", "\u03B3[AB]", "\u03B3[BA]",
        "\u03B5[A]", "\u03B5[B]", "\u03B5[AB]", "\u03B5[BA]",
        "\u0393[A]", "\u0393[B]", "\u0393[AB]", "\u0393[BA]",
        "\u0395[A]", "\u0395[B]", "\u0395[AB]", "\u0395[BA]" )  # spesify what order to plot the simulations

dat1 <- dat[1:320000,]
pos1 <- pos[1:8]

ggplot(data=dat1, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos1)+
  labs(y="", x="")+
  ggtitle("va_mustela_sdet_site_s161-203")

ggsave(
  "modperf_va_mustela_sdet_site_s161_203.png"
)


dat2 <- dat[320001:640000,]
pos2 <- pos[9:16]

ggplot(data=dat2, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos2)+
  labs(y="", x="")+
  ggtitle("va_mustela_sdet_block_s161-203")

ggsave(
  "modperf_va_mustela_sdet_block_s161_203.png"
)

###########################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   ko_mustela_sdet_s1_203        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
str(ko_mustela_sdet_s1_203_ni250k)
dat <- data.frame(sims=unlist(ko_mustela_sdet_s1_203_ni250k$sims.list[c(1:8,10:17)]),
                  par=c(rep("\u03B3[A]", times=40000),rep("\u03B3[B]", times=40000),rep("\u03B3[AB]", times=40000),rep("\u03B3[BA]", times=40000),
                        rep("\u03B5[A]", times=40000),rep("\u03B5[B]", times=40000),rep("\u03B5[AB]", times=40000),rep("\u03B5[BA]", times=40000),
                        rep("\u0393[A]", times=40000),rep("\u0393[B]", times=40000),rep("\u0393[AB]", times=40000),rep("\u0393[BA]", times=40000),
                        rep("\u0395[A]", times=40000),rep("\u0395[B]", times=40000),rep("\u0395[AB]", times=40000),rep("\u0395[BA]", times=40000)))

# make a seperate array with estimates from all model runs for all gammas, epsilons and ps
pos <-c("\u03B3[A]", "\u03B3[B]", "\u03B3[AB]", "\u03B3[BA]",
        "\u03B5[A]", "\u03B5[B]", "\u03B5[AB]", "\u03B5[BA]",
        "\u0393[A]", "\u0393[B]", "\u0393[AB]", "\u0393[BA]",
        "\u0395[A]", "\u0395[B]", "\u0395[AB]", "\u0395[BA]")  # spesify what order to plot the simulations

dat1 <- dat[1:320000,]
pos1 <- pos[1:8]

ggplot(data=dat1, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos1)+
  labs(y="", x="")+
  ggtitle("va_snowbed_mustela_sdet_block_s1-203")

setwd("../plot")

ggsave(
  "modperf_va_snowbed_mustela_sdet_site_s1_203.png"
)


dat2 <- dat[320001:640000,]
pos2 <- pos[9:16]

ggplot(data=dat2, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos2)+
  labs(y="", x="")+
  ggtitle("va_snowbed_mustela_sdet_block_s1-203")

ggsave(
  "modperf_va_snowbed_mustela_sdet_block_s1_203.png"
)
##############################################################################

###########################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   va_snowbed_mustela_vole_sdet_s1_203        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat <- data.frame(sims=unlist(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$sims.list[c(1:8,10:17)]),
                  par=c(rep("\u03B3[A]", times=40000),rep("\u03B3[B]", times=40000),rep("\u03B3[AB]", times=40000),rep("\u03B3[BA]", times=40000),
                        rep("\u03B5[A]", times=40000),rep("\u03B5[B]", times=40000),rep("\u03B5[AB]", times=40000),rep("\u03B5[BA]", times=40000),
                        rep("\u0393[A]", times=40000),rep("\u0393[B]", times=40000),rep("\u0393[AB]", times=40000),rep("\u0393[BA]", times=40000),
                        rep("\u0395[A]", times=40000),rep("\u0395[B]", times=40000),rep("\u0395[AB]", times=40000),rep("\u0395[BA]", times=40000)))

# make a seperate array with estimates from all model runs for all gammas, epsilons and ps
pos <-c("\u03B3[A]", "\u03B3[B]", "\u03B3[AB]", "\u03B3[BA]",
        "\u03B5[A]", "\u03B5[B]", "\u03B5[AB]", "\u03B5[BA]",
        "\u0393[A]", "\u0393[B]", "\u0393[AB]", "\u0393[BA]",
        "\u0395[A]", "\u0395[B]", "\u0395[AB]", "\u0395[BA]")  # spesify what order to plot the simulations

dat1 <- dat[1:320000,]
pos1 <- pos[1:8]

ggplot(data=dat1, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos1)+
  labs(y="", x="")+
  ggtitle("va_snowbed_mustela_rodent_sdet_site_s1-203")

setwd("../plot")

ggsave(
  "modperf_va_snowbed_mustela_rodent_sdet_site_s1_203.png"
)


dat2 <- dat[320001:640000,]
pos2 <- pos[9:16]

ggplot(data=dat2, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos2)+
  labs(y="", x="")+
  ggtitle("va_snowbed_mustela_rodent_sdet_block_s1-203")

ggsave(
  "modperf_va_snowbed_mustela_rodent_sdet_block_s1_203.png"
)
##############################################################################

###########################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   va_snowbed_mustela_vole_seasall_s1_203      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

vec <- c()
x <- 40000

vec[1:x] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[1]][,1]
vec[(x+1):(2*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[1]][,15]
vec[(2*x+1):(3*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[2]][,1]
vec[(3*x+1):(4*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[2]][,15]
vec[(4*x+1):(5*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[3]][,1]
vec[(5*x+1):(6*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[3]][,15]
vec[(6*x+1):(7*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[4]][,1]
vec[(7*x+1):(8*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[4]][,15]
vec[(8*x+1):(9*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[5]][,1]
vec[(9*x+1):(10*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[5]][,15]
vec[(10*x+1):(11*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[6]][,1]
vec[(11*x+1):(12*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[6]][,15]
vec[(12*x+1):(13*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[7]][,1]
vec[(13*x+1):(14*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[7]][,15]
vec[(14*x+1):(15*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[8]][,1]
vec[(15*x+1):(16*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[8]][,15]
vec[(16*x+1):(17*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[9]][,1]
vec[(17*x+1):(18*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[9]][,15]
vec[(18*x+1):(19*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[10]][,1]
vec[(19*x+1):(20*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[10]][,15]
vec[(20*x+1):(21*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[11]][,1]
vec[(21*x+1):(22*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[11]][,15]
vec[(22*x+1):(23*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[12]][,1]
vec[(23*x+1):(24*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[12]][,15]
vec[(24*x+1):(25*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[13]][,1]
vec[(25*x+1):(26*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[13]][,15]
vec[(26*x+1):(27*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[14]][,1]
vec[(27*x+1):(28*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[14]][,15]
vec[(28*x+1):(29*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[15]][,1]
vec[(29*x+1):(30*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[15]][,15]
vec[(30*x+1):(31*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[16]][,1]
vec[(31*x+1):(32*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:17)][[16]][,15]
vec[(32*x+1):(33*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:19)][[17]][,1]
vec[(33*x+1):(34*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:19)][[17]][,15]
vec[(34*x+1):(35*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:19)][[18]][,1]
vec[(35*x+1):(36*x)] <- va_mustela_rodent_seasall_s3_200_ni250k$sims.list[c(1:8,10:19)][[18]][,15]



dat <- data.frame(sims=vec,
                  par=c(rep("\u03B3[A_s]", times=40000), rep("\u03B3[A_w]", times=40000),
                        rep("\u03B3[B_s]", times=40000), rep("\u03B3[B_w]", times=40000),
                        rep("\u03B3[AB_s]", times=40000), rep("\u03B3[AB_w]", times=40000),
                        rep("\u03B3[BA_s]", times=40000), rep("\u03B3[BA_w]", times=40000),
                        rep("\u03B5[A_s]", times=40000), rep("\u03B5[A_w]", times=40000),
                        rep("\u03B5[B_s]", times=40000),rep("\u03B5[B_w]", times=40000),
                        rep("\u03B5[AB_s]", times=40000),rep("\u03B5[AB_w]", times=40000),
                        rep("\u03B5[BA_s]", times=40000),rep("\u03B5[BA_w]", times=40000),
                        rep("\u0393[A_s]", times=40000), rep("\u0393[A_w]", times=40000),
                        rep("\u0393[B_s]", times=40000),rep("\u0393[B_w]", times=40000),
                        rep("\u0393[AB_s]", times=40000),rep("\u0393[AB_w]", times=40000),
                        rep("\u0393[BA_s]", times=40000),rep("\u0393[BA_w]", times=40000),
                        rep("\u0395[A_s]", times=40000),rep("\u0395[A_w]", times=40000),
                        rep("\u0395[B_s]", times=40000),rep("\u0395[B_w]", times=40000),
                        rep("\u0395[AB_s]", times=40000),rep("\u0395[AB_w]", times=40000),
                        rep("\u0395[BA_s]", times=40000),rep("\u0395[BA_w]", times=40000),
                        rep("pA_s", times=40000),rep("pA_w", times=40000),
                        rep("pB_s", times=40000),rep("pB_w", times=40000)))

# make a seperate array with estimates from all model runs for all gammas, epsilons and ps
pos <-c("\u03B3[A_s]", "\u03B3[A_w]", "\u03B3[B_s]", "\u03B3[B_w]", "\u03B3[AB_s]", "\u03B3[AB_w]", "\u03B3[BA_s]", "\u03B3[BA_w]",
        "\u03B5[A_s]", "\u03B5[A_w]", "\u03B5[B_s]", "\u03B5[B_w]", "\u03B5[AB_s]", "\u03B5[AB_w]", "\u03B5[BA_s]", "\u03B5[BA_w]",
        "\u0393[A_s]", "\u0393[A_w]", "\u0393[B_s]", "\u0393[B_w]", "\u0393[AB_s]", "\u0393[AB_w]", "\u0393[BA_s]", "\u0393[BA_w]",
        "\u0395[A_s]", "\u0395[A_w]", "\u0395[B_s]", "\u0395[B_w]", "\u0395[AB_s]", "\u0395[AB_w]", "\u0395[BA_s]", "\u0395[BA_w]",
        "pA_s","pA_w","pB_s","pB_w")  # spesify what order to plot the simulations


dat1 <- dat[1:320000,]
pos1 <- pos[1:8]

ggplot(data=dat1, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos1)+
  labs(y="", x="")+
  ggtitle("va_seasall_site1_s1-203")

ggsave(
  "modperf_va_seasall_site1_s1_203.png"
)

dat2 <- dat[320001:640000,]
pos2 <- pos[9:16]

ggplot(data=dat2, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos2)+
  labs(y="", x="")+
  ggtitle("va_seasall_site2_s1-203")

ggsave(
  "modperf_va_seasall_site2_s1_203.png"
)


###################################################
dat4 <- dat[640001:960000,]
pos4 <- pos[17:24]

ggplot(data=dat4, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos4)+
  labs(y="", x="")+
  ggtitle("va_seasall_block1_s1-203")

ggsave(
  "modperf_va_seasall_block1_s1_203.png"
)

##########################################################
# make subplot to remove overly wide parameters so that it's easier to see the distribution of the other parameters
dat6 <- dat[c(640001:680000,720001:760000,800001:960000),]
pos6 <- pos[c(17,19,21:24)]

ggplot(data=dat6, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos6)+
  labs(y="", x="")+
  ggtitle("va_seasall_block1_s1-203")

ggsave(
  "modperf_va_seasall_block_subplot_1_s1_203.png"
)
###################################################

dat5 <- dat[960001:1280000),]
pos5 <- pos[25:32]

ggplot(data=dat5, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos5)+
  labs(y="", x="")+
  ggtitle("va_seasall_block2_s1-203")

ggsave("modperf_va_seasall_block2_s1_203.png")

########################################################




dat3 <- dat[1280001:1440000,]
pos3 <- pos[33:36]

ggplot(data=dat3, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  scale_x_discrete(limits = pos3)+
  labs(y="", x="")+
  ggtitle("va_seasall_p_s1-203")

ggsave("modperf_va_seasall_p_s1-203.png")
