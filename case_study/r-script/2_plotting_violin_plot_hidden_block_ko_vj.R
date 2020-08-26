###############################################################################################################################
##                    Plotting the model output of  the Spatial Dynamic two-species occupancy model                          ##
##                                                                                                                           ##
##                         by EFK and FB                                                                                     ##
###############################################################################################################################

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

 load("va_snowbed_mustela_rodent_sdet_s1_203_ni250k.rda")

#####################################
###   violin plot with ggplot2    ###
#####################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   va_snowbed_mustela_vole_sdet_s1_203     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##############################################################################

dat <- data.frame(sims=unlist(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$sims.list[c(1:8,10:17)]),
                  par=c(rep("gamA", times=40000),rep("gamB", times=40000),rep("gamAB", times=40000),rep("gamBA", times=40000),
                        rep("epsA", times=40000),rep("epsB", times=40000),rep("epsAB", times=40000),rep("epsBA", times=40000),
                        rep("GamA", times=40000),rep("GamB", times=40000),rep("GamAB", times=40000),rep("GamBA", times=40000),
                        rep("EpsA", times=40000),rep("EpsB", times=40000),rep("EpsAB", times=40000),rep("EpsBA", times=40000)))


 dat2 <- data.frame(sims=unlist(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$mean[c(1:8,10:17)]),
                    par=c("gamA", "gamB", "gamAB", "gamBA", "epsA", "epsB", "epsAB", "epsBA",
                          "GamA", "GamB", "GamAB", "GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
 
# make a seperate array with estimates from all model runs for all gammas, epsilons and ps
pos <-c("gamA", "gamB", "gamAB", "gamBA",
        "epsA", "epsB", "epsAB", "epsBA",
        "GamA", "GamB", "GamAB", "GamBA",
        "EpsA", "EpsB", "EpsAB", "EpsBA")  # spesify what order to plot the simulations

ggplot(data=dat, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  theme(axis.text.x = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=20, face="bold"), legend.position = "none")+
  scale_x_discrete(limits=pos, labels=c('gamA' = parse(text = TeX("$\\gamma_{A}$")), 'gamB' = parse(text = TeX("$\\gamma_{B}$")), 
                            'gamAB' = parse(text = TeX("$\\gamma_{AB}$")), 'gamBA' = parse(text = TeX("$\\gamma_{BA}$")),
                            'epsA' = parse(text = TeX("$\\epsilon_{A}$")), 'epsB' = parse(text = TeX("$\\epsilon_{B}$")),
                            'epsAB' = parse(text = TeX("$\\epsilon_{AB}$")), 'epsBA' = parse(text = TeX("$\\epsilon_{BA}$")),
                            'GamA' = parse(text = TeX("$\\Gamma_{A}$")), 'GamB' = parse(text = TeX("$\\Gamma_{B}$")),
                            'GamAB' = parse(text = TeX("$\\Gamma_{AB}$")), 'GamBA' = parse(text = TeX("$\\Gamma_{BA}$")),
                            'EpsA' = parse(text = TeX("$E_{A}$")), 'EpsB' = parse(text = TeX("$E_{B}$")), 
                            'EpsAB' = parse(text = TeX("$E_{AB}$")), 'EpsBA' = parse(text = TeX("$E_{BA}$")) ))

# you should add posterior mean/median to the figure as well

setwd("../plot")
ggsave("modperf_va_snowbed_mustela_rodent_sdet_s1_203_site&block.png", width = 60, height = 20, units="cm")

###########################################################################
## plot estimated detection probabilities


dat3 <- data.frame(sims=c(unlist(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$sims.list[18][[1]][,c(1,15)]),
                          unlist(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$sims.list[19][[1]][,c(1,15)])),
           par=c(rep("pA_S", times=40000), rep("pA_W", times=40000), rep("pB_S", times=40000), rep("pB_W", times=40000)) )

dat4 <- data.frame(sims=c(unlist(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$mean[18][[1]][c(1,15)]),
                          unlist(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$mean[19][[1]][c(1,15)])),
                   par=c("pA_S", "pA_W", "pB_S", "pB_W"))


pos3 <-c("pA_S", "pA_W", "pB_S", "pB_W")  # spesify what order to plot the simulations

ggplot(data=dat3, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_point(data=dat4, color="red", shape="-", size=15)+
  labs(y="", x="")+
  theme(axis.text.x = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=20, face="bold"), legend.position = "none")+
  scale_x_discrete(limits=pos3, labels=c('pA_S' = parse(text = TeX("$p_{A_S}$")), 'pA_W' = parse(text = TeX("$p_{A_W}$")), 
                                        'pB_S' = parse(text = TeX("$p_{B_S}$")), 'pB_W' = parse(text = TeX("$p_{B_W}$"))))


ggsave("modperf_va_snowbed_mustela_rodent_sdet_s1_203_p.png", width = 60, height = 20, units="cm")


###########################################################################
## plot estimated initial occupancy probability

dat5 <- data.frame(sims=c(as.vector(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$sims.list$psi[,,1]),
                          as.vector(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$sims.list$psi[,,2]),
                          as.vector(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$sims.list$psi[,,3])),
                   par=c(rep("psi1", times=320000), rep("psi2", times=320000), rep("psi3", times=320000)),
                  block=c(rep(c(rep("b1", times=40000),rep("b2",times=40000),rep("b3", times=40000),rep("b4", times=40000),rep("b5", times=40000),rep("b6", times=40000),rep("b7", times=40000),rep("b8", times=40000))
                              , times=3)))

dat6 <- data.frame(sims=c(va_snowbed_mustela_rodent_sdet_s1_203_ni250k$mean$psi[,1],
                          va_snowbed_mustela_rodent_sdet_s1_203_ni250k$mean$psi[,2],
                          va_snowbed_mustela_rodent_sdet_s1_203_ni250k$mean$psi[,3]),
                   par=c(rep("psi1", times=8),rep("psi2", times=8),rep("psi3", times=8)),
                   block=c(rep(c("b1","b2","b3","b4","b5","b6","b7","b8"), times=3)) )


pos5 <-c("psi1","psi2","psi3")  # spesify what order to plot the simulations

ggplot(data=dat5, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_point(data=dat6, color="red", shape="-", size=15)+
  labs(y="", x="")+
  facet_wrap(~block, labeller = label_parsed, ncol=4)+
    theme(strip.text = element_text(size=16,face="bold"), axis.text.x = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=20, face="bold"), legend.position = "none")+
  scale_x_discrete(limits=pos5, 
                   labels=c('psi1' = parse(text = TeX("$\\psi_{1}$")), 'psi2' = parse(text = TeX("$\\psi_{2}$")), 'psi3' = parse(text = TeX("$\\psi_{3}$"))))


ggsave("modperf_va_snowbed_mustela_rodent_sdet_s1_203_psi_2.png", width = 60, height = 30, units="cm")

