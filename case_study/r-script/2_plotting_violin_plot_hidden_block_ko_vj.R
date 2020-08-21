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

getwd()


