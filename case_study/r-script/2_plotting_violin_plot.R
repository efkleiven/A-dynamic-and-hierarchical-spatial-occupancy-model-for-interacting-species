####################################################################################################################
##   Plotting the model output from the dynamic occupancy model for interacting species with two spatial scales   ##
##   analyzing the small mammal data from the Varanger peninsula                                                  ##
##                         by EFK and FB                                                                          ##
####################################################################################################################

# empty environment
rm(list=ls())

# Call jags(and other packages)
library(jagsUI)
library(ggplot2)
library(latex2exp)

# set working directory
setwd("./case_study/model_output")

#load model

dir()

load("va_snowbed_mustela_rodent_gof_nestedloop_temp_pri1.rda")
##############################################################################

# traceplots
par(mfrow=c(2,4))
traceplot(out.gof, parameters = c("GamA","GamB","GamAB","GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
traceplot(out.gof, parameters = c("gamA","gamB","gamAB","gamBA", "epsA", "epsB", "epsAB", "epsBA"))

# create df in a format suitable for ggplot
n <- 5700

dat <- data.frame(sims=unlist(out.gof$sims.list[c(1:8,10:17)]),
                  par=c(rep("gamA", times=n),rep("gamB", times=n),rep("gamAB", times=n),rep("gamBA", times=n),
                        rep("epsA", times=n),rep("epsB", times=n),rep("epsAB", times=n),rep("epsBA", times=n),
                        rep("GamA", times=n),rep("GamB", times=n),rep("GamAB", times=n),rep("GamBA", times=n),
                        rep("EpsA", times=n),rep("EpsB", times=n),rep("EpsAB", times=n),rep("EpsBA", times=n)),
                  level=c(rep("Site parameters", times=n*8), rep("Block parameters", times=n*8)))


dat2 <- data.frame(sims=unlist(out.gof$mean[c(1:8,10:17)]),
                   par=c("gamA", "gamB", "gamAB", "gamBA", "epsA", "epsB", "epsAB", "epsBA",
                         "GamA", "GamB", "GamAB", "GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"),
                   level=c(rep("Site parameters", times=8), rep("Block parameters", times=8)))

# plot colonization and extinction probabilities
ggplot(data=dat, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  facet_wrap(~level, scales="free_x")+
   theme(strip.text = element_text(size=30,face="bold"),
         axis.text.x = element_text(size=40, face="bold"),
         axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete( labels=c('gamA' = parse(text = TeX("$\\gamma_{A}$")), 'gamB' = parse(text = TeX("$\\gamma_{B}$")), 
                                        'gamAB' = parse(text = TeX("$\\gamma_{A|B}$")), 'gamBA' = parse(text = TeX("$\\gamma_{B|A}$")),
                                        'epsA' = parse(text = TeX("$\\epsilon_{A}$")), 'epsB' = parse(text = TeX("$\\epsilon_{B}$")),
                                        'epsAB' = parse(text = TeX("$\\epsilon_{A|B}$")), 'epsBA' = parse(text = TeX("$\\epsilon_{B|A}$")),
                                        'GamA' = parse(text = TeX("$\\Gamma_{A}$")), 'GamB' = parse(text = TeX("$\\Gamma_{B}$")),
                                        'GamAB' = parse(text = TeX("$\\Gamma_{A|B}$")), 'GamBA' = parse(text = TeX("$\\Gamma_{B|A}$")),
                                        'EpsA' = parse(text = TeX("$E_{A}$")), 'EpsB' = parse(text = TeX("$E_{B}$")), 
                                        'EpsAB' = parse(text = TeX("$E_{A|B}$")), 'EpsBA' = parse(text = TeX("$E_{B|A}$")) ))

# save plot
setwd("../plot")
ggsave("modperf_va_snowbed_mustela_rodent_site&block.png", width = 60, height = 30, units="cm")

###########################################################################################################################
##    Plot estimated detection probabilities

# making df in a suitable format for ggplot
dat3 <- data.frame(sims=c(unlist(out.gof$sims.list[18][[1]][,1,1,c(1,15)]),
                          unlist(out.gof$sims.list[19][[1]][,1,1,c(1,15)])),
                   par=c(rep("pA_S", times=n), rep("pA_W", times=n), rep("pB_S", times=n), rep("pB_W", times=n)) )

dat4 <- data.frame(sims=c(unlist(out.gof$mean[18][[1]][1,1,c(1,15)]),
                          unlist(out.gof$mean[19][[1]][1,1,c(1,15)])),
                   par=c("pA_S", "pA_W", "pB_S", "pB_W"))

# specifing what order the parameters should be plotted
pos3 <-c("pA_S", "pA_W", "pB_S", "pB_W")  

#make plot
ggplot(data=dat3, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey", size=1.5)+
  geom_point(data=dat4, color="red", shape="-", size=40)+
  geom_boxplot(data=dat3, width=0.4, color="black", alpha=0, fill="white")+
  labs(y="", x="")+
  theme(axis.text.x = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete(limits=pos3, labels=c('pA_S' = parse(text = TeX("$p_{A_S}$")), 'pA_W' = parse(text = TeX("$p_{A_W}$")), 
                                         'pB_S' = parse(text = TeX("$p_{B_S}$")), 'pB_W' = parse(text = TeX("$p_{B_W}$"))))

# save plot
ggsave("modperf_va_snowbed_mustela_rodent_p.png", width = 60, height = 20, units="cm")


## plot estimated initial occupancy probability

# make df with a suitable structure for ggplot
dat5 <- data.frame(sims=c(as.vector(out.gof$sims.list$psi[,,1]),
                          as.vector(out.gof$sims.list$psi[,,2]),
                          as.vector(out.gof$sims.list$psi[,,3])),
                   par=c(rep("psi1", times=n*8), rep("psi2", times=n*8), rep("psi3", times=n*8)),
                   block=c(rep(c(rep("b1", times=n),rep("b2",times=n),rep("b3", times=n),rep("b4", times=n),rep("b5", times=n),rep("b6", times=n),rep("b7", times=n),rep("b8", times=n))
                               , times=3)))

# spesify what order to plot the simulations
pos5 <-c("psi1","psi2","psi3")  

# make plot
ggplot(data=dat5, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_boxplot(data=dat5, width=0.4, color="black", alpha=0.4, fill="white")+
  #geom_point(data=dat6, color="red", shape="-", size=15)+
  labs(y="", x="")+
  facet_wrap(~block, labeller = label_parsed, ncol=4)+
  theme(strip.text = element_text(size=30,face="bold"), axis.text.x = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete(limits=pos5, 
                   labels=c('psi1' = parse(text = TeX("$\\psi_{A}$")), 'psi2' = parse(text = TeX("$\\psi_{B}$")), 'psi3' = parse(text = TeX("$\\psi_{A|B}$"))))

# save the model
ggsave("modperf_va_snowbed_mustela_rodent_psi.png", width = 60, height = 30, units="cm")

####
str(out.gof$mean[c(1,3,2,4,5,7,6,8)])
str(out.gof$sd[c(1,3,2,4,5,7,6,8)])
str(out.gof$q2.5[c(1,3,2,4,5,7,6,8)])
str(out.gof$q97.5[c(1,3,2,4,5,7,6,8)])

str(out.gof$mean[c(10,12,11,13,14,16,15,17)])
str(out.gof$sd[c(10,12,11,13,14,16,15,17)])
str(out.gof$q2.5[c(10,12,11,13,14,16,15,17)])
str(out.gof$q97.5[c(10,12,11,13,14,16,15,17)])
##

# the interaction

# 95% CRI's

quantile(out.gof$sims.list$diff_gamA, probs=c(0.5,0.05,0.95)) 
quantile(out.gof$sims.list$diff_gamB, probs=c(0.5,0.05,0.95)) 
quantile(out.gof$sims.list$diff_epsA, probs=c(0.5,0.05,0.95)) 
quantile(out.gof$sims.list$diff_epsB, probs=c(0.5,0.05,0.95))

quantile(out.gof$sims.list$diff_GamA, probs=c(0.5,0.05,0.95)) 
quantile(out.gof$sims.list$diff_GamB, probs=c(0.5,0.05,0.95)) 
quantile(out.gof$sims.list$diff_EpsA, probs=c(0.5,0.05,0.95)) 
quantile(out.gof$sims.list$diff_EpsB, probs=c(0.5,0.05,0.95))

#~ End of Script