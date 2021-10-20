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
#setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models")
setwd("~/UiT/GitProjects/A-dynamic-and-hierarchical-spatial-occupancy-model-for-interacting-species/case_study/model_output")

#load model
#setwd("./hidden_block_ko_fromautoclass/ko_vj/model_output")
dir()

load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop_ni10k_pri1.rda")
pri1 <- out.gof

load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop_ni10k_pri2.rda")
pri2 <- out.gof

load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop_ni10k_pri3_2.rda")
pri3 <- out.gof

##############################################################################

# traceplots
par(mfrow=c(2,4))
traceplot(pri1, parameters = c("GamA","GamB","GamAB","GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
traceplot(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ, parameters = c("gamA","gamB","gamAB","gamBA", "epsA", "epsB", "epsAB", "epsBA"))

# create df in a format suitable for ggplot
dat <- data.frame(sims=c(as.numeric(unlist(pri1$sims.list[c(1:8,10:17)])),as.numeric(unlist(pri2$sims.list[c(1:8,10:17)])), as.numeric(unlist(pri3$sims.list[c(1:8,10:17)]))),
                  par=c(rep("gamA1", times=1600),rep("gamB1", times=1600), rep("gamAB1", times=1600),rep("gamBA1", times=1600),
                        rep("epsA1", times=1600),rep("epsB1", times=1600), rep("epsAB1", times=1600),rep("epsBA1", times=1600),
                        rep("GamA1", times=1600),rep("GamB1", times=1600), rep("GamAB1", times=1600),rep("GamBA1", times=1600),
                        rep("EpsA1", times=1600),rep("EpsB1", times=1600), rep("EpsAB1", times=1600),rep("EpsBA1", times=1600),
                        rep("gamA2", times=1600),rep("gamB2", times=1600),rep("gamAB2", times=1600),rep("gamBA2", times=1600),
                        rep("epsA2", times=1600),rep("epsB2", times=1600),rep("epsAB2", times=1600),rep("epsBA2", times=1600),
                        rep("GamA2", times=1600),rep("GamB2", times=1600),rep("GamAB2", times=1600),rep("GamBA2", times=1600),
                        rep("EpsA2", times=1600),rep("EpsB2", times=1600),rep("EpsAB2", times=1600),rep("EpsBA2", times=1600),
                        rep("gamA3", times=1600),rep("gamB3", times=1600),rep("gamAB3", times=1600),rep("gamBA3", times=1600),
                        rep("epsA3", times=1600),rep("epsB3", times=1600),rep("epsAB3", times=1600),rep("epsBA3", times=1600),
                        rep("GamA3", times=1600),rep("GamB3", times=1600),rep("GamAB3", times=1600),rep("GamBA3", times=1600),
                        rep("EpsA3", times=1600),rep("EpsB3", times=1600),rep("EpsAB3", times=1600),rep("EpsBA3", times=1600)),
                  par2= rep(c(rep("gamA",times=1600), rep("gamB",time=1600),  rep("gamAB",time=1600), rep("gamBA",time=1600),
                         rep("epsA",times=1600), rep("epsB",time=1600), rep("epsAB",time=1600), rep("epsBA",time=1600),
                         rep("GamA",times=1600), rep("GamB",time=1600), rep("GamAB",time=1600), rep("GamBA",time=1600),
                         rep("EpsA",times=1600), rep("EpsB",times=1600), rep("EpsAB",time=1600), rep("EpsBA",time=1600)), times=3),
                  level=c(rep("Site parameters", times=1600*8), rep("Block parameters", times=1600*8),
                          rep("Site parameters", times=1600*8), rep("Block parameters", times=1600*8),
                          rep("Site parameters", times=1600*8), rep("Block parameters", times=1600*8)))


dat2 <- data.frame(sims=c(as.numeric(unlist(pri1$mean[c(1:8,10:17)])),as.numeric(unlist(pri2$mean[c(1:8,10:17)])),as.numeric(unlist(pri3$mean[c(1:8,10:17)]))),
                   par=c("gamA1", "gamB1",  "gamAB1", "gamBA1","epsA1", "epsB1", "epsAB1", "epsBA1",
                         "GamA1", "GamB1", "GamAB1", "GamBA1", "EpsA1", "EpsB1", "EpsAB1", "EpsBA1",
                         "gamA2", "gamB2",  "gamAB2", "gamBA2","epsA2", "epsB2", "epsAB2", "epsBA2",
                         "GamA2", "GamB2", "GamAB2", "GamBA2", "EpsA2", "EpsB2", "EpsAB2", "EpsBA2",
                         "gamA3", "gamB3",  "gamAB3", "gamBA3","epsA3", "epsB3", "epsAB3", "epsBA3",
                         "GamA3", "GamB3", "GamAB3", "GamBA3", "EpsA3", "EpsB3", "EpsAB3", "EpsBA3"),
                   par2=c(rep(c("gamA", "gamB",  "gamAB", "gamBA","epsA", "epsB", "epsAB", "epsBA",
                          "GamA", "GamB", "GamAB", "GamBA","EpsA", "EpsB", "EpsAB", "EpsBA"), times=3)),
                   level=c(rep("Site parameters", times=8), rep("Block parameters", times=8),
                           rep("Site parameters", times=8), rep("Block parameters", times=8),
                           rep("Site parameters", times=8), rep("Block parameters", times=8)))

# plot colonization and extinction probabilities
ggplot(data=dat, aes(x=par, y=sims, fill=par2))+
  geom_violin(aes(colour=par2))+
  geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  facet_wrap(~level, scales="free_x")+
   theme(strip.text = element_text(size=30,face="bold"),
         axis.text.x = element_text(size=25, face="bold"),
         axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete( labels=c('gamA1' = parse(text = TeX("$\\gamma_{A1}$")), 'gamA2' = parse(text = TeX("$\\gamma_{A2}$")),'gamA3' = parse(text = TeX("$\\gamma_{A3}$")),
                             'gamB1' = parse(text = TeX("$\\gamma_{B1}$")),'gamB2' = parse(text = TeX("$\\gamma_{B2}$")), 'gamB3' = parse(text = TeX("$\\gamma_{B3}$")),
                              'gamAB1' = parse(text = TeX("$\\gamma_{AB1}$")),'gamAB2' = parse(text = TeX("$\\gamma_{AB2}$")),'gamAB3' = parse(text = TeX("$\\gamma_{AB3}$")),
                             'gamBA1' = parse(text = TeX("$\\gamma_{BA1}$")),'gamBA2' = parse(text = TeX("$\\gamma_{BA2}$")),'gamBA3' = parse(text = TeX("$\\gamma_{BA3}$")),
                              'epsA1' = parse(text = TeX("$\\epsilon_{A1}$")),'epsA2' = parse(text = TeX("$\\epsilon_{A2}$")),'epsA3' = parse(text = TeX("$\\epsilon_{A3}$")),
                             'epsB1' = parse(text = TeX("$\\epsilon_{B1}$")),'epsB2' = parse(text = TeX("$\\epsilon_{B2}$")),'epsB3' = parse(text = TeX("$\\epsilon_{B3}$")),
                              'epsAB1' = parse(text = TeX("$\\epsilon_{AB1}$")),'epsAB2' = parse(text = TeX("$\\epsilon_{AB2}$")),'epsAB3' = parse(text = TeX("$\\epsilon_{AB3}$")),
                             'epsBA1' = parse(text = TeX("$\\epsilon_{BA1}$")),'epsBA2' = parse(text = TeX("$\\epsilon_{BA2}$")),'epsBA3' = parse(text = TeX("$\\epsilon_{BA3}$")),
                             'GamA1' = parse(text = TeX("$\\Gamma_{A1}$")), 'GamA2' = parse(text = TeX("$\\Gamma_{A2}$")),'GamA3' = parse(text = TeX("$\\Gamma_{A3}$")),
                             'GamB1' = parse(text = TeX("$\\Gamma_{B1}$")),'GamB2' = parse(text = TeX("$\\Gamma_{B2}$")),'GamB3' = parse(text = TeX("$\\Gamma_{B3}$")),
                              'GamAB1' = parse(text = TeX("$\\Gamma_{AB1}$")),'GamAB2' = parse(text = TeX("$\\Gamma_{AB2}$")),'GamAB3' = parse(text = TeX("$\\Gamma_{AB3}$")),
                             'GamBA1' = parse(text = TeX("$\\Gamma_{BA1}$")),'GamBA2' = parse(text = TeX("$\\Gamma_{BA2}$")),'GamBA3' = parse(text = TeX("$\\Gamma_{BA3}$")),
                              'EpsA1' = parse(text = TeX("$E_{A1}$")),'EpsA2' = parse(text = TeX("$E_{A2}$")),'EpsA3' = parse(text = TeX("$E_{A3}$")),
                             'EpsB1' = parse(text = TeX("$E_{B1}$")),'EpsB2' = parse(text = TeX("$E_{B2}$")), 'EpsB3' = parse(text = TeX("$E_{B3}$")),
                              'EpsAB1' = parse(text = TeX("$E_{AB1}$")),'EpsAB2' = parse(text = TeX("$E_{AB2}$")),'EpsAB3' = parse(text = TeX("$E_{AB3}$")),
                             'EpsBA1' = parse(text = TeX("$E_{BA1}$")),'EpsBA2' = parse(text = TeX("$E_{BA2}$")),'EpsBA3' = parse(text = TeX("$E_{BA3}$")) ))

# save plot
setwd("../plot")
ggsave("modperf_site&block_prisens.png", width = 90, height = 30, units="cm")

###########################################################################################################################
##    Plot estimated detection probabilities

# making df in a suitable format for ggplot
dat3 <- data.frame(sims=c(unlist(pri1$sims.list[18][[1]][,c(1,15)]), unlist(pri1$sims.list[19][[1]][,c(1,15)]),
                          unlist(pri2$sims.list[18][[1]][,c(1,15)]), unlist(pri2$sims.list[19][[1]][,c(1,15)]),
                          unlist(pri3$sims.list[18][[1]][,c(1,15)]), unlist(pri3$sims.list[19][[1]][,c(1,15)])),
                   par=c(rep("pA_S1", times=1600), rep("pA_W1", times=1600), rep("pB_S1", times=1600), rep("pB_W1", times=1600),
                         rep("pA_S2", times=1600), rep("pA_W2", times=1600), rep("pB_S2", times=1600), rep("pB_W2", times=1600),
                         rep("pA_S3", times=1600), rep("pA_W3", times=1600), rep("pB_S3", times=1600), rep("pB_W3", times=1600)),
                   par2=rep(c(rep("pA_S", times=1600), rep("pA_W", times=1600), rep("pB_S", times=1600), rep("pB_W", times=1600)), times=3) )

dat4 <- data.frame(sims=c(unlist(pri1$mean[18][[1]][c(1,15)]), unlist(pri1$mean[19][[1]][c(1,15)]),
                          unlist(pri2$mean[18][[1]][c(1,15)]), unlist(pri2$mean[19][[1]][c(1,15)]),
                          unlist(pri3$mean[18][[1]][c(1,15)]), unlist(pri3$mean[19][[1]][c(1,15)])),
                   par=c("pA_S1", "pA_W1", "pB_S1", "pB_W1","pA_S2", "pA_W2", "pB_S2", "pB_W2","pA_S3", "pA_W3", "pB_S3", "pB_W3"),
                   par2=rep(c("pA_S", "pA_W", "pB_S", "pB_W"), times=3) )

# specifing what order the parameters should be plotted
#pos3 <-c("pA_S", "pA_W", "pB_S", "pB_W")  

#make plot
ggplot(data=dat3, aes(x=par, y=sims, fill=par2))+
  geom_violin(aes(colour=par2), size=1.5)+
  geom_point(data=dat4, color="red", shape="-", size=20)+
  labs(y="", x="")+
  theme(axis.text.x = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete( labels=c('pA_S1' = parse(text = TeX("$p_{A_S1}$")),'pA_S2' = parse(text = TeX("$p_{A_S2}$")),'pA_S3' = parse(text = TeX("$p_{A_S3}$")),
                                         'pA_W1' = parse(text = TeX("$p_{A_W1}$")), 'pA_W2' = parse(text = TeX("$p_{A_W2}$")),'pA_W3' = parse(text = TeX("$p_{A_W3}$")),
                                         'pB_S1' = parse(text = TeX("$p_{B_S1}$")),'pB_S2' = parse(text = TeX("$p_{B_S2}$")),'pB_S3' = parse(text = TeX("$p_{B_S3}$")),
                                         'pB_W1' = parse(text = TeX("$p_{B_W1}$")), 'pB_W2' = parse(text = TeX("$p_{B_W2}$")), 'pB_W3' = parse(text = TeX("$p_{B_W3}$"))))

# save plot
ggsave("modperf_p_priorsens.png", width = 60, height = 20, units="cm")


## plot estimated initial occupancy probability

# make df with a suitable structure for ggplot
dat5 <- data.frame(sims=c(as.vector(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$psi[,,1]),
                          as.vector(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$psi[,,2]),
                          as.vector(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$psi[,,3])),
                   par=c(rep("psi1", times=120000), rep("psi2", times=120000), rep("psi3", times=120000)),
                   block=c(rep(c(rep("b1", times=15000),rep("b2",times=15000),rep("b3", times=15000),rep("b4", times=15000),rep("b5", times=15000),rep("b6", times=15000),rep("b7", times=15000),rep("b8", times=15000))
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
                   labels=c('psi1' = parse(text = TeX("$\\psi_{A}$")), 'psi2' = parse(text = TeX("$\\psi_{B}$")), 'psi3' = parse(text = TeX("$\\psi_{AB}$"))))

# save the model
ggsave("modperf_va_snowbed_mustela_rodent_sdet_4stpm_s1_203_psi_rmBQ.png", width = 60, height = 30, units="cm")

####
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$mean[c(1,3,2,4,5,7,6,8)])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sd[c(1,3,2,4,5,7,6,8)])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$q2.5[c(1,3,2,4,5,7,6,8)])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$q97.5[c(1,3,2,4,5,7,6,8)])

str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$mean[c(10,12,11,13,14,16,15,17)])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sd[c(10,12,11,13,14,16,15,17)])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$q2.5[c(10,12,11,13,14,16,15,17)])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$q97.5[c(10,12,11,13,14,16,15,17)])
##

# the interaction

# 90% CI's

quantile(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$diff_gamA, probs=c(0.5,0.1,0.9)) 
quantile(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$diff_gamB, probs=c(0.5,0.1,0.9)) 
quantile(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$diff_epsA, probs=c(0.5,0.1,0.9)) 
quantile(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$diff_epsB, probs=c(0.5,0.1,0.9))

quantile(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$diff_GamA, probs=c(0.5,0.1,0.9)) 
quantile(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$diff_GamB, probs=c(0.5,0.1,0.9)) 
quantile(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$diff_EpsA, probs=c(0.5,0.1,0.9)) 
quantile(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ$sims.list$diff_EpsB, probs=c(0.5,0.1,0.9))

#~ End of Script