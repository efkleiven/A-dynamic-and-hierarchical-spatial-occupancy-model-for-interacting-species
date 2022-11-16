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

load("va_snowbed_mustela_rodent_gof_pri1_0005_2016_2021.rda")
#pri1 <- out.gof

load("va_snowbed_mustela_rodent_gof_pri2_0005_2016_2021.rda")
#pri2 <- out.gof

load("va_snowbed_mustela_rodent_gof_pri3_0005_2016_2021.rda")
#pri3 <- out.gof

##############################################################################

# traceplots
par(mfrow=c(2,4))
traceplot(pri1, parameters = c("GamA","GamB","GamAB","GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
traceplot(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_rmBQ, parameters = c("gamA","gamB","gamAB","gamBA", "epsA", "epsB", "epsAB", "epsBA"))

# create df in a format suitable for ggplot
n=11400
dat <- data.frame(sims=c(as.numeric(unlist(out.gof1$sims.list[c(1:8,10:17)])),as.numeric(unlist(out.gof2$sims.list[c(1:8,10:17)])), as.numeric(unlist(out.gof3$sims.list[c(1:8,10:17)]))),
                  par=c(rep("gamA1", times=n),rep("gamB1", times=n), rep("gamAB1", times=n),rep("gamBA1", times=n),
                        rep("epsA1", times=n),rep("epsB1", times=n), rep("epsAB1", times=n),rep("epsBA1", times=n),
                        rep("GamA1", times=n),rep("GamB1", times=n), rep("GamAB1", times=n),rep("GamBA1", times=n),
                        rep("EpsA1", times=n),rep("EpsB1", times=n), rep("EpsAB1", times=n),rep("EpsBA1", times=n),
                        rep("gamA2", times=n),rep("gamB2", times=n),rep("gamAB2", times=n),rep("gamBA2", times=n),
                        rep("epsA2", times=n),rep("epsB2", times=n),rep("epsAB2", times=n),rep("epsBA2", times=n),
                        rep("GamA2", times=n),rep("GamB2", times=n),rep("GamAB2", times=n),rep("GamBA2", times=n),
                        rep("EpsA2", times=n),rep("EpsB2", times=n),rep("EpsAB2", times=n),rep("EpsBA2", times=n),
                        rep("gamA3", times=n),rep("gamB3", times=n),rep("gamAB3", times=n),rep("gamBA3", times=n),
                        rep("epsA3", times=n),rep("epsB3", times=n),rep("epsAB3", times=n),rep("epsBA3", times=n),
                        rep("GamA3", times=n),rep("GamB3", times=n),rep("GamAB3", times=n),rep("GamBA3", times=n),
                        rep("EpsA3", times=n),rep("EpsB3", times=n),rep("EpsAB3", times=n),rep("EpsBA3", times=n)),
                  par2= rep(c(rep("gamA",times=n), rep("gamB",time=n),  rep("gamAB",time=n), rep("gamBA",time=n),
                         rep("epsA",times=n), rep("epsB",time=n), rep("epsAB",time=n), rep("epsBA",time=n),
                         rep("GamA",times=n), rep("GamB",time=n), rep("GamAB",time=n), rep("GamBA",time=n),
                         rep("EpsA",times=n), rep("EpsB",times=n), rep("EpsAB",time=n), rep("EpsBA",time=n)), times=3),
                  level=c(rep("Site parameters", times=n*8), rep("Block parameters", times=n*8),
                          rep("Site parameters", times=n*8), rep("Block parameters", times=n*8),
                          rep("Site parameters", times=n*8), rep("Block parameters", times=n*8)))


dat2 <- data.frame(sims=c(as.numeric(unlist(out.gof1$mean[c(1:8,10:17)])),as.numeric(unlist(out.gof2$mean[c(1:8,10:17)])),as.numeric(unlist(out.gof3$mean[c(1:8,10:17)]))),
                   par=c("gamA1", "gamB1",  "gamAB1", "gamBA1","epsA1", "epsB1", "epsAB1", "epsBA1",
                         "GamA1", "GamB1", "GamAB1", "GamBA1", "EpsA1", "EpsB1", "EpsAB1", "EpsBA1",
                         "gamA2", "gamB2",  "gamAB2", "gamBA2","epsA2", "epsB2", "epsAB2", "epsBA2",
                         "GamA2", "GamB2", "GamAB2", "GamBA2", "EpsA2", "EpsB2", "EpsAB2", "EpsBA2",
                         "gamA3", "gamB3",  "gamAB3", "gamBA3","epsA3", "epsB3", "epsAB3", "epsBA3",
                         "GamA3", "GamB3", "GamAB3", "GamBA3", "EpsA3", "EpsB3", "EpsAB3", "EpsBA3"),
                   par2=c(rep(c("gamA", "gamB",  "gamAB", "gamBA","EpsA", "EpsB", "EpsAB", "EpsBA",
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
         axis.text.x = element_text(size=28, face="bold"),
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
setwd("../../plot/revision2")
ggsave("modperf_site&block_prisens.png", width = 90, height = 30, units="cm")


# prior posterior check
# first for site parameters
dat11 <- data.frame(sims=c(as.numeric(unlist(out.gof1$sims.list[c(1:8)])),as.numeric(unlist(out.gof2$sims.list[c(1:8)])), as.numeric(unlist(out.gof3$sims.list[c(1:8)]))),
                  par=c(rep("gamA1", times=n),rep("gamB1", times=n), rep("gamAB1", times=n),rep("gamBA1", times=n),
                        rep("epsA1", times=n),rep("epsB1", times=n), rep("epsAB1", times=n),rep("epsBA1", times=n),
                        rep("gamA2", times=n),rep("gamB2", times=n),rep("gamAB2", times=n),rep("gamBA2", times=n),
                        rep("epsA2", times=n),rep("epsB2", times=n),rep("epsAB2", times=n),rep("epsBA2", times=n),
                        rep("gamA3", times=n),rep("gamB3", times=n),rep("gamAB3", times=n),rep("gamBA3", times=n),
                        rep("epsA3", times=n),rep("epsB3", times=n),rep("epsAB3", times=n),rep("epsBA3", times=n)),
                  par2= rep(c(rep("gamA",times=n), rep("gamB",time=n),  rep("gamAB",time=n), rep("gamBA",time=n),
                              rep("epsA",times=n), rep("epsB",time=n), rep("epsAB",time=n), rep("epsBA",time=n)), times=3) )


dat12 <- data.frame(sims=c(rep(runif(n, 0,1), times=8), rep(rbeta(n, 4,4), times=8), rbeta(n, 4,2), rbeta(n, 2,4), rbeta(n, 2,4), rbeta(n, 4,2), rbeta(n, 2,4), rbeta(n, 4,2), rbeta(n, 4,2), rbeta(n, 2,4)),
                    par=c(rep("gamA1T", times=n),rep("gamB1T", times=n), rep("gamAB1T", times=n),rep("gamBA1T", times=n),
                          rep("epsA1T", times=n),rep("epsB1T", times=n), rep("epsAB1T", times=n),rep("epsBA1T", times=n),
                          rep("gamA2T", times=n),rep("gamB2T", times=n),rep("gamAB2T", times=n),rep("gamBA2T", times=n),
                          rep("epsA2T", times=n),rep("epsB2T", times=n),rep("epsAB2T", times=n),rep("epsBA2T", times=n),
                          rep("gamA3T", times=n),rep("gamB3T", times=n),rep("gamAB3T", times=n),rep("gamBA3T", times=n),
                          rep("epsA3T", times=n),rep("epsB3T", times=n),rep("epsAB3T", times=n),rep("epsBA3T", times=n)),
                    par2= rep(c(rep("gamA",times=n), rep("gamB",time=n),  rep("gamAB",time=n), rep("gamBA",time=n),
                                rep("epsA",times=n), rep("epsB",time=n), rep("epsAB",time=n), rep("epsBA",time=n)), times=3) )

dat13 <- rbind(dat11,dat12)

par_names <- list('gamA' = parse(text = TeX("$\\gamma_{A}$")), 'gamB' = parse(text = TeX("$\\gamma_{B}$")), 
                  'gamAB' = parse(text = TeX("$\\gamma_{A|B}$")), 'gamBA' = parse(text = TeX("$\\gamma_{B|A}$")),
                  'epsA' = parse(text = TeX("$\\epsilon_{A}$")), 'epsB' = parse(text = TeX("$\\epsilon_{B}$")),
                  'epsAB' = parse(text = TeX("$\\epsilon_{A|B}$")), 'epsBA' = parse(text = TeX("$\\epsilon_{B|A}$")),
                  'GamA' = parse(text = TeX("$\\Gamma_{A}$")), 'GamB' = parse(text = TeX("$\\Gamma_{B}$")),
                  'GamAB' = parse(text = TeX("$\\Gamma_{A|B}$")), 'GamBA' = parse(text = TeX("$\\Gamma_{B|A}$")),
                  'EpsA' = parse(text = TeX("$E_{A}$")), 'EpsB' = parse(text = TeX("$E_{B}$")), 
                  'EpsAB' = parse(text = TeX("$E_{A|B}$")), 'EpsBA' = parse(text = TeX("$E_{B|A}$")) )

par_labeller <- function(variable,value){
  return(par_names[value])
}


ggplot(data=dat13, aes(x=par, y=sims, fill=par2))+
  geom_violin(aes(colour=par2))+
  #geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  facet_wrap(~par2, labeller = par_labeller, ncol=4, scales="free_x")+
  theme(strip.text = element_text(size=30,face="bold"),
        axis.text.x = element_text(size=28, face="bold"),
        axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete( labels=c('gamA1' = parse(text = TeX("$\\gamma_{A_1}$")),       'gamA2' = parse(text = TeX("$\\gamma_{A_2}$")),           'gamA3' = parse(text = TeX("$\\gamma_{A_3}$")),
                             'gamA1T' = parse(text = TeX("$\\gamma_{A_{1T}}$")),   'gamA2T' = parse(text = TeX("$\\gamma_{A_{2T}}$")),       'gamA3T' = parse(text = TeX("$\\gamma_{A_{3T}}$")),
                             'gamB1' = parse(text = TeX("$\\gamma_{B_1}$")),       'gamB2' = parse(text = TeX("$\\gamma_{B_2}$")),           'gamB3' = parse(text = TeX("$\\gamma_{B_3}$")),
                             'gamB1T' = parse(text = TeX("$\\gamma_{B_{1T}}$")),   'gamB2T' = parse(text = TeX("$\\gamma_{B_{2T}$")),        'gamB3T' = parse(text = TeX("$\\gamma_{B_{3T}}$")),
                             'gamAB1' = parse(text = TeX("$\\gamma_{A|B_1}$")),    'gamAB2' = parse(text = TeX("$\\gamma_{A|B_2}$")),        'gamAB3' = parse(text = TeX("$\\gamma_{A|B_3}$")),
                             'gamAB1T' = parse(text = TeX("$\\gamma_{A|B_{1T}}$")),'gamAB2T' = parse(text = TeX("$\\gamma_{A|B_{2T}}$")),    'gamAB3T' = parse(text = TeX("$\\gamma_{A|B_{3T}}$")),
                             'gamBA1' = parse(text = TeX("$\\gamma_{B|A_1}$")),    'gamBA2' = parse(text = TeX("$\\gamma_{B|A_2}$")),        'gamBA3' = parse(text = TeX("$\\gamma_{B|A_3}$")),
                             'gamBA1T' = parse(text = TeX("$\\gamma_{B|A_{1T}}$")),'gamBA2T' = parse(text = TeX("$\\gamma_{B|A_{2T}}$")),    'gamBA3T' = parse(text = TeX("$\\gamma_{B1A_{3T}}$")),
                             'epsA1' = parse(text = TeX("$\\epsilon_{A_1}$")),     'epsA2' = parse(text = TeX("$\\epsilon_{A_2}$")),         'epsA3' = parse(text = TeX("$\\epsilon_{A_3}$")),
                             'epsA1T' = parse(text = TeX("$\\epsilon_{A_{1T}}$")), 'epsA2T' = parse(text = TeX("$\\epsilon_{A_{2T}}$")),     'epsA3T' = parse(text = TeX("$\\epsilon_{A_{3T}}$")),
                             'epsB1' = parse(text = TeX("$\\epsilon_{B_1}$")),     'epsB2' = parse(text = TeX("$\\epsilon_{B_2}$")),         'epsB3' = parse(text = TeX("$\\epsilon_{B_3}$")),
                             'epsB1T' = parse(text = TeX("$\\epsilon_{B_{1T}}$")), 'epsB2T' = parse(text = TeX("$\\epsilon_{B_{2T}}$")),     'epsB3T' = parse(text = TeX("$\\epsilon_{B_{3T}}$")),
                             'epsAB1' = parse(text = TeX("$\\epsilon_{A|B_1}$")),  'epsAB2' = parse(text = TeX("$\\epsilon_{A|B_2}$")),      'epsAB3' = parse(text = TeX("$\\epsilon_{A|B_3}$")),
                             'epsAB1T' = parse(text = TeX("$\\epsilon_{A|B_{1T}}$")),'epsAB2T' = parse(text = TeX("$\\epsilon_{A|B_{2T}}$")),'epsAB3T' = parse(text = TeX("$\\epsilon_{A|B_{3T}}$")),
                             'epsBA1' = parse(text = TeX("$\\epsilon_{B|A_{1}}$")), 'epsBA2' = parse(text = TeX("$\\epsilon_{B|A_2}$")),     'epsBA3' = parse(text = TeX("$\\epsilon_{B|A_3}$")),
                             'epsBA1T' = parse(text = TeX("$\\epsilon_{B|A_{1T}}$")),'epsBA2T' = parse(text = TeX("$\\epsilon_{B|A_{2T}}$")),'epsBA3T' = parse(text = TeX("$\\epsilon_{B|A_{3T}}$")) ))

ggsave("modperf_site_pripost.png", width = 60, height = 30, units="cm")


# then for block parameters

dat21 <- data.frame(sims=c(as.numeric(unlist(out.gof1$sims.list[c(10:17)])),as.numeric(unlist(out.gof2$sims.list[c(10:17)])), as.numeric(unlist(out.gof3$sims.list[c(10:17)]))),
                    par=c(rep("GamA1", times=n),rep("GamB1", times=n), rep("GamAB1", times=n),rep("GamBA1", times=n),
                          rep("EpsA1", times=n),rep("EpsB1", times=n), rep("EpsAB1", times=n),rep("EpsBA1", times=n),
                          rep("GamA2", times=n),rep("GamB2", times=n),rep("GamAB2", times=n),rep("GamBA2", times=n),
                          rep("EpsA2", times=n),rep("EpsB2", times=n),rep("EpsAB2", times=n),rep("EpsBA2", times=n),
                          rep("GamA3", times=n),rep("GamB3", times=n),rep("GamAB3", times=n),rep("GamBA3", times=n),
                          rep("EpsA3", times=n),rep("EpsB3", times=n),rep("EpsAB3", times=n),rep("EpsBA3", times=n)),
                    par2= rep(c(rep("GamA",times=n), rep("GamB",time=n),  rep("GamAB",time=n), rep("GamBA",time=n),
                                rep("EpsA",times=n), rep("EpsB",time=n), rep("EpsAB",time=n), rep("EpsBA",time=n)), times=3) )

dat22 <- data.frame(sims=c(rep(runif(n, 0,1), times=8), rep(rbeta(n, 4,4), times=8), rep(rbeta(n, 2,4), times=5), rep(rbeta(n, 4,2), times=2), rbeta(n, 2,4)),
                    par=c(rep("GamA1T", times=n),rep("GamB1T", times=n), rep("GamAB1T", times=n),rep("GamBA1T", times=n),
                          rep("EpsA1T", times=n),rep("EpsB1T", times=n), rep("EpsAB1T", times=n),rep("EpsBA1T", times=n),
                          rep("GamA2T", times=n),rep("GamB2T", times=n),rep("GamAB2T", times=n),rep("GamBA2T", times=n),
                          rep("EpsA2T", times=n),rep("EpsB2T", times=n),rep("EpsAB2T", times=n),rep("EpsBA2T", times=n),
                          rep("GamA3T", times=n),rep("GamB3T", times=n),rep("GamAB3T", times=n),rep("GamBA3T", times=n),
                          rep("EpsA3T", times=n),rep("EpsB3T", times=n),rep("EpsAB3T", times=n),rep("EpsBA3T", times=n)),
                    par2= rep(c(rep("GamA",times=n), rep("GamB",time=n),  rep("GamAB",time=n), rep("GamBA",time=n),
                                rep("EpsA",times=n), rep("EpsB",time=n), rep("EpsAB",time=n), rep("EpsBA",time=n)), times=3) )

dat23 <- rbind(dat21,dat22)


par_names <- list('gamA' = parse(text = TeX("$\\gamma_{A}$")), 'gamB' = parse(text = TeX("$\\gamma_{B}$")), 
                  'gamAB' = parse(text = TeX("$\\gamma_{A|B}$")), 'gamBA' = parse(text = TeX("$\\gamma_{B|A}$")),
                  'epsA' = parse(text = TeX("$\\epsilon_{A}$")), 'epsB' = parse(text = TeX("$\\epsilon_{B}$")),
                  'epsAB' = parse(text = TeX("$\\epsilon_{A|B}$")), 'epsBA' = parse(text = TeX("$\\epsilon_{B|A}$")),
                  'GamA' = parse(text = TeX("$\\Gamma_{A}$")), 'GamB' = parse(text = TeX("$\\Gamma_{B}$")),
                  'GamAB' = parse(text = TeX("$\\Gamma_{A|B}$")), 'GamBA' = parse(text = TeX("$\\Gamma_{B|A}$")),
                  'EpsA' = parse(text = TeX("$E_{A}$")), 'EpsB' = parse(text = TeX("$E_{B}$")), 
                  'EpsAB' = parse(text = TeX("$E_{A|B}$")), 'EpsBA' = parse(text = TeX("$E_{B|A}$")) )

par_labeller <- function(variable,value){
  return(par_names[value])
}


ggplot(data=dat23, aes(x=par, y=sims, fill=par2))+
  geom_violin(aes(colour=par2))+
  #geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  facet_wrap(~par2, labeller = par_labeller, ncol=4, scales="free_x")+
  theme(strip.text = element_text(size=30,face="bold"),
        axis.text.x = element_text(size=28, face="bold"),
        axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete( labels=c('GamA1' = parse(text = TeX("$\\Gamma_{A_1}$")),       'GamA2' = parse(text = TeX("$\\Gamma_{A_2}$")),           'GamA3' = parse(text = TeX("$\\Gamma_{A_3}$")),
                             'GamA1T' = parse(text = TeX("$\\Gamma_{A_{1T}}$")),   'GamA2T' = parse(text = TeX("$\\Gamma_{A_{2T}}$")),       'GamA3T' = parse(text = TeX("$\\Gamma_{A_{3T}}$")),
                             'GamB1' = parse(text = TeX("$\\Gamma_{B_1}$")),       'GamB2' = parse(text = TeX("$\\Gamma_{B_2}$")),           'GamB3' = parse(text = TeX("$\\Gamma_{B_3}$")),
                             'GamB1T' = parse(text = TeX("$\\Gamma_{B_{1T}}$")),   'GamB2T' = parse(text = TeX("$\\Gamma_{B_{2T}$")),        'GamB3T' = parse(text = TeX("$\\Gamma_{B_{3T}}$")),
                             'GamAB1' = parse(text = TeX("$\\Gamma_{A|B_1}$")),    'GamAB2' = parse(text = TeX("$\\Gamma_{A|B_2}$")),        'GamAB3' = parse(text = TeX("$\\Gamma_{A|B_3}$")),
                             'GamAB1T' = parse(text = TeX("$\\Gamma_{A|B_{1T}}$")),'GamAB2T' = parse(text = TeX("$\\Gamma_{A|B_{2T}}$")),    'GamAB3T' = parse(text = TeX("$\\Gamma_{A|B_{3T}}$")),
                             'GamBA1' = parse(text = TeX("$\\Gamma_{B|A_1}$")),    'GamBA2' = parse(text = TeX("$\\Gamma_{B|A_2}$")),        'GamBA3' = parse(text = TeX("$\\Gamma_{B|A_3}$")),
                             'GamBA1T' = parse(text = TeX("$\\Gamma_{B|A_{1T}}$")),'GamBA2T' = parse(text = TeX("$\\Gamma_{B|A_{2T}}$")),    'GamBA3T' = parse(text = TeX("$\\Gamma_{B1A_{3T}}$")),
                             'EpsA1' = parse(text = TeX("$\\Epsilon_{A_1}$")),     'EpsA2' = parse(text = TeX("$\\Epsilon_{A_2}$")),         'EpsA3' = parse(text = TeX("$\\Epsilon_{A_3}$")),
                             'EpsA1T' = parse(text = TeX("$\\Epsilon_{A_{1T}}$")), 'EpsA2T' = parse(text = TeX("$\\Epsilon_{A_{2T}}$")),     'EpsA3T' = parse(text = TeX("$\\Epsilon_{A_{3T}}$")),
                             'EpsB1' = parse(text = TeX("$\\Epsilon_{B_1}$")),     'EpsB2' = parse(text = TeX("$\\Epsilon_{B_2}$")),         'EpsB3' = parse(text = TeX("$\\Epsilon_{B_3}$")),
                             'EpsB1T' = parse(text = TeX("$\\Epsilon_{B_{1T}}$")), 'EpsB2T' = parse(text = TeX("$\\Epsilon_{B_{2T}}$")),     'EpsB3T' = parse(text = TeX("$\\Epsilon_{B_{3T}}$")),
                             'EpsAB1' = parse(text = TeX("$\\Epsilon_{A|B_1}$")),  'EpsAB2' = parse(text = TeX("$\\Epsilon_{A|B_2}$")),      'EpsAB3' = parse(text = TeX("$\\Epsilon_{A|B_3}$")),
                             'EpsAB1T' = parse(text = TeX("$\\Epsilon_{A|B_{1T}}$")),'EpsAB2T' = parse(text = TeX("$\\Epsilon_{A|B_{2T}}$")),'EpsAB3T' = parse(text = TeX("$\\Epsilon_{A|B_{3T}}$")),
                             'EpsBA1' = parse(text = TeX("$\\Epsilon_{B|A_{1}}$")), 'EpsBA2' = parse(text = TeX("$\\Epsilon_{B|A_2}$")),     'EpsBA3' = parse(text = TeX("$\\Epsilon_{B|A_3}$")),
                             'EpsBA1T' = parse(text = TeX("$\\Epsilon_{B|A_{1T}}$")),'EpsBA2T' = parse(text = TeX("$\\Epsilon_{B|A_{2T}}$")),'EpsBA3T' = parse(text = TeX("$\\Epsilon_{B|A_{3T}}$")) ))

ggsave("modperf_block_pripost.png", width = 60, height = 30, units="cm")

###########################################################################################################################
##    Plot estimated detection probabilities

dim(pri1$sims.list[18][[1]])

# making df in a suitable format for ggplot
dat3 <- data.frame(sims=c(unlist(out.gof1$sims.list[18][[1]][,1,1,c(1,15)]), unlist(out.gof1$sims.list[19][[1]][,1,1,c(1,15)]),
                          unlist(out.gof2$sims.list[18][[1]][,1,1,c(1,15)]), unlist(out.gof2$sims.list[19][[1]][,1,1,c(1,15)]),
                          unlist(out.gof3$sims.list[18][[1]][,1,1,c(1,15)]), unlist(out.gof3$sims.list[19][[1]][,1,1,c(1,15)])),
                   par=c(rep("pA_S1", times=n), rep("pA_W1", times=n), rep("pB_S1", times=n), rep("pB_W1", times=n),
                         rep("pA_S2", times=n), rep("pA_W2", times=n), rep("pB_S2", times=n), rep("pB_W2", times=n),
                         rep("pA_S3", times=n), rep("pA_W3", times=n), rep("pB_S3", times=n), rep("pB_W3", times=n)),
                   par2=rep(c(rep("pA_S", times=n), rep("pA_W", times=n), rep("pB_S", times=n), rep("pB_W", times=n)), times=3) )

dat4 <- data.frame(sims=c(unlist(out.gof1$mean[18][[1]][1,1,c(1,15)]), unlist(out.gof1$mean[19][[1]][1,1,c(1,15)]),
                          unlist(out.gof2$mean[18][[1]][1,1,c(1,15)]), unlist(out.gof2$mean[19][[1]][1,1,c(1,15)]),
                          unlist(out.gof3$mean[18][[1]][1,1,c(1,15)]), unlist(out.gof3$mean[19][[1]][1,1,c(1,15)])),
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


# plot prior posterior overlap

dat31 <- data.frame(sims=c(unlist(out.gof1$sims.list[18][[1]][,1,1,c(1,15)]), unlist(out.gof1$sims.list[19][[1]][,1,1,c(1,15)]),
                          unlist(out.gof2$sims.list[18][[1]][,1,1,c(1,15)]), unlist(out.gof2$sims.list[19][[1]][,1,1,c(1,15)]),
                          unlist(out.gof3$sims.list[18][[1]][,1,1,c(1,15)]), unlist(out.gof3$sims.list[19][[1]][,1,1,c(1,15)])),
                   par=c(rep("pA_S1", times=n), rep("pA_W1", times=n), rep("pB_S1", times=n), rep("pB_W1", times=n),
                         rep("pA_S2", times=n), rep("pA_W2", times=n), rep("pB_S2", times=n), rep("pB_W2", times=n),
                         rep("pA_S3", times=n), rep("pA_W3", times=n), rep("pB_S3", times=n), rep("pB_W3", times=n)),
                   par2=rep(c(rep("pA_S", times=n), rep("pA_W", times=n), rep("pB_S", times=n), rep("pB_W", times=n)), times=3) )

dat32 <- data.frame(sims=c(rep(plogis(rnorm(n,0,1)), times=4), rep(plogis(rlogis(n,0,1)), times=4), rep(plogis(rbeta(n,2,4)), times=4) ),
                    par=c(rep("pA_S1T", times=n), rep("pA_W1T", times=n), rep("pB_S1T", times=n), rep("pB_W1T", times=n),
                          rep("pA_S2T", times=n), rep("pA_W2T", times=n), rep("pB_S2T", times=n), rep("pB_W2T", times=n),
                          rep("pA_S3T", times=n), rep("pA_W3T", times=n), rep("pB_S3T", times=n), rep("pB_W3T", times=n)),
                    par2=rep(c(rep("pA_S", times=n), rep("pA_W", times=n), rep("pB_S", times=n), rep("pB_W", times=n)), times=3) )


dat33 <- rbind(dat31,dat32)


par_names <- list('pA_S' = parse(text = TeX("$\\p_{A_S}$")), 'pA_W' = parse(text = TeX("$\\p_{A_W}$")), 
                  'pB_S' = parse(text = TeX("$\\p_{B_S}$")), 'pB_W' = parse(text = TeX("$\\p_{B_W}$")) )

par_labeller <- function(variable,value){
  return(par_names[value])
}


ggplot(data=dat33, aes(x=par, y=sims, fill=par2))+
  geom_violin(aes(colour=par2))+
  labs(y="", x="")+
  facet_wrap(~par2, labeller = par_labeller, ncol=4, scales="free_x")+
  theme(strip.text = element_text(size=30,face="bold"),
        axis.text.x = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete( labels=c('pA_S1' = parse(text = TeX("$\\p_{A_S_{1}}$")), 'pA_S2' = parse(text = TeX("$\\p_{A_S_{2}}$")), 'pA_S3' = parse(text = TeX("$\\p_{A_S{3}}$")),
                             'pA_W1' = parse(text = TeX("$\\p_{A_W_{1}}$")), 'pA_W2' = parse(text = TeX("$\\p_{A_W_{2}}$")), 'pA_W3' = parse(text = TeX("$\\p_{A_W{3}}$")),
                             'pB_S1' = parse(text = TeX("$\\p_{B_S_{1}}$")), 'pB_S2' = parse(text = TeX("$\\p_{B_S_{2}}$")), 'pB_S3' = parse(text = TeX("$\\p_{B_S{3}}$")),
                             'pB_W1' = parse(text = TeX("$\\p_{B_W_{1}}$")), 'pB_W2' = parse(text = TeX("$\\p_{B_W_{2}}$")), 'pB_W3' = parse(text = TeX("$\\p_{B_W{3}}$")),
                             'pA_S1T' = parse(text = TeX("$\\p_{A_S_{1T}}$")), 'pA_S2T' = parse(text = TeX("$\\p_{A_S_{2T}}$")), 'pA_S3T' = parse(text = TeX("$\\p_{A_S_{3T}}$")),
                             'pA_W1T' = parse(text = TeX("$\\p_{A_W_{1T}}$")), 'pA_W2T' = parse(text = TeX("$\\p_{A_W_{2T}}$")), 'pA_W3T' = parse(text = TeX("$\\p_{A_W_{3T}}$")),
                             'pB_S1T' = parse(text = TeX("$\\p_{B_S_{1T}}$")), 'pB_S2T' = parse(text = TeX("$\\p_{B_S_{2T}}$")), 'pB_S3T' = parse(text = TeX("$\\p_{B_S_{3T}}$")),
                             'pB_W1T' = parse(text = TeX("$\\p_{B_W_{1T}}$")), 'pB_W2T' = parse(text = TeX("$\\p_{B_W_{2T}}$")), 'pB_W3T' = parse(text = TeX("$\\p_{B_W_{3T}}$")) ))
                             
                            
ggsave("modperf_p_pripost.png", width = 60, height = 30, units="cm")
#~ End of Script