######################################################################################################################################################
##  Plots comparing model estimates from the A dynamic occupancy model for interacting species with two spatial scales with true parameters values   ##
##              by EFK                                                                                                                              ##
######################################################################################################################################################

rm(list=ls())

# Call jags(and other packages)
library(jagsUI)
library(ggplot2)
library(latex2exp)

# set working directory
#setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models")   

setwd("./simulation_study/model_output")

#load model outputs
load("mod_4stpm_sim_low_det_2.rda")
load("mod_4stpm_sim_mid_det_2.rda")
load("mod_4stpm_sim_hig_det.rda")
load("mod_4stpm_sim_low_occ_2.rda")
load("mod_4stpm_sim_mid_occ_2.rda")
load("mod_4stpm_sim_hig_occ.rda")

####
# checking convergence
par(mfrow=c(1,1))
hist(mod_4stpm_sim_low_det[[1]]$Rhat$x)

#trace plots
par(mfrow=c(2,4))
traceplot(mod_4stpm_sim_low_det[[1]], parameters = c("gamA","gamAB","gamB","gamBA","epsA","epsAB","epsB","epsBA"))
traceplot(mod_4stpm_sim_low_det[[1]], parameters = c("GamA","GamAB","GamB","GamBA","EpsA","EpsAB","EpsB","EpsBA"))

traceplot(mod_4stpm_sim_hig_det[[1]], parameters = c("gamA","gamAB","gamB","gamBA","epsA","epsAB","epsB","epsBA"))
traceplot(mod_4stpm_sim_hig_det[[1]], parameters = c("GamA","GamAB","GamB","GamBA","EpsA","EpsAB","EpsB","EpsBA"))

# check some simulations
mod_4stpm_sim_mid_det[[1]]$mean$EpsAB
mod_4stpm_sim_mid_det[[2]]$mean$EpsAB
mod_4stpm_sim_mid_det[[3]]$mean$EpsAB
mod_4stpm_sim_mid_det[[4]]$mean$EpsAB

mod_4stpm_sim_hig_det[[1]]$mean$EpsAB
mod_4stpm_sim_hig_det[[2]]$mean$EpsAB
mod_4stpm_sim_hig_det[[3]]$mean$EpsAB
mod_4stpm_sim_hig_det[[4]]$mean$EpsAB

###################################
### violin plot with ggplot2

# make a array with estimates from all model runs for all gammas and epsilons
dat<-c()

# make a vector of the colonization and extinction estimates from all 50 simulations
for (i in 1:50){
  dat[i]     <- mod_4stpm_sim_low_det[[i]]$mean$GamA
  dat[i+50]  <- mod_4stpm_sim_mid_det[[i]]$mean$GamA 
  dat[i+100] <- mod_4stpm_sim_hig_det[[i]]$mean$GamA
  dat[i+150] <- mod_4stpm_sim_low_occ[[i]]$mean$GamA
  dat[i+200] <- mod_4stpm_sim_mid_occ[[i]]$mean$GamA
  dat[i+250] <- mod_4stpm_sim_hig_occ[[i]]$mean$GamA
  
  dat[i+300] <- mod_4stpm_sim_low_det[[i]]$mean$GamB
  dat[i+350] <- mod_4stpm_sim_mid_det[[i]]$mean$GamB
  dat[i+400] <- mod_4stpm_sim_hig_det[[i]]$mean$GamB
  dat[i+450] <- mod_4stpm_sim_low_occ[[i]]$mean$GamB
  dat[i+500] <- mod_4stpm_sim_mid_occ[[i]]$mean$GamB
  dat[i+550] <- mod_4stpm_sim_hig_occ[[i]]$mean$GamB
  
  dat[i+600] <- mod_4stpm_sim_low_det[[i]]$mean$GamAB
  dat[i+650] <- mod_4stpm_sim_mid_det[[i]]$mean$GamAB
  dat[i+700] <- mod_4stpm_sim_hig_det[[i]]$mean$GamAB
  dat[i+750] <- mod_4stpm_sim_low_occ[[i]]$mean$GamAB
  dat[i+800] <- mod_4stpm_sim_mid_occ[[i]]$mean$GamAB
  dat[i+850] <- mod_4stpm_sim_hig_occ[[i]]$mean$GamAB
  
  dat[i+900] <- mod_4stpm_sim_low_det[[i]]$mean$GamBA
  dat[i+950] <- mod_4stpm_sim_mid_det[[i]]$mean$GamBA
  dat[i+1000] <- mod_4stpm_sim_hig_det[[i]]$mean$GamBA
  dat[i+1050] <- mod_4stpm_sim_low_occ[[i]]$mean$GamBA
  dat[i+1100] <- mod_4stpm_sim_mid_occ[[i]]$mean$GamBA
  dat[i+1150] <- mod_4stpm_sim_hig_occ[[i]]$mean$GamBA
  
  dat[i+1200] <- mod_4stpm_sim_low_det[[i]]$mean$EpsA
  dat[i+1250] <- mod_4stpm_sim_mid_det[[i]]$mean$EpsA
  dat[i+1300] <- mod_4stpm_sim_hig_det[[i]]$mean$EpsA
  dat[i+1350] <- mod_4stpm_sim_low_occ[[i]]$mean$EpsA
  dat[i+1400] <- mod_4stpm_sim_mid_occ[[i]]$mean$EpsA
  dat[i+1450] <- mod_4stpm_sim_hig_occ[[i]]$mean$EpsA
  
  dat[i+1500] <- mod_4stpm_sim_low_det[[i]]$mean$EpsB
  dat[i+1550] <- mod_4stpm_sim_mid_det[[i]]$mean$EpsB
  dat[i+1600] <- mod_4stpm_sim_hig_det[[i]]$mean$EpsB
  dat[i+1650] <- mod_4stpm_sim_low_occ[[i]]$mean$EpsB
  dat[i+1700] <- mod_4stpm_sim_mid_occ[[i]]$mean$EpsB
  dat[i+1750] <- mod_4stpm_sim_hig_occ[[i]]$mean$EpsB
  
  dat[i+1800] <- mod_4stpm_sim_low_det[[i]]$mean$EpsAB
  dat[i+1850] <- mod_4stpm_sim_mid_det[[i]]$mean$EpsAB
  dat[i+1900] <- mod_4stpm_sim_hig_det[[i]]$mean$EpsAB
  dat[i+1950] <- mod_4stpm_sim_low_occ[[i]]$mean$EpsAB
  dat[i+2000] <- mod_4stpm_sim_mid_occ[[i]]$mean$EpsAB
  dat[i+2050] <- mod_4stpm_sim_hig_occ[[i]]$mean$EpsAB
  
  dat[i+2100] <- mod_4stpm_sim_low_det[[i]]$mean$EpsBA
  dat[i+2150] <- mod_4stpm_sim_mid_det[[i]]$mean$EpsBA
  dat[i+2200] <- mod_4stpm_sim_hig_det[[i]]$mean$EpsBA
  dat[i+2250] <- mod_4stpm_sim_low_occ[[i]]$mean$EpsBA
  dat[i+2300] <- mod_4stpm_sim_mid_occ[[i]]$mean$EpsBA
  dat[i+2350] <- mod_4stpm_sim_hig_occ[[i]]$mean$EpsBA
  
  dat[i+2400] <- mod_4stpm_sim_low_det[[i]]$mean$gamA
  dat[i+2450] <- mod_4stpm_sim_mid_det[[i]]$mean$gamA
  dat[i+2500] <- mod_4stpm_sim_hig_det[[i]]$mean$gamA
  dat[i+2550] <- mod_4stpm_sim_low_occ[[i]]$mean$gamA
  dat[i+2600] <- mod_4stpm_sim_mid_occ[[i]]$mean$gamA
  dat[i+2650] <- mod_4stpm_sim_hig_occ[[i]]$mean$gamA
  
  dat[i+2700] <- mod_4stpm_sim_low_det[[i]]$mean$gamB
  dat[i+2750] <- mod_4stpm_sim_mid_det[[i]]$mean$gamB
  dat[i+2800] <- mod_4stpm_sim_hig_det[[i]]$mean$gamB
  dat[i+2850] <- mod_4stpm_sim_low_occ[[i]]$mean$gamB
  dat[i+2900] <- mod_4stpm_sim_mid_occ[[i]]$mean$gamB
  dat[i+2950] <- mod_4stpm_sim_hig_occ[[i]]$mean$gamB
  
  dat[i+3000] <- mod_4stpm_sim_low_det[[i]]$mean$gamAB
  dat[i+3050] <- mod_4stpm_sim_mid_det[[i]]$mean$gamAB
  dat[i+3100] <- mod_4stpm_sim_hig_det[[i]]$mean$gamAB
  dat[i+3150] <- mod_4stpm_sim_low_occ[[i]]$mean$gamAB
  dat[i+3200] <- mod_4stpm_sim_mid_occ[[i]]$mean$gamAB
  dat[i+3250] <- mod_4stpm_sim_hig_occ[[i]]$mean$gamAB
  
  dat[i+3300] <- mod_4stpm_sim_low_det[[i]]$mean$gamBA
  dat[i+3350] <- mod_4stpm_sim_mid_det[[i]]$mean$gamBA
  dat[i+3400] <- mod_4stpm_sim_hig_det[[i]]$mean$gamBA
  dat[i+3450] <- mod_4stpm_sim_low_occ[[i]]$mean$gamBA
  dat[i+3500] <- mod_4stpm_sim_mid_occ[[i]]$mean$gamBA
  dat[i+3550] <- mod_4stpm_sim_hig_occ[[i]]$mean$gamBA
  
  dat[i+3600] <- mod_4stpm_sim_low_det[[i]]$mean$epsA
  dat[i+3650] <- mod_4stpm_sim_mid_det[[i]]$mean$epsA
  dat[i+3700] <- mod_4stpm_sim_hig_det[[i]]$mean$epsA
  dat[i+3750] <- mod_4stpm_sim_low_occ[[i]]$mean$epsA
  dat[i+3800] <- mod_4stpm_sim_mid_occ[[i]]$mean$epsA
  dat[i+3850] <- mod_4stpm_sim_hig_occ[[i]]$mean$epsA
  
  dat[i+3900] <- mod_4stpm_sim_low_det[[i]]$mean$epsB
  dat[i+3950] <- mod_4stpm_sim_mid_det[[i]]$mean$epsB
  dat[i+4000] <- mod_4stpm_sim_hig_det[[i]]$mean$epsB
  dat[i+4050] <- mod_4stpm_sim_low_occ[[i]]$mean$epsB
  dat[i+4100] <- mod_4stpm_sim_mid_occ[[i]]$mean$epsB
  dat[i+4150] <- mod_4stpm_sim_hig_occ[[i]]$mean$epsB
  
  dat[i+4200] <- mod_4stpm_sim_low_det[[i]]$mean$epsAB
  dat[i+4250] <- mod_4stpm_sim_mid_det[[i]]$mean$epsAB
  dat[i+4300] <- mod_4stpm_sim_hig_det[[i]]$mean$epsAB
  dat[i+4350] <- mod_4stpm_sim_low_occ[[i]]$mean$epsAB
  dat[i+4400] <- mod_4stpm_sim_mid_occ[[i]]$mean$epsAB
  dat[i+4450] <- mod_4stpm_sim_hig_occ[[i]]$mean$epsAB
  
  dat[i+4500] <- mod_4stpm_sim_low_det[[i]]$mean$epsBA
  dat[i+4550] <- mod_4stpm_sim_mid_det[[i]]$mean$epsBA
  dat[i+4600] <- mod_4stpm_sim_hig_det[[i]]$mean$epsBA
  dat[i+4650] <- mod_4stpm_sim_low_occ[[i]]$mean$epsBA
  dat[i+4700] <- mod_4stpm_sim_mid_occ[[i]]$mean$epsBA
  dat[i+4750] <- mod_4stpm_sim_hig_occ[[i]]$mean$epsBA
}

dat <- as.data.frame(dat) # make into data.frame to be able to add more columns
"A \U2223 B"

# add a column describing the parameter name with greek letters
dat$param <- factor(c(rep("GamA", times=300),rep("GamB", times=300),rep("GamAB", times=300),rep("GamBA", times=300),
                      rep("EpsA", times=300), rep("EpsB", times=300),rep("EpsAB", times=300),rep("EpsBA", times=300),
                      rep("gamA", times=300),rep("gamB", times=300),rep("gamAB", times=300),rep("gamBA", times=300),
                      rep("epsA", times=300), rep("epsB", times=300),rep("epsAB", times=300),rep("epsBA", times=300)),
                    levels=c("GamA","GamB","GamAB","GamBA",
                             "EpsA","EpsB","EpsAB","EpsBA",
                             "gamA","gamB","gamAB","gamBA",
                             "epsA","epsB","epsAB","epsBA"))

# adding a column to describe what set of true parameter values was used
dat$sim <- as.factor(c(rep("ld", times=50),rep("md", times=50),rep("hd", times=50), rep("lo", times=50),rep("mo", times=50),rep("ho", times=50)))

# make data.frame of true parameter values
true <- c()
true[1:16] <- c(0.5, 0.1, 0.05, 0.4, 0.05, 0.6, 0.5, 0.2, 0.5, 0.3, 0.1, 0.7, 0.3, 0.8, 0.9, 0.1)
true[17:32] <- c(0.5, 0.1, 0.05, 0.4, 0.05, 0.6, 0.5, 0.2, 0.5, 0.3, 0.1, 0.7, 0.3, 0.8, 0.9, 0.1)
true[33:48] <- c(0.5, 0.1, 0.05, 0.4, 0.05, 0.6, 0.5, 0.2, 0.5, 0.3, 0.1, 0.7, 0.3, 0.8, 0.9, 0.1)
true[49:64] <- c(0.1, 0.05,0.05, 0.2, 0.1,  0.6, 0.5, 0.2, 0.3, 0.3, 0.1, 0.6, 0.3, 0.8, 0.9, 0.1)
true[65:80] <- c(0.5, 0.1, 0.05, 0.4, 0.05, 0.6, 0.5, 0.2, 0.5, 0.3, 0.1, 0.7, 0.3, 0.8, 0.9, 0.1)
true[81:96] <- c(0.8, 0.2, 0.2,  0.7, 0.05, 0.8, 0.4, 0.2, 0.8, 0.3, 0.1, 0.7, 0.1, 0.6, 0.6, 0.1)

true <- as.data.frame(true) # turn into df to be able to add columns

# estimate mean of simulation means
avr <- aggregate(dat$dat, by=list(dat$param, dat$sim), mean)
std <- aggregate(dat$dat, by=list(dat$param, dat$sim), sd)

#<- "A\u007CB"
"A_\u007C_B"
# add a column describing the parameter name with greek letters
true$param <- factor(rep(c("GamA","GamB","GamAB","GamBA",
                           "EpsA","EpsB","EpsAB","EpsBA",
                           "gamA","gamB","gamAB","gamBA",
                           "epsA","epsB","epsAB","epsBA"), times=6), 
                     levels=c("GamA","GamB","GamAB","GamBA",
                              "EpsA","EpsB","EpsAB","EpsBA",
                              "gamA","gamB","gamAB","gamBA",
                              "epsA","epsB","epsAB","epsBA") )

# adding a column to describe what set of true parameter values was used                 
true$sim <- as.factor(c(rep("ld", times=16), rep("md", times=16), rep("hd", times=16), rep("lo", times=16), rep("mo", times=16), rep("ho", times=16)))

names(true)[1] <- "dat"                       # change col name to match the other data frame(dat)

# change col names for avg to match the others
names(avr) <- c(names(true)[2],names(true)[3],names(true)[1])

positions <-c("ld","md","hd", "lo","mo","ho") # spesify what order to plot the simulations

par_names <- list('GamA' = parse(text = TeX("$\\Gamma_{A}$")), 'GamB' = parse(text = TeX("$\\Gamma_{B}$")),
                  'GamAB' = parse(text = TeX("$\\Gamma_{A|B}$")), 'GamBA' = parse(text = TeX("$\\Gamma_{B|A}$")),
                  'EpsA' = parse(text = TeX("$E_{A}$")), 'EpsB' = parse(text = TeX("$E_{B}$")), 
                  'EpsAB' = parse(text = TeX("$E_{A|B}$")), 'EpsBA' = parse(text = TeX("$E_{B|A}$")),
                  'gamA' = parse(text = TeX("$\\gamma_{A}$")), 'gamB' = parse(text = TeX("$\\gamma_{B}$")), 
                  'gamAB' = parse(text = TeX("$\\gamma_{A|B}$")), 'gamBA' = parse(text = TeX("$\\gamma_{B|A}$")),
                  'epsA' = parse(text = TeX("$\\epsilon_{A}$")), 'epsB' = parse(text = TeX("$\\epsilon_{B}$")),
                  'epsAB' = parse(text = TeX("$\\epsilon_{A|B}$")), 'epsBA' = parse(text = TeX("$\\epsilon_{B|A}$"))
 )

par_labeller <- function(variable,value){
  return(par_names[value])
}

positions <-c("ld","md","hd", "lo","mo","ho")
               
ggplot(data=dat, aes(x=sim, y=dat, fill="grey"))+
  geom_violin(fill="grey")+
  geom_boxplot(width=0.4, color="black", alpha=0.4, fill="white")+
  geom_point(data=avr, color="grey40", shape="-", size=20)+
  geom_point(data=true, color="red",  shape="-", size=20, alpha=0.8)+
  labs(y="", x="")+
  facet_wrap(~param, labeller = par_labeller,  ncol=4)+
  theme(strip.text = element_text(size=12,face="bold"), 
        axis.text.x = element_text( size = 12, face="bold"), 
        axis.text.y = element_text( size = 12, face="bold"),
        legend.position = "none")+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  scale_x_discrete(limits = positions)


# save plot
setwd("../plot")
#setwd("~/UiT/GitProjects/A-dynamic-and-hierarchical-spatial-occupancy-model-for-interacting-species/simulation_study/plot")
ggsave( "modelperformance_4stpm_sim_gameps_4.png", width = 24, height = 20,units="cm", dpi = 300)


###########################################################
## violin plot for detection probabilities (pA and pB)  ###
###########################################################

dat_p <- c()

# make a vector of the detection probability estimates from all 50 simulations
for (i in 1:50){
  dat_p[i]     <- mod_4stpm_sim_low_det[[i]]$mean$pA
  dat_p[i+50]  <- mod_4stpm_sim_mid_det[[i]]$mean$pA 
  dat_p[i+100] <- mod_4stpm_sim_hig_det[[i]]$mean$pA
  dat_p[i+150] <- mod_4stpm_sim_low_occ[[i]]$mean$pA
  dat_p[i+200] <- mod_4stpm_sim_mid_occ[[i]]$mean$pA
  dat_p[i+250] <- mod_4stpm_sim_hig_occ[[i]]$mean$pA

  dat_p[i+300] <- mod_4stpm_sim_low_det[[i]]$mean$pB
  dat_p[i+350] <- mod_4stpm_sim_mid_det[[i]]$mean$pB
  dat_p[i+400] <- mod_4stpm_sim_hig_det[[i]]$mean$pB
  dat_p[i+450] <- mod_4stpm_sim_low_occ[[i]]$mean$pB
  dat_p[i+500] <- mod_4stpm_sim_mid_occ[[i]]$mean$pB
  dat_p[i+550] <- mod_4stpm_sim_hig_occ[[i]]$mean$pB
  }

dat_p <- as.data.frame(dat_p) # turn into df to add columns

# add a column with the parameter name
dat_p$param <- c(rep("pA", times=300),rep("pB", times=300))

#addin a column describing what true parameter set that was used in the simulation
dat_p$sim <- rep(c(rep("ld", times=50), rep("md", times=50), rep("hd", times=50), rep("lo", times=50), rep("mo", times=50), rep("ho", times=50)),times=2)

# spesifying true parameter values
true_p <- c(0.2,0.1,0.5,0.5,0.9,0.8,0.5,0.5,0.5,0.5,0.5,0.5)

true_p <- as.data.frame(true_p) # turn into df to add columns
true_p$param <- rep(c("pA","pB"), times=6) # add parameter names
true_p$sim <- c("ld","ld","md","md","hd","hd","lo","lo","mo","mo","ho","ho") # add simulation ID

names(true_p)[1] <- "dat_p"                       # change col name to match the other data frame(dat)
positions <-c("ld","md","hd", "lo","mo","ho") # spesify what order to plot the simulations

ggplot(data=dat_p, aes(x=sim, y=dat_p, fill="grey"))+
  geom_violin(fill="grey")+
  geom_boxplot(width=0.4, color="black", alpha=0.4, fill="white")+
  geom_point(data=true_p, color="red", shape="-", size=12)+
  labs(y="", x="")+
  facet_wrap(~param, ncol=2)+
  theme(strip.text = element_text(size=12,face="bold"), 
        axis.text.x = element_text( size = 12, face="bold"), 
        axis.text.y = element_text( size = 12, face="bold"),
        legend.position = "none")+
  scale_x_discrete(limits = positions)

ggsave("modelperformance_4stpm_sim_p_2.png")

###########################################################
## violin plot for initial occupancy probabilities (psi)  ###
###########################################################

dat_psi <- c()
for (i in 1:50){
  dat_psi[i]     <- mean(mod_4stpm_sim_low_det[[i]]$mean$psi[,1])
  dat_psi[i+50]  <- mean(mod_4stpm_sim_mid_det[[i]]$mean$psi[,1]) 
  dat_psi[i+100] <- mean(mod_4stpm_sim_hig_det[[i]]$mean$psi[,1])
  dat_psi[i+150] <- mean(mod_4stpm_sim_low_occ[[i]]$mean$psi[,1])
  dat_psi[i+200] <- mean(mod_4stpm_sim_mid_occ[[i]]$mean$psi[,1])
  dat_psi[i+250] <- mean(mod_4stpm_sim_hig_occ[[i]]$mean$psi[,1])
  
  dat_psi[i+300] <- mean(mod_4stpm_sim_low_det[[i]]$mean$psi[,2])
  dat_psi[i+350] <- mean(mod_4stpm_sim_mid_det[[i]]$mean$psi[,2]) 
  dat_psi[i+400] <- mean(mod_4stpm_sim_hig_det[[i]]$mean$psi[,2])
  dat_psi[i+450] <- mean(mod_4stpm_sim_low_occ[[i]]$mean$psi[,2])
  dat_psi[i+500] <- mean(mod_4stpm_sim_mid_occ[[i]]$mean$psi[,2])
  dat_psi[i+550] <- mean(mod_4stpm_sim_hig_occ[[i]]$mean$psi[,2])
  
  dat_psi[i+600] <- mean(mod_4stpm_sim_low_det[[i]]$mean$psi[,3])
  dat_psi[i+650] <- mean(mod_4stpm_sim_mid_det[[i]]$mean$psi[,3])
  dat_psi[i+700] <- mean(mod_4stpm_sim_hig_det[[i]]$mean$psi[,3])
  dat_psi[i+750] <- mean(mod_4stpm_sim_low_occ[[i]]$mean$psi[,3])
  dat_psi[i+800] <- mean(mod_4stpm_sim_mid_occ[[i]]$mean$psi[,3])
  dat_psi[i+850] <- mean(mod_4stpm_sim_hig_occ[[i]]$mean$psi[,3])
}

dat_psi <- as.data.frame(dat_psi)

dat_psi$param <- c(rep("\u03C8[A]", times=300),rep("\u03C8[B]", times=300),rep("\u03C8[AB]", times=300))
dat_psi$sim <- rep(c(rep("ld", times=50), rep("md", times=50), rep("hd", times=50), rep("lo", times=50), rep("mo", times=50), rep("ho", times=50)),times=3)

true_psi <- c(0.25, 0.15, 0.1,  0.25, 0.15, 0.1,  0.25, 0.15, 0.1,
              0.2, 0.15, 0.05,  0.25, 0.15, 0.1,  0.3, 0.2, 0.1)

true_psi <- as.data.frame(true_psi)

true_psi$param <- rep(c("\u03C8[A]","\u03C8[B]","\u03C8[AB]"), times=6)
true_psi$sim <- c("ld","ld","ld","md","md","md","hd","hd","hd","lo","lo","lo","mo","mo","mo","ho","ho","ho")

names(true_psi)[1] <- "dat_psi"                   # change col name to match the other data frame(dat)
positions <-c("ld","md","hd", "lo","mo","ho")     # spesify the order the simulations will be plotted in

ggplot(data=dat_psi, aes(x=sim, y=dat_psi, fill="grey"))+
  geom_violin(fill="grey")+
  geom_boxplot(width=0.4, color="black", alpha=0.4, fill="white")+
  geom_point(data=true_psi, color="red", shape="-", size=12)+
  labs(y="", x="")+
  facet_wrap(~param, ncol=3, labeller = label_parsed)+
  theme(strip.text = element_text(size=12,face="bold"), 
        axis.text.x = element_text( size = 12, face="bold"), 
        axis.text.y = element_text( size = 12, face="bold"),
        legend.position = "none")+
  scale_x_discrete(limits = positions)

ggsave("modelperformance_4stpm_sim_psi_2.png")

#~ end of script


str(mod_4stpm_sim_low_det)
str(dat)
library(dplyr)

sub <- filter(dat, param == "Ε[B]", sim == "ld") 

mean(sub$dat)# mean of simulation means = 0.72
true         # true EB under ld = 0.6

0.1288/0.6 # 21% bias


sub2 <- filter(dat, param == "Ε[AB]", sim == "ld") 

mean(sub2$dat)# mean of simulation means = 0.6
true         # true E_AB under ld = 0.5

0.1/0.5 # 20% bias

unique(dat$param)
