#set wd
getwd()

setwd("../../case_study/model_output")
setwd("./simulation_study/model_output")
setwd("../model_output")

#import model
dir()
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop_ni10k_pri1.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop_ni10k_pri2.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop_ni10k_pri2_4.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop_ni10k_pri2_5.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop_ni10k_pri2_6.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop11_ni500_pri1.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop_ni10k_pri3_2.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop11_ni500_pri1_2.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop11_ni500_pri1_3.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop11_ni500_pri1_3.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop11_ni500_pri1_3e.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop11_ni500_pri1_4.rda")
load("va_snowbed_mustela_rodent_outgof_jagui_nestedloop12_ni1k_pri1.rda")

# traceplots and convergence check
traceplot(out.gof, parameters=c("pA", "pB"))
print(out.gof, 2)

hist(out.gof$Rhat$z)
max(out.gof$Rhat$z, na.rm=T)

# Plots of expected versus observed value of fit stats
# Open part
# species A
par(mfrow=c(2,2))
pl <- range(c(out.gof$sims.list$Chi2Open_A, out.gof$sims.list$Chi2repOpen_A))
plot(out.gof$sims.list$Chi2Open_A, out.gof$sims.list$Chi2repOpen_A,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Open part of model Species A ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(450, 1430, paste('Bpv = ', round(mean(out.gof$sims.list$Chi2repOpen_A >
                                             out.gof$sims.list$Chi2Open_A), 2)), cex = 2)

# species B
pl <- range(c(out.gof$sims.list$Chi2Open_B, out.gof$sims.list$Chi2repOpen_B))
plot(out.gof$sims.list$Chi2Open_B, out.gof$sims.list$Chi2repOpen_B,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Open part of model Species B ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(450, 1430, paste('Bpv = ', round(mean(out.gof$sims.list$Chi2repOpen_B >
                                             out.gof$sims.list$Chi2Open_B), 2)), cex = 2)

# Closed part of model: Chi-squared
pl <- range(c(out.gof$sims.list$Chi2Closed_A, out.gof$sims.list$Chi2repClosed_A))
plot(out.gof$sims.list$Chi2Closed_A, out.gof$sims.list$Chi2repClosed_A,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Closed part of model Species A", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(525, 720, paste('Bpv = ', round(mean(out.gof$sims.list$Chi2repClosed_A >
                                            out.gof$sims.list$Chi2Closed_A), 2)), cex = 2)

# Closed part of model: Chi-squared
pl <- range(c(out.gof$sims.list$Chi2Closed_B, out.gof$sims.list$Chi2repClosed_B))
plot(out.gof$sims.list$Chi2Closed_B, out.gof$sims.list$Chi2repClosed_B,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Closed part of model Species B", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(525, 720, paste('Bpv = ', round(mean(out.gof$sims.list$Chi2repClosed_B >
                                            out.gof$sims.list$Chi2Closed_B), 2)), cex = 2)


# Closed part of model: Freeman-Tukey
pl <- range(c(out.gof$sims.list$FTClosed_B, out.gof$sims.list$FTrepClosed_B))
plot(out.gof$sims.list$FTClosed_B, out.gof$sims.list$FTrepClosed_B,
     xlab = "FT observed data", ylab = "FT expected data",
     main = "Closed part of model (Freeman-Tukey) Species B", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(240, 335, paste('Bpv = ', round(mean(out.gof$sims.list$FTrepClosed_B >
                                            out.gof$sims.list$FTClosed_B), 2)), cex = 2)

pl <- range(c(out.gof$sims.list$FTClosed_A, out.gof$sims.list$FTrepClosed_A))
plot(out.gof$sims.list$FTClosed_A, out.gof$sims.list$FTrepClosed_A,
     xlab = "FT observed data", ylab = "FT expected data",
     main = "Closed part of model (Freeman-Tukey) Species A", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(240, 335, paste('Bpv = ', round(mean(out.gof$sims.list$FTrepClosed_B >
                                            out.gof$sims.list$FTClosed_B), 2)), cex = 2)

#par(op)

# Inspect observed and expected numbers of transitions for each interval
obs.trans <- sum(out.gof$mean$tm_A[])
rep.trans <- out.gof$mean$tmrep_A
exp.trans <- out.gof$mean$Etm_A
dimnames(obs.trans) <- dimnames(rep.trans) <- dimnames(exp.trans) <-  list(c("From Non-occ", "From Occ"), c("To Non-occ", "To Occ"))

for(t in 1:(data$nseason-1)){
  cat("\n*** Interval", t, " ***\n")
  cat("* Observed transitions*\n")
  print(obs.trans[,,t])
  
  cat("\n* Replicate data transitions*\n")
  print(rep.trans[,,t])
  
  cat("\n* Expected transitions*\n")
  print(exp.trans[,,t])
}

table(out.gof$mean$ext_A)
table(out.gof$mean$nonext_A)
table(out.gof$mean$colo_A)
table(out.gof$mean$noncolo_A)

dim(out.gof$mean$ext_A)

sum(out.gof$mean$ext_A[,1])

sum(out.gof$mean$tm_A[1,1,])
sum(out.gof$mean$Etm_A[1,1,])

t=dim(out.gof$mean$tm_A)[3]

dat1 <- c()
dat2 <- c()
dat3 <- c()
dat4 <- c()

e = 0.01

for(i in 1:t){
 dat1[i] <- (out.gof$mean$tm_A[1,1,i] - out.gof$mean$Etm_A[1,1,i])^2 / (out.gof$mean$Etm_A[1,1,i] + e)
 dat2[i] <- (out.gof$mean$tm_A[1,2,i] - out.gof$mean$Etm_A[1,2,i])^2 / (out.gof$mean$Etm_A[1,2,i] + e)
 dat3[i] <- (out.gof$mean$tm_A[2,1,i] - out.gof$mean$Etm_A[2,1,i])^2 / (out.gof$mean$Etm_A[2,1,i] + e)
 dat4[i] <- (out.gof$mean$tm_A[2,2,i] - out.gof$mean$Etm_A[2,2,i])^2 / (out.gof$mean$Etm_A[2,2,i] + e)
}

sum(dat1)
sum(dat2)
sum(dat3)
sum(dat4)

sum(dat1,dat2,dat3,dat4)

################################################

sum(out.gof$mean$tmp_A[6,,1])
out.gof$mean$tmp_A[6,,1]

table(out.gof$mean$E_A)

out.gof$mean$pB
out.gof$mean$pA


table(out.gof$mean$E_B)



str(out.gof$mean$detfreq_A)
str(out.gof$mean$z_A)

out.gof$mean$detfreq_A[1,1:10]
out.gof$mean$z_A[1,1:10]


dim(out.gof$mean$tmp_A)
sum(out.gof$mean$tmp_A[1,,1])
out.gof$mean$detfreq_A[1,1]

sum(out.gof$mean$tmp_A[6,,1])
out.gof$mean$detfreq_A[6,1]

detfreq_A <- out.gof$mean$detfreq_A
E_A <- out.gof$mean$E_A
e = 0.0001

dat<-c()

for(i in 1:dim(detfreq_A)[1]){
dat[i] <- (detfreq_A[i,2] - E_A[i,2])^2 / (E_A[i,2] + e)
}

sum(dat)

detfreqrep_A <- out.gof$mean$detfreqrep_A
E_A <- out.gof$mean$E_A
dat2<-c()

for(i in 1:dim(detfreqrep_A)[1]){
  dat2[i] <- (detfreqrep_A[i,2] - E_A[i,2])^2 / (E_A[i,2] + e)
}

sum(dat2)



# how many mustelids are observed during winter?
getwd()
setwd("../data")

library(dplyr)
library(tibble)

dir()
load("occm_mustela_rodent_var_snowbed_rmBQ.rda")
load("season.rda")

dim(occm_ko3)


occ_w <- occm_ko3[,,c(which(season==1)),]
occ_s <- occm_ko3[,,c(which(season==0)[-length(which(season==0))]),]

table(occ_w)
table(occ_s)
