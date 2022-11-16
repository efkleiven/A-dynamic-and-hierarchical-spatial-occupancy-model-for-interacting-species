#set wd
getwd()

#setwd("../../case_study/model_output")
setwd("../../model_output/revision2")

#import model
dir()

load("va_snowbed_mustela_rodent_gof_nestedloop_temp_pri1_005.rda")
load("va_snowbed_mustela_rodent_gof_nestedloop_temp_pri2_005.rda")
load("va_snowbed_mustela_rodent_gof_nestedloop_temp_pri3_005.rda")

# traceplots and convergence check
#traceplot(out.gof, parameters=c("pA", "pB"))
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots of expected versus observed value of fit stats
# Open part
# species A

setwd("../plot")

pdf(file="GOF_casestudy.pdf")
par(mfrow=c(2,2))
pl <- range(c(out.gof$sims.list$Chi2Open_A, out.gof$sims.list$Chi2repOpen_A))
pl[1]<- 0
pl[2]<- 200000

plot(out.gof$sims.list$Chi2Open_A, out.gof$sims.list$Chi2repOpen_A,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Open part of model Species A ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(50000, 180000, paste('Bpv = ', round(mean(out.gof$sims.list$Chi2repOpen_A >
                                             out.gof$sims.list$Chi2Open_A), 2)), cex = 1.5)

# species B
pl <- range(c(out.gof$sims.list$Chi2Open_B, out.gof$sims.list$Chi2repOpen_B))
plot(out.gof$sims.list$Chi2Open_B, out.gof$sims.list$Chi2repOpen_B,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Open part of model Species B ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(300000, 20000,  paste('Bpv = ', round(mean(out.gof$sims.list$Chi2repOpen_B >
                                             out.gof$sims.list$Chi2Open_B), 3)), cex = 1.5)

# Closed part of model: Chi-squared
pl <- range(c(out.gof$sims.list$Chi2Closed_A, out.gof$sims.list$Chi2repClosed_A))
plot(out.gof$sims.list$Chi2Closed_A, out.gof$sims.list$Chi2repClosed_A,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Closed part of model Species A", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(250000,1000000, paste('Bpv = ', round(mean(out.gof$sims.list$Chi2repClosed_A >
                                            out.gof$sims.list$Chi2Closed_A), 2)), cex = 1.5)

# Closed part of model: Chi-squared
pl <- range(c(out.gof$sims.list$Chi2Closed_B, out.gof$sims.list$Chi2repClosed_B))
plot(out.gof$sims.list$Chi2Closed_B, out.gof$sims.list$Chi2repClosed_B,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Closed part of model Species B", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(350000, 80000, paste('Bpv = ', round(mean(out.gof$sims.list$Chi2repClosed_B >
                                            out.gof$sims.list$Chi2Closed_B), 2)), cex = 1.5)
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# end of script

