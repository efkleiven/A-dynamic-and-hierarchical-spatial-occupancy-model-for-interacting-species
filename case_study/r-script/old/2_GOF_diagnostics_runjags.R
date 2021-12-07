#set wd
getwd()
setwd("./model_output")

#import model
dir()
load("va_snowbed_mustela_rodent_outgof_rj_ni5k.rda")

# add summary
# add.summary(out.gof)
# traceplot(out.gof)
print(out.gof, 2)

out.gof$mcmc[[1]][,1]
tail(names(out.gof$mcmc[[1]][1,]))

which(names(out.gof$mcmc[[1]][1,]) == "Chi2Open_A")
which(names(out.gof$mcmc[[1]][1,]) == "Chi2repOpen_A")

class(out.gof[[1]][,11835][[1]])

as.matrix(out.gof[[1]][,11835][[1]])
extract()

# Plots of expected versus observed value of fit stats
# Open part
# species A
pl <- range(c(as.matrix(out.gof[[1]][,11835][[2]]), as.matrix(out.gof[[1]][,11836][[2]])))
plot(as.matrix(out.gof[[1]][,11835][[2]]), as.matrix(out.gof[[1]][,11836][[2]]),
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Open part of model ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(450, 1430, paste('Bpv = ', round(mean(as.matrix(out.gof[[1]][,11836][[1]]) >
                                             as.matrix(out.gof[[1]][,11835][[1]])), 2)), cex = 2)


# species B

which(names(out.gof$mcmc[[1]][1,]) == "Chi2Open_B")
which(names(out.gof$mcmc[[1]][1,]) == "Chi2repOpen_B")

pl <- range(c(as.matrix(out.gof[[1]][,14265][[1]]), as.matrix(out.gof[[1]][,14266][[1]])))
plot(as.matrix(out.gof[[1]][,14265][[1]]), as.matrix(out.gof[[1]][,14266][[1]]),
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Open part of model ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(450, 1430, paste('Bpv = ', round(mean(as.matrix(out.gof[[1]][,14266][[1]]) >
                                             as.matrix(out.gof[[1]][,14265][[1]])), 2)), cex = 2)

# Closed part of model: Chi-squared

which(names(out.gof$mcmc[[1]][1,]) == "Chi2Closed_A")
which(names(out.gof$mcmc[[1]][1,]) == "Chi2repClosed_A")

pl <- range(c(as.matrix(out.gof[[1]][,11838][[1]]), as.matrix(out.gof[[1]][,11839][[1]])))
plot(as.matrix(out.gof[[1]][,11838][[1]]), as.matrix(out.gof[[1]][,11839][[1]]),
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Closed part of model (Chi-squared)", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(525, 720, paste('Bpv = ', round(mean(as.matrix(out.gof[[1]][,11839][[1]]) >
                                            as.matrix(out.gof[[1]][,11838][[1]])), 2)), cex = 2)




which(names(out.gof$mcmc[[1]][1,]) == "Chi2ratioOpen_A")
which(names(out.gof$mcmc[[1]][1,]) == "Chi2ratioOpen_B")

as.matrix(out.gof[[1]][,11837][[1]])
as.matrix(out.gof[[1]][,14267][[1]])

which(names(out.gof$mcmc[[1]][1,]) == "Chi2ratioClosed_A")
which(names(out.gof$mcmc[[1]][1,]) == "Chi2ratioClosed_B")

hist(as.matrix(out.gof[[1]][,11840][[1]]))
hist(as.matrix(out.gof[[1]][,14270][[1]]))





# Closed part of model: Freeman-Tukey

pl <- range(c(out.gof$sims.list$FTClosed_B, out.gof$sims.list$FTrepClosed_B))
plot(out.gof$sims.list$FTClosed_B, out.gof$sims.list$FTrepClosed_B,
     xlab = "FT observed data", ylab = "FT expected data",
     main = "Closed part of model (Freeman-Tukey)", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(240, 335, paste('Bpv = ', round(mean(out.gof$sims.list$FTrepClosed_B >
                                            out.gof$sims.list$FTClosed_B), 2)), cex = 2)
#par(op)

# Inspect observed and expected numbers of transitions for each interval
obs.trans <- out.gof$mean$tm_A
rep.trans <- out.gof$mean$tmrep_A
exp.trans <- out.gof$mean$Etm_A
dimnames(obs.trans) <- dimnames(rep.trans) <- dimnames(exp.trans) <-
  list(c("From Non-occ", "From Occ"), c("To Non-occ", "To Occ"))

for(t in 1:(data$nseason-1)){
  cat("\n*** Interval", t, " ***\n")
  cat("* Observed transitions*\n")
  print(obs.trans[,,t])
  
  cat("\n* Replicate data transitions*\n")
  print(rep.trans[,,t])
  
  cat("\n* Expected transitions*\n")
  print(exp.trans[,,t])
}