## Started 7 September 2020 ##
## By EM Wolkovich, based off decsensSims.R ##

## What happens with no warming where only fstar (aka beta) varies? ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

set.seed(113)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("ailene", getwd()))>0) { 
} else setwd("~/Documents/git/lab/labgit/projects/decsenspost")

##########################
# The below sets up data #
##########################

# Make some data ... note that this runs simsnum times for 45 sites, via a loop

# Step 1: Set up years, days per year, temperatures, required GDD (fstar)
daysperyr <- 60
yearz <- 30
sitez <- 45 # reps
simsnum <- 40
fstarsims <- seq(100, 300, length.out=simsnum)
sigma <- 4
basetemp <- 4 # alpha_0
dailytempchange <- 0.1 # alpha_1

# Step 2: Build the data and calculate sensitivities 
df <- data.frame(fstar=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric(),
    simplelm.trunc=numeric(), loglm.trunc=numeric(), varx=numeric(), covarxy=numeric(), covarlogxy=numeric(),
    varx.trunc=numeric(), covarxy.trunc=numeric(), covarlogxy.trunc=numeric())

for (i in fstarsims){
   for (j in 1:sitez){
       daily_temp <- sapply(rep(NA, yearz), function(x) rnorm(daysperyr, basetemp, sigma))
       daily_temp <- daily_temp + c(1:daysperyr)*dailytempchange 
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > i)))
       yearly_temp <- colMeans(daily_temp)
       yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x]))
       per_leafout_date <- leafout_date/mean(leafout_date)
       per_yearly_temp <- yearly_temp/mean(yearly_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfadd <- data.frame(fstar=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
           perlm=coef(lm(per_leafout_date~per_yearly_temp))[2],
           simplelm.trunc=coef(lm(leafout_date~yearly_temp_trunc))[2],
           loglm.trunc=coef(lm(log(leafout_date)~log(yearly_temp_trunc)))[2],
           varx=var(yearly_temp), covarxy=cov(leafout_date, yearly_temp),
           covarlogxy=cov(log(leafout_date), log(yearly_temp)),
           varx.trunc=var(yearly_temp_trunc), covarxy.trunc=cov(leafout_date, yearly_temp_trunc),
           covarlogxy.trunc=cov(log(leafout_date), log(yearly_temp_trunc))
           )
       df <- rbind(df, dfadd)
    }
}

dfsm <- data.frame(fstar=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric())

for (i in fstarsims){
   for (j in 1:sitez){
       daily_temp <- sapply(rep(NA, yearz), function(x) rnorm(daysperyr, basetemp, sigma))
       daily_temp <- daily_temp + c(1:daysperyr)*dailytempchange 
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > i)))
       yearly_temp <- colMeans(daily_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfadd <- data.frame(fstar=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2])
       dfsm <- rbind(dfsm, dfadd)
    }
}

mean.sims.sm <- aggregate(dfsm[c("simplelm", "loglm")], dfsm["fstar"], FUN=mean)


##############
## Plotting ##
##############
pdf(file.path("figures/basicsims_fstaronly_varcov.pdf"), width = 10, height = 5)
par(mfrow=c(2,4))
plot(simplelm.trunc~fstar, data=df, xlab="thermal sum", ylab="lm sensitivity (temp until leafout)")
plot(varx.trunc~fstar, data=df, xlab="thermal sum", ylab="var(temp until leafout)")
plot(covarxy.trunc~fstar, data=df, xlab="thermal sum", ylab="covar(leafout day, temp until leafout)")
plot(covarlogxy.trunc~fstar, data=df, xlab="thermal sum", ylab="covar(log(leafout day), log(temp until leafout))")
plot(simplelm~fstar, data=df, xlab="thermal sum", ylab="lm sensitivity (temp over window)")
plot(varx~fstar, data=df, xlab="thermal sum", ylab="var(temp over window)")
plot(covarxy~fstar, data=df, xlab="thermal sum", ylab="covar(leafout day, temp over window)")
plot(covarlogxy~fstar, data=df, xlab="thermal sum", ylab="covar(log(leafout day), (log(temp over window))")
dev.off()

# Summarize the fstar sims
mean.sims <- aggregate(df[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], df["fstar"], FUN=mean)
sd.sims <- aggregate(df[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], df["fstar"], FUN=sd)

colz <- c("blue4", "salmon")
colzalpha <- adjustcolor(colz, alpha.f = 0.7)
cexhere <- 0.75
cexhere <- 1.2
cextext <- 0.75
jitterpep <- -0.04
pdf(file.path("figures/basicsims_fstaronly.pdf"), width = 7.5, height = 11)
par(xpd=FALSE)
par(mfrow=c(2,1))
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(90, 310), ylim=c(-8, 0.1), yaxt="n",
     ylab=expression(paste("Estimated sensitivity"), sep=""),
     xlab=expression(paste("Thermal sum required (", degree, "C)")), main="", cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.5,0), tck=-.01)
axis(2,seq(-8,0,1),las=2)
tempsteps <- simsnum
tempdiffplot <- c(0,1)
for(i in 1:tempsteps){
  pos.x <- mean.sims$fstar[i]
  pos.y <- mean.sims$simplelm.trunc[i]
  sdhere <- sd.sims$simplelm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[1])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[1])
}
for(i in 1:tempsteps){
  pos.x <- mean.sims$fstar[i]
  pos.y <- mean.sims$loglm.trunc[i]
  sdhere <- sd.sims$loglm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[2])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[2])
}
legend("topright", pch=c(19, 19), col=c(colzalpha[2], colzalpha[1]), legend=c("Using logged variables", "Simple linear regression"),
   cex=1, bty="n")
plot(x=NULL,y=NULL, xlim=c(90, 310), ylim=c(-6.6, 0.1), yaxt="n",
     ylab=expression(paste("Estimated sensitivity"), sep=""),
     xlab=expression(paste("Thermal sum required (", degree, "C)")), main="", cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.5, 0), tck=-.01)
axis(2,seq(-6,0,1),las=2)
tempsteps <- simsnum
tempdiffplot <- c(0,1)
for(i in 1:tempsteps){
  pos.x <- mean.sims$fstar[i]
  pos.y <- mean.sims$simplelm[i]
  sdhere <- sd.sims$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[1])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[1])
}
for(i in 1:tempsteps){
  pos.x <- mean.sims$fstar[i]
  pos.y <- mean.sims$loglm[i]
  sdhere <- sd.sims$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[2])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[2])
}
dev.off()
