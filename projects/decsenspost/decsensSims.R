## Started 4 January 2020 ##
## By EM Wolkovich, help from J Auerbach ##

## Simulation of the declining sensitivities ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

set.seed(113)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("ailene", getwd()))>0) { 
} else setwd("~/Documents/git/lab/labgit/projects/decsenspost")

##########################
# The below sets up data #
# for Fig 1 in main text #
##########################

# Make some data ... note that this runs 100 times for 45 sites, via a loop (slow)
# Set simsnum.maintext=500 for the exact figure in the main text   

# Step 1: Set up years, days per year, temperatures, required GDD (fstar)
daysperyr <- 60
yearz <- 30
sitez <- 45 # reps
simsnum.maintext <- 100
degreez.maintext <- seq(0, 2, length.out=simsnum.maintext)
degreez.forsupp <- c(0, 0.5, 1, 1.5, 2, 2.5, 4, 7)
sigma <- 4
basetemp <- 4 # alpha_0
dailytempchange <- 0.1 # alpha_1
fstar <- 150 # beta

# Step 2: Build the data and calculate sensitivities (note that alpha_1 -- spring temperature increase is set to 0)
df <- data.frame(degwarm=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric(),
    simplelm.trunc=numeric(), loglm.trunc=numeric())

for (i in degreez.maintext){
   for (j in 1:sitez){
       yearly_expected_temp <- rep(basetemp, yearz)
       daily_temp <- sapply(yearly_expected_temp, function(x) rnorm(daysperyr, basetemp + i, sigma))
       daily_temp <- daily_temp + c(1:daysperyr)*dailytempchange # add daily temp increase
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > fstar)))
       yearly_temp <- colMeans(daily_temp)
       yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x]))
       per_leafout_date <- leafout_date/mean(leafout_date)
       per_yearly_temp <- yearly_temp/mean(yearly_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfadd <- data.frame(degwarm=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
           perlm=coef(lm(per_leafout_date~per_yearly_temp))[2],
           simplelm.trunc=coef(lm(leafout_date~yearly_temp_trunc))[2],
           loglm.trunc=coef(lm(log(leafout_date)~log(yearly_temp_trunc)))[2])
       df <- rbind(df, dfadd)
    }
}


dfsupp <- data.frame(degwarm=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric(),
    simplelm.trunc=numeric(), loglm.trunc=numeric())

for (i in degreez.forsupp){
   for (j in 1:sitez){
       yearly_expected_temp <- rep(basetemp, yearz)
       daily_temp <- sapply(yearly_expected_temp, function(x) rnorm(daysperyr, basetemp + i, sigma))
       daily_temp <- daily_temp + c(1:daysperyr)*dailytempchange  # add daily temp increase
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > fstar)))
       yearly_temp <- colMeans(daily_temp)
       yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x]))
       per_leafout_date <- leafout_date/mean(leafout_date)
       per_yearly_temp <- yearly_temp/mean(yearly_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfsuppadd <- data.frame(degwarm=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
           perlm=coef(lm(per_leafout_date~per_yearly_temp))[2],
           simplelm.trunc=coef(lm(leafout_date~yearly_temp_trunc))[2],
           loglm.trunc=coef(lm(log(leafout_date)~log(yearly_temp_trunc)))[2])
       dfsupp <- rbind(dfsupp, dfsuppadd)
    }
}

plot(simplelm.trunc~degwarm, data=dfsupp, ylab="Sensitivity (days/C or log(days)/log(C)", xlab="Degree warming")
points(simplelm~degwarm, data=dfsupp, pch=16, col="gray")
points(loglm.trunc~degwarm, data=dfsupp, col="dodgerblue")
points(loglm~degwarm, data=dfsupp, col="dodgerblue", pch=16)
plot(perlm~degwarm, data=dfsupp, col="firebrick")

plot(abs(simplelm)~degwarm, data=dfsupp, col="lightgrey",
    ylab="Abs(Sensitivity (days/C or log(days)/log(C))", xlab="Degree warming")
dfsupp$degwarmJitter <- dfsupp$degwarm + 0.05
points(abs(loglm)~degwarmJitter, data=dfsupp, col="dodgerblue", cex=0.8)


##############
## Plotting ##
##############

# Summarize the sims
mean.sims <- aggregate(dfsupp[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], dfsupp["degwarm"], FUN=mean)
sd.sims <- aggregate(dfsupp[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], dfsupp["degwarm"], FUN=sd)

cexhere <- 0.95
pdf(file.path("figures/basicsims.pdf"), width = 6, height = 8)
par(mfrow=c(2,1))
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-6, -0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C to leafout)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="",
     bty="l", mgp=c(1.5,.5,0), tck=-.01)
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$simplelm.trunc[i]
  sdhere <- sd.sims$simplelm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$loglm.trunc[i]
  sdhere <- sd.sims$loglm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("salmon", "darkblue"), legend=c("Using logged variables", "Simple linear regression"),
   cex=1, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-6, -0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C over window)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="",
     bty="l", mgp=c(1.5,.5,0), tck=-.01)
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$simplelm[i]
  sdhere <- sd.sims$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$loglm[i]
  sdhere <- sd.sims$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("salmon", "darkblue"), legend=c("Using logged variables", "Simple linear regression"),
   cex=1, bty="n")
dev.off()

########################
## Plotting, plus PEP ##
########################

# Summarize the main text sims
mean.sims <- aggregate(df[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], df["degwarm"], FUN=mean)
sd.sims <- aggregate(df[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], df["degwarm"], FUN=sd)

# Get data ...
dfpep <- read.csv("pep_analyses/output/bpenestimates_withlog_1950_2000.csv", header=TRUE)

# Using 60 d as main one, so change names
names(dfpep)[names(dfpep)=="mat60slope"] <- "matslope"
names(dfpep)[names(dfpep)=="mat60slopelog"] <- "matslopelog"
names(dfpep)[names(dfpep)=="meanmat60"] <- "meanmat"
    
# Get means and SD
mean.betpen <- aggregate(dfpep[c("matslope", "matslopelog", "meanmat")], dfpep["cc"], FUN=mean)
sd.betpen <- aggregate(dfpep[c("matslope", "matslopelog", "meanmat")], dfpep["cc"],  FUN=sd)

tempdiff <- mean.betpen$meanmat[which(mean.betpen$cc=="2000-2010")]-
    mean.betpen$meanmat[which(mean.betpen$cc=="1950-1960")]
tempdiffplot <- c(0, tempdiff)

library(grDevices)
colz <- c("blue4", "violetred4", "blue1", "violetred1")
colzalpha <- adjustcolor(colz, alpha.f = 0.7)

cexhere <- 0.75
cexhere <- 1.2
cextext <- 0.75
jitterpep <- -0.04

# Set up for shading the sims
# THIS is the figure currently in the main text ...
simsrange <- seq(min(mean.sims$degwarm), max(mean.sims$degwarm), length.out=simsnum.maintext)
sdupp <- mean.sims$simplelm + sd.sims$simplelm
sdlow <- mean.sims$simplelm - sd.sims$simplelm

sdupplog <- mean.sims$loglm + sd.sims$loglm
sdlowlog <- mean.sims$loglm - sd.sims$loglm


pdf(file.path("figures/basicsimsandpepalt1.pdf"), width = 9, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
par(mfrow=c(1,2))
plot(x=NULL,y=NULL, xlim=c(-0.25, 2.25), ylim=c(-6.6, -0.5), 
     ylab=expression(paste("Estimated sensitivity"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="Linear (untransformed)", font.main = 1, cex.main = 0.9, 
     cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
legend("bottomright", pch=c(19, 19), col=c(colzalpha[2], colzalpha[3]), legend=c("simulations", "observations"),
   cex=1, bty="n")
# axis(2,seq(-6,0,1), las=2)
tempsteps <- simsnum.maintext
tempdiffplot <- c(0,1)
polygon(c(rev(simsrange), simsrange), c(rev(sdupp), sdlow), col = colzalpha[2], border = NA)
for(i in 1:tempsteps){
  lines(pos.x, pos.y, cex=cexhere, col=colzalpha[2])
}
for(i in 1:length(unique(mean.betpen$cc))){ # i=2
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslope[i]
  sdhere <- sd.betpen$matslope[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[3])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[3])
  text(pos.x + 0.35, pos.y-1.1, labels=unique(mean.betpen$cc)[i], 
       cex=cextext, col=colzalpha[3])
}
plot(x=NULL,y=NULL, xlim=c(-0.25, 2.25), ylim=c(-1.5, 0.2), 
     ylab=expression(paste("Estimated sensitivity"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="Non-linear (logged)", font.main = 1, cex.main = 0.9,
     cex.lab=1.2,
     bty="l", mgp=c(1.5,.25,0), tck=-.01)
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslopelog[i]
  sdhere <- sd.betpen$matslopelog[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[3])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[3])
  text(pos.x + 0.37, pos.y, labels=unique(mean.betpen$cc)[i], cex=cextext, col=colzalpha[3])
}
polygon(c(rev(simsrange), simsrange), c(rev(sdupplog), sdlowlog), col = colzalpha[2], border = NA)
for(i in 1:tempsteps){
  lines(pos.x, pos.y, cex=cexhere,col=colzalpha[2])
}
dev.off()

############################
# Sims and plotting Fig S1 #
# Once in decsensAuerbach.R #
############################


# First, build up some simulated climate and leafout (after 200 GDD) data
yearly_expected_temp <- c(rep(0,20), rep(5,20), rep(10,20), rep(20, 20))  # create holder for three sets of temps, 20 per set (note, these are bigger than we have done -- we have done 1 to 7 degrees, this does 0, 5, 10, and 20 degrees)

daily_temp <- sapply(yearly_expected_temp, function(x) rnorm(30, 10 + x, 4)) # now make some simple daily temps
leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > 200))) # set leafout date as whenever 200 GDD is reached
yearly_temp <- colMeans(daily_temp) # estimate the mean temp of each simulated dataset
yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x])) # estimate the mean temp of each simulated dataset only until leafout 

plot(yearly_temp, leafout_date, pch=20)
points(yearly_temp_trunc, leafout_date, pch=20, col = "red")
plot(yearly_temp_trunc, leafout_date, pch=20, col = "red")

# Figure S1 currently
cexhere <- 0.5

plot(log(yearly_temp_trunc), log(leafout_date), pch=20, col = "dodgerblue") 
pdf(file.path("figures/simslogging.pdf"), width = 9, height = 5)
par(mfrow=c(2,3))
plot(yearly_temp_trunc, leafout_date, pch=20, xlab="Simulated spring temperature to leafout",
     ylab="Leafout date", main="", cex=cexhere, bty="l", mgp=c(1.5,.5,0), tck=-.01)
plot(yearly_temp_trunc, log(leafout_date), pch=20, xlab="Simulated spring temperature to leafout",
     ylab="log(Leafout date)", main="", cex=cexhere, bty="l", mgp=c(1.5,.5,0), tck=-.01)
plot(log(yearly_temp_trunc), log(leafout_date), pch=20, xlab="log(Simulated spring temperature to leafout)",
     ylab="log(Leafout date)", main="", cex=cexhere, bty="l", mgp=c(1.5,.5,0), tck=-.01)
plot(yearly_temp, leafout_date, pch=20, xlab="Simulated spring temperature",
    ylab="Leafout date", main="", cex=cexhere, bty="l", mgp=c(1.5,.5,0), tck=-.01)
plot(yearly_temp, log(leafout_date), pch=20, xlab="Simulated spring temperature",
    ylab="log(Leafout date)", main="", cex=cexhere, bty="l", mgp=c(1.5,.5,0), tck=-.01)
plot(log(yearly_temp), log(leafout_date), pch=20, xlab="log(Simulated spring temperature)",
     ylab="log(Leafout date)", main="", cex=cexhere, bty="l", mgp=c(1.5,.5,0), tck=-.01)
dev.off()


