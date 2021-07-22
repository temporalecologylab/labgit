## Started 28 December 2019 ##
## By EM Wolkovich, with help from AK Ettinger #

## Simulations of winter to spring temperatures with fstar (GDD, beta in first-hitting time model) required increasing when chilling is low ##
## ... and! delayed leafout due to photoperiod cue ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

set.seed(113)

# Setting working directory (only needed to write out figures)
if(length(grep("ailene", getwd()))>0) { 
} else setwd("~/Documents/git/lab/labgit/projects/decsenspost")

# libraries
require(ggplot2)

#############################
## Chilling x forcing sims ##
#############################
# This code requires higher forcing (aka GDD, aka fstar) for leafout when chilling is low #
# It does not vary when GDD accumulates though, it always just starts accumulating on a certain day (daystostart) #

# Step 1: Set up years, days per year, temperatures, required GDD (fstar), required chill (cstar) and how much higher fstar is when cstar is not met
dayswinter <- 120
daysspring <- 90
wintertemp <- 1
springtemp <- 2
springtempinc <- 0.1
sigma <- 4
fstar <- 200 # GDD (beta in first-hitting-time model)
cstar <- 110 # threshold chilling
fstaradjforchill <- 3 # how much more GDD to you need based on diff from cstar at end of daystostart
yearz <- 50 # used to be 20, adjusted to 50 to compare to Jonathan's code 
sitez <- 45 # number of sites
simsnum <- 30 # how many sims to run (goes with below)
degreez <- seq(0, 7, length.out=simsnum) # warming -- applied evenly across `winter' and `spring'


## Step 2: Now, put together the seasonal temps, varying fstar (increases when chill is low) and calculate the sensitivities

df <- data.frame(degwarm=numeric(), rep=numeric(), chill=numeric(), fstar=numeric(), simplelm=numeric(),
    loglm=numeric(), perlm=numeric(), propryrschillmet=numeric(), meangddsum=numeric()) 

yearlytemp <- "alltemps"
for (i in degreez){
   for (j in 1:sitez){
       yearly_temp <- rep(0, yearz) # set up for sapply
       daily_temp <- sapply(yearly_temp, function(x) c(rnorm(dayswinter, wintertemp + i, sigma),
           (rnorm(daysspring, springtemp + i , sigma) + c(1:daysspring)*springtempinc))) # plot(daily_temp[,1] ~ c(1:length(daily_temp[,1])))
       chill <- daily_temp
       chill[(chill)<0] <- 0
       chill[(chill)>5] <- 0
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       gddreq <- c()
       leafout_date <- c()
       for (k in 1:ncol(chill)){
           chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
           gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[dayswinter:nrow(gdd),x])))
           if (chillsum[dayswinter,k]>cstar) {
           gddreq[k] <- fstar
           } else {
           gddreq[k] <- fstar + (cstar-chillsum[dayswinter,k])*fstaradjforchill
           }
           leafout_date[k] <- min(which(gddsum[,k] > gddreq[k])) # leafout_date unit of days *after* dayswinter
           # calculating a few means to include in dataframe
           meanchill <- mean(chillsum[dayswinter,]) 
           meanfstar <- mean(gddreq)
           chillmet<-length(which(chillsum[dayswinter,]>cstar))/yearz
           meangddsum<- mean(gddsum)
           }
           yearly_tempall <- colMeans(daily_temp)
           yearly_temppostwinter <- colMeans(daily_temp[dayswinter:nrow(daily_temp),])
           if(yearlytemp=="postwinter") {
               yearly_temp <- yearly_temppostwinter
               } else {
               yearly_temp <- yearly_tempall
               }
           per_leafout_date <- leafout_date/mean(leafout_date)
           per_yearly_temp <- yearly_temp/mean(yearly_temp)
           # plot(yearly_temp, leafout_date, pch=20)
           dfadd <- data.frame(degwarm=i, rep=j, chill=meanchill, fstar=meanfstar,     
               simplelm=coef(lm(leafout_date~yearly_temp))[2],
               loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
               perlm=coef(lm(per_leafout_date~per_yearly_temp))[2],
               propryrschillmet = chillmet,meangddsum= meangddsum)
           df <- rbind(df, dfadd)     
       }
   }



plot(simplelm~degwarm, data=df, pch=16, ylab="Sensitivity (days/C or log(days)/log(C)", xlab="Degree warming")
points(loglm~degwarm, data=df, col="dodgerblue")
points(perlm~degwarm, data=df, col="firebrick")

plot(propryrschillmet~degwarm, data=df, pch=16, ylab="Number of years (out of 30) when chilling is met", xlab="Degree warming")
plot(fstar~degwarm, data=df, pch=16, ylab="GDD", xlab="Degree warming")
plot(meangddsum~degwarm, data=df, pch=16, ylab="GDD", xlab="Degree warming")


################################
## Plotting of chilling sims ##
###############################

mean.sims <- aggregate(df[c("simplelm", "loglm", "perlm","propryrschillmet", "fstar","meangddsum")], df["degwarm"], FUN=mean)
sd.sims <- aggregate(df[c("simplelm", "loglm", "perlm","propryrschillmet", "fstar","meangddsum")], df["degwarm"], FUN=sd)

cexhere <- 0.95
cextext <- 0.8
pdf(file.path("figures/shiftingcuessims_2panels.pdf"), width = 6, height = 8)
par(mfrow=c(2,1),mar=c(5,5,2,5))
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-15, 5),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="", bty="l", mgp=c(1.5,.5,0), tck=-.01)
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
legend("bottomright", pch=c(19, 19), col=c("darkblue","salmon"), legend=c("Linear (untransformed)", "Non-linear (logged)"),
   cex=cextext, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-0.1, 1),mgp=c(1.5,.5,0), tck=-.01,xaxs="i",yaxs = "i",
     ylab="Proportion years when chilling is met",
     xlab=expression(paste("Warming (", degree, "C)")), bty="u",main="")
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$propryrschillmet[i]
  sdhere <- sd.sims$propryrschillmet[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkgray")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkgray")
}
par(new = TRUE)
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(200,300),yaxt="n", ylab="",xaxt="n", xlab="", bty="u",mgp=c(1.5,.5,0), tck=-.01)
axis(side = 4,mgp=c(1.5,.5,0), tck=-.01)
mtext(expression(paste("Thermal sum required for leafout (", degree, "C)"), sep=""),side=4, adj=.5, line=2)

for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$fstar[i]
  sdhere <- sd.sims$fstar[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkred")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkred")
}
text(mean.sims$degwarm[1] + 1.6, mean.sims$fstar[1] + 6,
     labels=expression(paste("Thermal sum (", degree, "C)"), sep=""), cex=cextext,
      col="darkred")
text(mean.sims$degwarm[1] + 0.8, mean.sims$fstar[1] + 80,
     labels=expression(paste("Chilling met"), sep=""), cex=cextext,
      col="darkgray")

dev.off()

##############################
## Daylength x forcing sims ##
##############################

# There are two oft-mentioned scenarios where photoperiod cues kick in with climate change
# (1) Plants reach threshold GDD too soon and thus wait for a certain photoperiod
# (2) Chilling is not met, so required GDD goes way up and plants never reach that, so they just leaf out at a certain photoperiod
# I model a simplified version of this below (effectively (1), but it's pretty similar to (2))
# Here I modeled a simple threshold photoperiod (I did this to make it more different than the chill scenarios), but ...
# could easily make it interactive similar to chilling code above.

# Step 1: Set up years, days per year, temperatures, required GDD (fstar), required photoperiod (pstar)

library(geosphere) # for daylengths

dayswinter <- 60 # sims are set up as though starting on Jan 1 (so 60 here means 'winter temps' end in early March)
daysspring <- 90
wintertemp <- 0
springtemp <- 4
springtempinc <- 0.1
photoper <- daylength(45, 1:(dayswinter+daysspring)) # latitude=45
sigma <- 4
fstar <- 200
pstar <- 12 # threshold to leafout
pstarday <- min(which(photoper > pstar))
yearz <- 50
sitez <- 45
simsnum <- 30
degreez <- round(seq(0, 7, length.out=simsnum), digits=1) # warming -- applied evenly across `winter' and `spring'

# Step 2: Run the simulations
dfphoto <- data.frame(degwarm=numeric(), rep=numeric(), meantemp=numeric(), leafoutdoy=numeric(), gddmetday=numeric(), 
    propyrsphoto=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric()) 

for (i in degreez){
   for (j in 1:sitez){
       yearly_temp <- rep(0, yearz) 
       daily_temp <- sapply(yearly_temp, function(x) c(rnorm(dayswinter, wintertemp + i, sigma),
           (rnorm(daysspring, springtemp + i , sigma)+ c(1:daysspring)*springtempinc)))
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       leafout_date <- c()
       gddmetday <- c()
       for (k in 1:ncol(daily_temp)){
           gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[dayswinter:nrow(gdd),x])))
           gddmetday[k] <- min(which(gddsum[,k] > fstar))
           if (photoper[(gddmetday[k] + dayswinter)] > pstar) { # check if the photoperiod threshold is met by the GDD-met day
            leafout_date[k] <- gddmetday[k]
           } else {
            leafout_date[k] <- pstarday-dayswinter # keep on same day scale as gdd
           }
           }
           yearly_temp <- colMeans(daily_temp)
           per_leafout_date <- leafout_date/mean(leafout_date)
           per_yearly_temp <- yearly_temp/mean(yearly_temp)
           driverbyyear <- NA
           driverbyyear[which(leafout_date==gddmetday)] <- "forcing"
           driverbyyear[which(leafout_date!=gddmetday)] <- "photo"
           driverbyyear[which(leafout_date==(pstarday-dayswinter))] <- "forcing/photo"
           photodriver <- length(which(driverbyyear!="forcing"))/yearz # here I pick only definitive forcing years!
           dfadd <- data.frame(degwarm=i, rep=j, meantemp=mean(yearly_temp),
               leafoutdoy=mean(leafout_date), gddmetday= mean(gddmetday),
               propyrsphoto=photodriver,
               simplelm=coef(lm(leafout_date~yearly_temp))[2],
               loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
               perlm=coef(lm(per_leafout_date~per_yearly_temp))[2])
           dfphoto <- rbind(dfphoto, dfadd)     
       }
   }

# plot by whether leafout was due to forcing and/or photoperiod (averaged over years)
# propyrsphoto for better resolution
dfphoto$leafoutdriver <- NA
dfphoto$leafoutdriver[which(dfphoto$propyrsphoto==0)] <- "all forcing"
dfphoto$leafoutdriver[which(dfphoto$propyrsphoto==1)] <- "all photo" # in current sims does not always occur
dfphoto$leafoutdriver[which(dfphoto$propyrsphoto>0 & dfphoto$propyrsphoto<1)] <- "photo/forcing"
dfphoto$degwarmtext <- paste("warming:", as.factor(dfphoto$degwarm), "C")

if(FALSE){ # extra plots, these are a little slow
ggplot(data=dfphoto, aes(x=meantemp, y=leafoutdoy, group=leafoutdriver, color=leafoutdriver)) +
   geom_point() +
   geom_smooth(method = "lm", linetype = 2, lwd=0.5, se = FALSE) +
   facet_wrap(.~degwarmtext, scales="free") +
        ylab("Mean day of year") +
    xlab(expression(paste("Mean daily temperature (", degree, "C)"), sep=""))

ggplot(data=dfphoto, aes(x=meantemp, y=leafoutdoy, color=propyrsphoto)) +
   geom_point() +
   geom_smooth(method = "lm", linetype = 2, lwd=0.5, se = FALSE) +
   facet_wrap(.~degwarmtext, scales="free") +
        ylab("Mean day of year") +
    xlab(expression(paste("Mean daily temperature (", degree, "C)"), sep=""))
}

# plot the means
mean.simsphoto <- aggregate(dfphoto[c("simplelm", "loglm", "perlm","leafoutdoy", "gddmetday", "propyrsphoto")], dfphoto["degwarm"], FUN=mean)
sd.simsphoto <- aggregate(dfphoto[c("simplelm", "loglm", "perlm","leafoutdoy", "gddmetday", "propyrsphoto")], dfphoto["degwarm"], FUN=sd)

cexhere <- 0.95
cextext <- 0.8

pdf(file.path("figures/shiftingcuessims_photo2panel.pdf"), width = 6, height = 8)
par(mfrow=c(2,1), mar=c(5,5,2,5))
plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez) + 0.5)), ylim=c(-8, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="", bty="l", mgp=c(1.5,.5,0), tck=-.01)
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$simplelm[i]
  sdhere <- sd.simsphoto$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$loglm[i]
  sdhere <- sd.simsphoto$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue","salmon"), legend=c("Linear (untransformed)", "Non-linear (logged)"),
   cex=cextext, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez)+ 0.5)), ylim=c(-0.1, (max(mean.simsphoto$propyrsphoto)+0.1)),
     mgp=c(1.5,.5,0), tck=-.01,xaxs="i",yaxs = "i",
     ylab="Proportion years photoperiod drives leafout",
     xlab=expression(paste("Warming (", degree, "C)")), bty="l",main="")
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$propyrsphoto[i]
  sdhere <- sd.simsphoto$propyrsphoto[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkgray")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkgray")
}
dev.off()



# Same as above, but for just one site ...
dfphoto_onesite <- data.frame(degwarm=numeric(), rep=numeric(), yearly_temp=numeric(), leafoutdoy=numeric(),
    gddmetday=numeric(), driver=numeric()) 

for (i in degreez){
       yearly_temp <- rep(0, yearz) 
       daily_temp <- sapply(yearly_temp, function(x) c(rnorm(dayswinter, wintertemp + i, sigma),
           (rnorm(daysspring, springtemp + i , sigma)+ c(1:daysspring)*springtempinc)))
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       leafout_date <- c()
       gddmetday <- c()
       for (k in 1:ncol(daily_temp)){
           gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[dayswinter:nrow(gdd),x])))
           gddmetday[k] <- min(which(gddsum[,k] > fstar))
           if (photoper[(gddmetday[k] + dayswinter)] > pstar) { # check if the photoperiod threshold is met by the GDD-met day
            leafout_date[k] <- gddmetday[k]
           } else {
            leafout_date[k] <- pstarday-dayswinter # keep on same day scale as gdd
           }
           }
           yearly_temp <- colMeans(daily_temp)
           per_leafout_date <- leafout_date/mean(leafout_date)
           per_yearly_temp <- yearly_temp/mean(yearly_temp)
           driverbyyear <- NA
           driverbyyear[which(leafout_date==gddmetday)] <- "forcing"
           driverbyyear[which(leafout_date!=gddmetday)] <- "photo"
           driverbyyear[which(gddmetday==(pstarday-dayswinter))] <- "forcing/photo"
           dfadd <- data.frame(degwarm=rep(i, yearz), yearly_temp=yearly_temp,
               leafoutdoy=leafout_date, gddmetday=gddmetday,
               driver=driverbyyear)
           dfphoto_onesite <- rbind(dfphoto_onesite, dfadd)     
       }


dfphoto_onesite$degwarmtext <- paste("warming:", as.factor(dfphoto_onesite$degwarm), "C")

ggplot(data=dfphoto_onesite, aes(x=yearly_temp, y=leafoutdoy, group=driver, color=driver)) +
   geom_point() +
   geom_smooth(method = "lm", linetype = 2, lwd=0.5, se = FALSE) +
   facet_wrap(.~degwarmtext, scales="free") +
        ylab("Day of year") +
    xlab(expression(paste("Mean daily temperature (", degree, "C)"), sep=""))


###################################
## Plots chilling and daylength  ##
##          sims together        ##
###################################


pdf(file.path("figures/shiftingcuessims_4panels.pdf"), width = 12, height = 8)
par(mfrow=c(2,2), mar=c(5,5,2,5))
plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez) + 0.5)), ylim=c(-10, 5),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="", bty="l", mgp=c(1.5,.5,0), tck=-.01)
mtext("(a)", side = 3, adj=0.05)
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
legend("bottomright", pch=c(19, 19), col=c("darkblue","salmon"), legend=c("Linear (untransformed)", "Non-linear (logged)"),
   cex=cexhere, bty="n")


plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez) + 0.5)), ylim=c(-6, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="", bty="l", mgp=c(1.5,.5,0), tck=-.01)
mtext("(b)", side = 3, adj=0.05)
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$simplelm[i]
  sdhere <- sd.simsphoto$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$loglm[i]
  sdhere <- sd.simsphoto$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue","salmon"), legend=c("Linear (untransformed)", "Non-linear (logged)"),
   cex=cexhere, bty="n")

plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez) + 0.5)), ylim=c(-0.1, 1),mgp=c(1.5,.5,0), tck=-.01,xaxs="i",yaxs = "i",
     ylab="Proportion years when chilling is met",
     xlab=expression(paste("Warming (", degree, "C)")), bty="u",main="")
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$propryrschillmet[i]
  sdhere <- sd.sims$propryrschillmet[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkgray")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkgray")
}
par(new = TRUE)
plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez) + 0.5)), ylim=c(200,300),yaxt="n", ylab="",xaxt="n", xlab="", bty="u",mgp=c(1.5,.5,0), tck=-.01)
axis(side = 4,mgp=c(1.5,.5,0), tck=-.01)
mtext(expression(paste("Thermal sum required for leafout (", degree, "C)"), sep=""), side=4, adj=.5, line=2, cex=cexhere)
mtext("(c)", side = 3, adj=0.05)
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$fstar[i]
  sdhere <- sd.sims$fstar[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkred")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkred")
}
text(mean.sims$degwarm[1] + 1.6, mean.sims$fstar[1] + 6,
     labels=expression(paste("Thermal sum (", degree, "C)"), sep=""), cex=cextext,
      col="darkred")
text(mean.sims$degwarm[1] + 0.8, mean.sims$fstar[1] + 80,
     labels=expression(paste("Chilling met"), sep=""), cex=cextext,
     col="darkgray")

plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez)+ 0.5)), ylim=c(-0.1, (max(mean.simsphoto$propyrsphoto)+0.1)),
     mgp=c(1.5,.5,0), tck=-.01,xaxs="i",yaxs = "i",
     ylab="Proportion years photoperiod drives leafout",
     xlab=expression(paste("Warming (", degree, "C)")), bty="l",main="")
mtext("(d)", side = 3, adj=0.05)
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$propyrsphoto[i]
  sdhere <- sd.simsphoto$propyrsphoto[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkgray")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkgray")
}
dev.off()
