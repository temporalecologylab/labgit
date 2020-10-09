## Started 18 January 2020 ##
## By Lizzie ##

## PEP data ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Setting working directory. 
setwd("~/Documents/git/projects/treegarden/decsens/analyses")
#setwd("~/Documents/git/decsens/analyses")

# Get data ...
df.20yr <- read.csv("pep_analyses/output/bpenestimates_withlog_1950to2010.csv", header=TRUE) 
df.10yr <- read.csv("pep_analyses/output/bpenestimates_withlog_1950_2000.csv", header=TRUE) # 10-year estimates (1950-1960 and 2000-2010) 
df.old <- read.csv("pep_analyses/output/zarchive/bpenestimates_withlog.csv", header=TRUE)
# dfswa <- read.csv("pep_analyses/output/swaestimates_withlog.csv", header=TRUE)
fs.20yr <- read.csv("pep_analyses/output/fsestimates_withlog_1950to2010.csv", header=TRUE) # older version: fsylestimates_withlog_1950-2010
fs.10yr <- read.csv("pep_analyses/output/fsestimates_withlog_1950_2000.csv", header=TRUE) # 10-year estimates (1950-1960 and 2000-2010)
qr.20yr <- read.csv("pep_analyses/output/qrestimates_withlog_1950to2010.csv", header=TRUE) # older version: fsylestimates_withlog_1950-2010
qr.10yr <- read.csv("pep_analyses/output/qrestimates_withlog_1950_2000.csv", header=TRUE) # 10-year estimates (1950-1960 and 2000-2010)


#################
## Format data ##
#################

length(unique(df.20yr$siteslist))
length(unique(fs.20yr$siteslist))
length(unique(df.10yr$siteslist))
length(unique(fs.10yr$siteslist))


mean.betpen.20yr <- aggregate(df.20yr[c("mat30slope", "mat30slopelog", "meanmat30", "varmat30",
    "mat45slope", "mat45slopelog", "meanmat45", "varmat45","mat60slope", "mat60slopelog", "meanmat60", 
    "varmat60", "varlo", "meangdd", "meanmatlo", "mat30slopeconfint11", "mat30slopeconfint89",
    "mat30slopelogconfint11", "mat30slopelogconfint89", "mat45slopeconfint11", "mat45slopeconfint89",
    "mat45slopelogconfint11", "mat45slopelogconfint89", "mat60slopeconfint11", "mat60slopeconfint89",
    "mat60slopelogconfint11", "mat60slopelogconfint89")], df.20yr["cc"], FUN=mean)

mean.betpen.10yr <- aggregate(df.10yr[c("mat30slope", "mat30slopelog", "meanmat30", "varmat30",
     "mat45slope", "mat45slopelog", "meanmat45", "varmat45","mat60slope", "mat60slopelog", "meanmat60", 
     "varmat60", "varlo", "meangdd", "meanmatlo", "mat30slopeconfint11", "mat30slopeconfint89",
     "mat30slopelogconfint11", "mat30slopelogconfint89",  "mat45slopeconfint11", "mat45slopeconfint89",
     "mat45slopelogconfint11", "mat45slopelogconfint89", "mat60slopeconfint11", "mat60slopeconfint89",
     "mat60slopelogconfint11", "mat60slopelogconfint89")], df.10yr["cc"], FUN=mean)

tempdiff1.20yr <- mean.betpen.20yr$meanmat60[which(mean.betpen.20yr$cc=="1970-1990")]-
    mean.betpen.20yr$meanmat60[which(mean.betpen.20yr$cc=="1950-1970")]
tempdiff2.20yr <- mean.betpen.20yr$meanmat60[which(mean.betpen.20yr$cc=="1990-2010")]-
    mean.betpen.20yr$meanmat60[which(mean.betpen.20yr$cc=="1950-1970")]
tempdiffplot.20yr <- c(0, tempdiff1.20yr, tempdiff2.20yr)

tempdiff1.10yr <- mean.betpen.10yr$meanmat60[which(mean.betpen.10yr$cc=="2000-2010")]-
    mean.betpen.10yr$meanmat60[which(mean.betpen.10yr$cc=="1950-1960")]
tempdiffplot.10yr <- c(0, tempdiff1.10yr)

## For Fagus sylvatica
mean.fs.20yr <- aggregate(fs.20yr[c("mat30slope", "mat30slopelog", "meanmat30", "varmat30",
     "mat45slope", "mat45slopelog", "meanmat45", "varmat45","mat60slope", "mat60slopelog", "meanmat60", 
     "varmat60", "varlo", "meangdd", "meanmatlo", "mat30slopeconfint11", "mat30slopeconfint89",
     "mat30slopelogconfint11", "mat30slopelogconfint89",  "mat45slopeconfint11", "mat45slopeconfint89",
     "mat45slopelogconfint11", "mat45slopelogconfint89", "mat60slopeconfint11", "mat60slopeconfint89",
     "mat60slopelogconfint11", "mat60slopelogconfint89")],
     fs.20yr["cc"], FUN=mean)

mean.fs.10yr <- aggregate(fs.10yr[c("mat30slope", "mat30slopelog", "meanmat30", "varmat30",
     "mat45slope", "mat45slopelog", "meanmat45", "varmat45","mat60slope", "mat60slopelog", "meanmat60", 
     "varmat60", "varlo", "meangdd", "meanmatlo", "mat30slopeconfint11", "mat30slopeconfint89",
     "mat30slopelogconfint11", "mat30slopelogconfint89",  "mat45slopeconfint11", "mat45slopeconfint89",
     "mat45slopelogconfint11", "mat45slopelogconfint89", "mat60slopeconfint11", "mat60slopeconfint89",
     "mat60slopelogconfint11", "mat60slopelogconfint89")],
      fs.10yr["cc"], FUN=mean)

tempdiff1fs.20yr <- mean.fs.20yr$meanmat60[which(mean.fs.20yr$cc=="1970-1990")]-
    mean.fs.20yr$meanmat60[which(mean.fs.20yr$cc=="1950-1970")]
tempdiff2fs.20yr <- mean.fs.20yr$meanmat60[which(mean.fs.20yr$cc=="1990-2010")]-
    mean.fs.20yr$meanmat60[which(mean.fs.20yr$cc=="1950-1970")]

tempdiffplotfs.20yr <- c(0, tempdiff1fs.20yr, tempdiff2fs.20yr)

tempdiff1fs.10yr <- mean.fs.10yr$meanmat60[which(mean.fs.10yr$cc=="2000-2010")]-
    mean.fs.10yr$meanmat60[which(mean.fs.10yr$cc=="1950-1960")]
tempdiffplotfs.10yr <- c(0, tempdiff1fs.10yr)

## For Quercus robur - prep for tables only
mean.qr.20yr <- aggregate(qr.20yr[c("mat30slope", "mat30slopelog", "meanmat30", "varmat30",
    "mat45slope", "mat45slopelog", "meanmat45", "varmat45","mat60slope", "mat60slopelog", "meanmat60", 
    "varmat60", "varlo", "meangdd", "meanmatlo", "mat30slopeconfint11", "mat30slopeconfint89",
    "mat30slopelogconfint11", "mat30slopelogconfint89",  "mat45slopeconfint11", "mat45slopeconfint89",
    "mat45slopelogconfint11", "mat45slopelogconfint89", "mat60slopeconfint11", "mat60slopeconfint89",
    "mat60slopelogconfint11", "mat60slopelogconfint89")],
    qr.20yr["cc"], FUN=mean)

mean.qr.10yr <- aggregate(qr.10yr[c("mat30slope", "mat30slopelog", "meanmat30", "varmat30",
    "mat45slope", "mat45slopelog", "meanmat45", "varmat45","mat60slope", "mat60slopelog", "meanmat60", 
    "varmat60", "varlo", "meangdd", "meanmatlo", "mat30slopeconfint11", "mat30slopeconfint89",
    "mat30slopelogconfint11", "mat30slopelogconfint89",  "mat45slopeconfint11", "mat45slopeconfint89",
    "mat45slopelogconfint11", "mat45slopelogconfint89", "mat60slopeconfint11", "mat60slopeconfint89",
    "mat60slopelogconfint11", "mat60slopelogconfint89")],
    qr.10yr["cc"], FUN=mean)

##############
## Plotting ##
##############

cexhere <- 0.95
cextext <- 0.75
pdf(file.path("figures/basicpep1950to2000.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.20yr$cc))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slope[i]
  ciherelo <- mean.betpen.20yr$mat60slopeconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext, col="darkblue")
  }
for(i in 1:length(unique(mean.betpen.20yr$cc))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slopelog[i]
  ciherelo <- mean.betpen.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Using raw x, y", "Using logged x, y"),
   cex=1, bty="n")
dev.off()


cexhere <- 0.95
pdf(file.path("figures/basicpep1950to2000fs.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.fs.20yr$cc))){
  pos.x <- tempdiffplotfs.20yr[i]
  pos.y <- mean.fs.20yr$mat60slope[i]
  ciherelo <- mean.fs.20yr$mat60slopeconfint11[i]
  cihereup <- mean.fs.20yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext, col="darkblue")
  }
for(i in 1:length(unique(mean.fs.20yr$cc))){
  pos.x <- tempdiffplotfs.20yr[i]
  pos.y <- mean.fs.20yr$mat60slopelog[i]
  ciherelo <- mean.fs.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.fs.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Using raw x, y", "Using logged x, y"),
   cex=1, bty="n")
dev.off()


##
## 2 panel with FagSyl 10 yr -- showing raw on top panels, and just logged on bottom panels
fagsyljitter <- 0.04
cexhere <- 0.5
cexhereleg <- 0.7
pdf(file.path("figures/basicpep1950to20102spp2panel.pdf"), width = 5, height = 7)
par(xpd=FALSE)
par(mfrow=c(2,1))
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")))
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.10yr$cc))){
  pos.x <- tempdiffplot.10yr[i]
  pos.y <- mean.betpen.10yr$mat60slope[i]
  ciherelo <- mean.betpen.10yr$mat60slopeconfint11[i]
  cihereup <- mean.betpen.10yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.fs.10yr$cc))){
  pos.x <- tempdiffplotfs.10yr[i]
  pos.y <- mean.fs.10yr$mat60slope[i]
  ciherelo <- mean.fs.10yr$mat60slopeconfint11[i]
  cihereup <- mean.fs.10yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x+fagsyljitter, 2), y=c(ciherelo, cihereup), col="dodgerblue")
  points(pos.x+fagsyljitter, pos.y, cex=cexhere, pch=19, col="dodgerblue")
  text(pos.x + 0.03, pos.y, labels=unique(mean.fs.10yr$cc)[i], cex=cextext, col="black")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "dodgerblue"),
   legend=c("Betula pendula", "Fagus sylvatica"), cex=cexhereleg, bty="n")

# Log-log 
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-0.9, 0.1),
     ylab=expression(paste("Estimated sensitivity (log(days)/log(", degree, "C))"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")))
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.10yr$cc))){
  pos.x <- tempdiffplot.10yr[i]
  pos.y <- mean.betpen.10yr$mat60slopelog[i]
  ciherelo <- mean.betpen.10yr$mat60slopelogconfint11[i]
  cihereup <- mean.betpen.10yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
for(i in 1:length(unique(mean.fs.10yr$cc))){
  pos.x <- tempdiffplotfs.10yr[i]
  pos.y <- mean.fs.10yr$mat60slopelog[i]
  ciherelo <- mean.fs.10yr$mat60slopelogconfint11[i]
  cihereup <- mean.fs.10yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x+fagsyljitter, 2), y=c(ciherelo, cihereup), col="pink")
  points(pos.x+fagsyljitter, pos.y, cex=cexhere, pch=19, col="pink")
  text(pos.x + 0.03, pos.y, labels=unique(mean.fs.10yr$cc)[i], cex=cextext, col="black")
  }
legend("bottomright", pch=c(19, 19), col=c("salmon", "pink"), legend=c("Betula pendula", "Fagus sylvatica"),
   cex=cexhereleg, bty="n")
dev.off()
## END: Two panel with FagSyl 10 yr
## 


##
## 2 panel with FagSyl 20 yr -- showing raw on top panels, and just logged on bottom panels
fagsyljitter <- 0.04
cexhere <- 0.5
cexhereleg <- 0.7
pdf(file.path("figures/basicpep1950to20002spp2panel.pdf"), width = 5, height = 7)
par(xpd=FALSE)
par(mfrow=c(2,1))
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")))
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.20yr$cc))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slope[i]
  ciherelo <- mean.betpen.20yr$mat60slopeconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.fs.20yr$cc))){
  pos.x <- tempdiffplotfs.20yr[i]
  pos.y <- mean.fs.20yr$mat60slope[i]
  ciherelo <- mean.fs.20yr$mat60slopeconfint11[i]
  cihereup <- mean.fs.20yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x+fagsyljitter, 2), y=c(ciherelo, cihereup), col="dodgerblue")
  points(pos.x+fagsyljitter, pos.y, cex=cexhere, pch=19, col="dodgerblue")
  text(pos.x + 0.03, pos.y, labels=unique(mean.fs.20yr$cc)[i], cex=cextext, col="black")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "dodgerblue"),
   legend=c("Betula pendula", "Fagus sylvatica"), cex=cexhereleg, bty="n")

# Log-log 
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-0.9, 0.1),
     ylab=expression(paste("Estimated sensitivity (log(days)/log(", degree, "C))"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")))
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.20yr$cc))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slopelog[i]
  ciherelo <- mean.betpen.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
for(i in 1:length(unique(mean.fs.20yr$cc))){
  pos.x <- tempdiffplotfs.20yr[i]
  pos.y <- mean.fs.20yr$mat60slopelog[i]
  ciherelo <- mean.fs.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.fs.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x+fagsyljitter, 2), y=c(ciherelo, cihereup), col="pink")
  points(pos.x+fagsyljitter, pos.y, cex=cexhere, pch=19, col="pink")
  text(pos.x + 0.03, pos.y, labels=unique(mean.fs.20yr$cc)[i], cex=cextext, col="black")
  }
legend("bottomright", pch=c(19, 19), col=c("salmon", "pink"), legend=c("Betula pendula", "Fagus sylvatica"),
   cex=cexhereleg, bty="n")
dev.off()
## END: Two panel with FagSyl 20 yr
## 


##
## Four panel with FagSyl 10 yr -- showing logged and raw on left panels, and just logged on right panels
cexhere <- 0.5
cexhereleg <- 0.7
pdf(file.path("figures/basicpep195020102spp4panel.pdf"), width = 9, height = 6)
par(xpd=FALSE)
par(mfrow=c(2,2))
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Betpen")
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.10yr$cc))){
  pos.x <- tempdiffplot.10yr[i]
  pos.y <- mean.betpen.10yr$mat60slope[i]
  ciherelo <- mean.betpen.10yr$mat60slopeconfint11[i]
  cihereup <- mean.betpen.10yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.10yr$cc)[i], cex=cextext, col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.betpen.10yr$cc))){
  pos.x <- tempdiffplot.10yr[i]
  pos.y <- mean.betpen.10yr$mat60slopelog[i]
  ciherelo <- mean.betpen.10yr$mat60slopelogconfint11[i]
  cihereup <- mean.betpen.10yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"),
   legend=c("Using raw x, y", "Using logged x, y"), cex=cexhereleg, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5),  ylim=c(-0.9, 0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Betpen (logged only)")
for(i in 1:length(unique(mean.betpen.10yr$cc))){
  pos.x <- tempdiffplot.10yr[i]
  pos.y <- mean.betpen.10yr$mat60slopelog[i]
  ciherelo <- mean.betpen.10yr$mat60slopelogconfint11[i]
  cihereup <- mean.betpen.10yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.10yr$cc)[i], cex=cextext, col="salmon")
  }
# FagSyl 
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Fagsyl")
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.fs.10yr$cc))){
  pos.x <- tempdiffplotfs.10yr[i]
  pos.y <- mean.fs.10yr$mat60slope[i]
  ciherelo <- mean.fs.10yr$mat60slopeconfint11[i]
  cihereup <- mean.fs.10yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.10yr$cc)[i], cex=cextext, col="darkblue")
  }
for(i in 1:length(unique(mean.fs.10yr$cc))){
  pos.x <- tempdiffplotfs.10yr[i]
  pos.y <- mean.fs.10yr$mat60slopelog[i]
  ciherelo <- mean.fs.10yr$mat60slopelogconfint11[i]
  cihereup <- mean.fs.10yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Using raw x, y", "Using logged x, y"),
   cex=cexhereleg, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5),  ylim=c(-0.9, 0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Fagsyl (logged only)")
for(i in 1:length(unique(mean.fs.10yr$cc))){
  pos.x <- tempdiffplotfs.10yr[i]
  pos.y <- mean.fs.10yr$mat60slopelog[i]
  ciherelo <- mean.fs.10yr$mat60slopelogconfint11[i]
  cihereup <- mean.fs.10yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.10yr$cc)[i], cex=cextext, col="salmon")
  }
dev.off()
## END: Four panel with FagSyl 10 yr -- showing logged and raw on left panels, and just logged on right panels
## 




##
## Four panel with FagSyl 20 yr -- showing logged and raw on left panels, and just logged on right panels
cexhere <- 0.5
cexhereleg <- 0.7
pdf(file.path("figures/basicpep1950to20002spp4panel.pdf"), width = 9, height = 6)
par(xpd=FALSE)
par(mfrow=c(2,2))
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Betpen")
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.20yr$cc))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slope[i]
  ciherelo <- mean.betpen.20yr$mat60slopeconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext, col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.betpen.20yr$cc))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slopelog[i]
  ciherelo <- mean.betpen.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"),
   legend=c("Using raw x, y", "Using logged x, y"), cex=cexhereleg, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5),  ylim=c(-0.9, 0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Betpen (logged only)")
for(i in 1:length(unique(mean.betpen.20yr$cc))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slopelog[i]
  ciherelo <- mean.betpen.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext, col="salmon")
  }
# FagSyl 
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Fagsyl")
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.fs.20yr$cc))){
  pos.x <- tempdiffplotfs.20yr[i]
  pos.y <- mean.fs.20yr$mat60slope[i]
  ciherelo <- mean.fs.20yr$mat60slopeconfint11[i]
  cihereup <- mean.fs.20yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext, col="darkblue")
  }
for(i in 1:length(unique(mean.fs.20yr$cc))){
  pos.x <- tempdiffplotfs.20yr[i]
  pos.y <- mean.fs.20yr$mat60slopelog[i]
  ciherelo <- mean.fs.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.fs.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Using raw x, y", "Using logged x, y"),
   cex=cexhereleg, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5),  ylim=c(-0.9, 0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Fagsyl (logged only)")
for(i in 1:length(unique(mean.fs.20yr$cc))){
  pos.x <- tempdiffplotfs.20yr[i]
  pos.y <- mean.fs.20yr$mat60slopelog[i]
  ciherelo <- mean.fs.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.fs.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext, col="salmon")
  }
dev.off()
## END: Four panel with FagSyl 20yr -- showing logged and raw on left panels, and just logged on right panels
##


##
## Four panel with FagSyl 10 yr -- showing just raw on left panels, and just logged on right panels
cexhere <- 0.5
cexhereleg <- 0.7
cex.mainhere <- 0.85
pdf(file.path("figures/basicpep195020102spp4paneladj.pdf"), width = 7.5, height = 6)
par(xpd=FALSE)
par(mfrow=c(2,2))
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-9, 1), bty="l", mgp=c(1.5, 0.25, 0), tck=-.01,
     ylab=expression(paste("Sensitivity (", italic("Betula pendula"), ")"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Untransformed", cex.main=cex.mainhere)
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.10yr$cc))){
  pos.x <- tempdiffplot.10yr[i]
  pos.y <- mean.betpen.10yr$mat60slope[i]
  ciherelo <- mean.betpen.10yr$mat60slopeconfint11[i]
  cihereup <- mean.betpen.10yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  text(pos.x + 0.2, pos.y, labels=unique(mean.betpen.10yr$cc)[i], cex=cextext)
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5),  ylim=c(-0.7, 0.1), bty="l", mgp=c(1.5, 0.25, 0), tck=-.01,
     ylab=expression(paste("Sensitivity (", italic("Betula pendula"), ")"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Logged", cex.main=cex.mainhere)# (days/", degree, "C)
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.10yr$cc))){
  pos.x <- tempdiffplot.10yr[i]
  pos.y <- mean.betpen.10yr$mat60slopelog[i]
  ciherelo <- mean.betpen.10yr$mat60slopelogconfint11[i]
  cihereup <- mean.betpen.10yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  text(pos.x + 0.2, pos.y, labels=unique(mean.betpen.10yr$cc)[i], cex=cextext)
  }
# legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"),
  #  legend=c("BP: untransformed", "BP: logged"), cex=cexhereleg, bty="n")
# FagSyl 
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-9, 1), bty="l", mgp=c(1.5, 0.25, 0), tck=-.01,
     ylab=expression(paste("Sensitivity (", italic("Fagus sylvatica"), ")"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Untransformed", cex.main=cex.mainhere)
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.fs.10yr$cc))){
  pos.x <- tempdiffplotfs.10yr[i]
  pos.y <- mean.fs.10yr$mat60slope[i]
  ciherelo <- mean.fs.10yr$mat60slopeconfint11[i]
  cihereup <- mean.fs.10yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  text(pos.x + 0.2, pos.y, labels=unique(mean.betpen.10yr$cc)[i], cex=cextext)
  }
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5),  ylim=c(-0.7, 0.1), bty="l", mgp=c(1.5, 0.25, 0), tck=-.01,
     ylab=expression(paste("Sensitivity (", italic("Fagus sylvatica"), ")"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")),  main="Logged", cex.main=cex.mainhere)
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.fs.10yr$cc))){
  pos.x <- tempdiffplotfs.10yr[i]
  pos.y <- mean.fs.10yr$mat60slopelog[i]
  ciherelo <- mean.fs.10yr$mat60slopelogconfint11[i]
  cihereup <- mean.fs.10yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  text(pos.x + 0.2, pos.y, labels=unique(mean.betpen.10yr$cc)[i], cex=cextext)
  }
dev.off()
## END: Four panel with FagSyl 10 yr -- showing just raw on left panels, and just logged on right panels
## 


##
## Four panel with FagSyl 20 yr -- showing just raw on left panels, and just logged on right panels
cexhere <- 0.5
cexhereleg <- 0.7
pdf(file.path("figures/basicpep1950to20002spp4paneladj.pdf"), width = 7.5, height = 6)
par(xpd=FALSE)
par(mfrow=c(2,2))
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-9, 1), bty="l", mgp=c(1.5, 0.25, 0), tck=-.01,
      ylab=expression(paste("Sensitivity (", italic("Betula pendula"), ")"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Untransformed", cex.main=cex.mainhere)
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.20yr$cc))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slope[i]
  ciherelo <- mean.betpen.20yr$mat60slopeconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  text(pos.x + 0.2, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext)
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5),  ylim=c(-0.7, 0.1), bty="l", mgp=c(1.5, 0.25, 0), tck=-.01,
     ylab=expression(paste("Sensitivity (", italic("Betula pendula"), ")"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Logged", cex.main=cex.mainhere)
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen.20yr$cc))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slopelog[i]
  ciherelo <- mean.betpen.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  text(pos.x + 0.2, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext)
  }
# FagSyl 
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-9, 1), bty="l", mgp=c(1.5, 0.25, 0), tck=-.01,
     ylab=expression(paste("Sensitivity (", italic("Fagus sylvatica"), ")"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")),  main="Untransformed", cex.main=cex.mainhere)
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.fs.20yr$cc))){
  pos.x <- tempdiffplotfs.20yr[i]
  pos.y <- mean.fs.20yr$mat60slope[i]
  ciherelo <- mean.fs.20yr$mat60slopeconfint11[i]
  cihereup <- mean.fs.20yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  text(pos.x + 0.2, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext)
  }
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5),  ylim=c(-0.7, 0.1), bty="l", mgp=c(1.5, 0.25, 0), tck=-.01,
    ylab=expression(paste("Sensitivity (", italic("Fagus sylvatica"), ")"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")),  main="Logged", cex.main=cex.mainhere)
abline(h=0, lty=2, col="darkgrey")
mean.fs.20yr$mat60slopelog[2] <- mean.fs.20yr$mat60slopelog[2] - 0.04 # move text label down
for(i in 1:length(unique(mean.fs.20yr$cc))){
  pos.x <- tempdiffplotfs.20yr[i]
  pos.y <- mean.fs.20yr$mat60slopelog[i]
  ciherelo <- mean.fs.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.fs.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  text(pos.x + 0.2, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext)
  }
dev.off()
## END: Four panel with FagSyl 20 yr -- showing just raw on left panels, and just logged on right panels
## 


## Now with the moving windows ...
if(FALSE){
mean.betpenswa <- aggregate(dfswa[c("matslope", "matslopelog", "meanmat")], dfswa["cc"], FUN=mean)
sd.betpenswa <- aggregate(dfswa[c("matslope", "matslopelog", "meanmat")], dfswa["cc"],  FUN=sd)

tempdiffswa <- mean.betpenswa$meanmat[which(mean.betpenswa$cc=="post")]-
    mean.betpenswa$meanmat[which(mean.betpenswa$cc=="pre")]
tempdiffplotswa <- c(0, tempdiffswa)

cexhere <- 0.95
pdf(file.path("figures/basicpepswa.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-3, 1), ylim=c(-8, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpenswa$cc))){
  pos.x <- tempdiffplotswa[i]
  pos.y <- mean.betpenswa$matslope[i]
  sdhere <- sd.betpenswa$matslope[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.betpenswa$cc))){
  pos.x <- tempdiffplotswa[i]
  pos.y <- mean.betpenswa$matslopelog[i]
  sdhere <- sd.betpenswa$matslopelog[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
# par(xpd=TRUE) # so I can plot legend outside
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Using raw x, y", "Using logged x, y"),
   cex=1, bty="n")
dev.off()
}

##############
## Tables ##
##############

mean.betpen.10yr$species <- rep("Betula pendula", nrow(mean.betpen.10yr))
mean.fs.10yr$species <- rep("Fagus sylvatica", nrow(mean.fs.10yr))
mean.qr.10yr$species <- rep("Quercus robur", nrow(mean.qr.10yr))

mean.betpen.forpaper10yr <- subset(mean.betpen.10yr, select=c("cc", "species", 
     "meanmat30", "meanmat45", "meanmat60", "meanmatlo", "varmat30",
     "varmat45", "varmat60", "varlo", "meangdd", "mat30slope", "mat45slope", "mat60slope",
     "mat30slopelog", "mat45slopelog",  "mat60slopelog"))
mean.fs.forpaper10yr <- subset(mean.fs.10yr,  select=c("cc", "species", "meanmat30", "meanmat45", 
    "meanmat60", "meanmatlo", "varmat30", "varmat45", "varmat60", "varlo",
    "meangdd", "mat30slope", "mat45slope", "mat60slope", "mat30slopelog", "mat45slopelog",
    "mat60slopelog"))
mean.qr.forpaper10yr <- subset(mean.qr.10yr,  select=c("cc", "species", "meanmat30", "meanmat45",
    "meanmat60", "meanmatlo", "varmat30", "varmat45", "varmat60", "varlo",
    "meangdd", "mat30slope", "mat45slope", "mat60slope", "mat30slopelog", "mat45slopelog",
     "mat60slopelog"))


mean2spp.forpaper10yr <- rbind(mean.betpen.forpaper10yr, mean.fs.forpaper10yr, mean.qr.forpaper10yr)
mean2spp.forpaper10yr$species <- gsub("Betula pendula", paste0("\\\\emph{","Betula","}"), mean2spp.forpaper10yr$species)
mean2spp.forpaper10yr$species <- gsub("Fagus sylvatica", paste0("\\\\emph{","Fagus","}"), mean2spp.forpaper10yr$species)
mean2spp.forpaper10yr$species <- gsub("Quercus robur", paste0("\\\\emph{","Quercus","}"), mean2spp.forpaper10yr$species)

mean2spp.forpaper10yr$mat30slope <- ifelse(mean2spp.forpaper10yr$mat30slope>=-0.01 & 
    mean2spp.forpaper10yr$mat30slope<=0, 0.0, mean2spp.forpaper10yr$mat30slope)
mean2spp.forpaper10yr$mat45slope <- ifelse(mean2spp.forpaper10yr$mat45slope>=-0.01 & 
    mean2spp.forpaper10yr$mat45slope<=0, 0.0, mean2spp.forpaper10yr$mat45slope)
mean2spp.forpaper10yr$mat60slope <- ifelse(mean2spp.forpaper10yr$mat60slope>=-0.01 & 
    mean2spp.forpaper10yr$mat60slope<=0, 0.0, mean2spp.forpaper10yr$mat60slope)

mean2spp.forpaper10yr$mat30slopelog <- ifelse(mean2spp.forpaper10yr$mat30slopelog>=-0.01 & 
    mean2spp.forpaper10yr$mat30slopelog<=0, 0.0, mean2spp.forpaper10yr$mat30slopelog)
mean2spp.forpaper10yr$mat45slopelog <- ifelse(mean2spp.forpaper10yr$mat45slopelog>=-0.01 & 
    mean2spp.forpaper10yr$mat45slopelog<=0, 0.0, mean2spp.forpaper10yr$mat45slopelog)
mean2spp.forpaper10yr$mat60slopelog <- ifelse(mean2spp.forpaper10yr$mat60slopelog>=-0.01 & 
    mean2spp.forpaper10yr$mat60slopelog<=0, 0.0, mean2spp.forpaper10yr$mat60slopelog)

names(mean2spp.forpaper10yr) <-  c("years", "species", "31", "45", "60", "mean ST.lo", 
    "31", "45", "60", "var (lo)", "GDD", "31", "45", "60", "31", "45", "60")

mean3spp.forpaper10yr <- mean2spp.forpaper10yr
mean2spp.forpaper10yr <- mean2spp.forpaper10yr[1:4,]

mean.betpen.20yr$species <- rep("Betula pendula", nrow(mean.betpen.20yr))
mean.fs.20yr$species <- rep("Fagus sylvatica", nrow(mean.fs.20yr))
mean.qr.20yr$species <- rep("Quercus robur", nrow(mean.qr.20yr))

mean.betpen.forpaper <- subset(mean.betpen.20yr, select=c("cc", "species", 
     "meanmat30", "meanmat45", "meanmat60", "meanmatlo", "varmat30", "varmat45", "varmat60", "varlo",
     "meangdd", "mat30slope", "mat45slope", "mat60slope", "mat30slopelog", "mat45slopelog",  "mat60slopelog"))
mean.fs.forpaper <- subset(mean.fs.20yr,  select=c("cc", "species", 
    "meanmat30", "meanmat45", "meanmat60", "meanmatlo", "varmat30", "varmat45", "varmat60", "varlo",
    "meangdd", "mat30slope", "mat45slope", "mat60slope", "mat30slopelog", "mat45slopelog",  "mat60slopelog"))
mean.qr.forpaper <- subset(mean.qr.20yr,  select=c("cc", "species", 
     "meanmat30", "meanmat45", "meanmat60", "meanmatlo", "varmat30", "varmat45", "varmat60", "varlo",
     "meangdd", "mat30slope", "mat45slope", "mat60slope", "mat30slopelog", "mat45slopelog",  "mat60slopelog"))


mean2spp.forpaper20yr <- rbind(mean.betpen.forpaper, mean.fs.forpaper, mean.qr.forpaper)
mean2spp.forpaper20yr$species <- gsub("Betula pendula", paste0("\\\\emph{","Betula","}"), mean2spp.forpaper20yr$species)
mean2spp.forpaper20yr$species <- gsub("Fagus sylvatica", paste0("\\\\emph{","Fagus","}"), mean2spp.forpaper20yr$species)
mean2spp.forpaper20yr$species <- gsub("Quercus robur", paste0("\\\\emph{","Quercus","}"), mean2spp.forpaper20yr$species)


mean2spp.forpaper20yr$mat30slope <- ifelse(mean2spp.forpaper20yr$mat30slope>=-0.01 & 
    mean2spp.forpaper20yr$mat30slope<=0, 0.0, mean2spp.forpaper20yr$mat30slope)
mean2spp.forpaper20yr$mat45slope <- ifelse(mean2spp.forpaper20yr$mat45slope>=-0.01 & 
    mean2spp.forpaper20yr$mat45slope<=0, 0.0, mean2spp.forpaper20yr$mat45slope)
mean2spp.forpaper20yr$mat60slope <- ifelse(mean2spp.forpaper20yr$mat60slope>=-0.01 & 
    mean2spp.forpaper20yr$mat60slope<=0, 0.0, mean2spp.forpaper20yr$mat60slope)

mean2spp.forpaper20yr$mat30slopelog <- ifelse(mean2spp.forpaper20yr$mat30slopelog>=-0.01 & 
    mean2spp.forpaper20yr$mat30slopelog<=0, 0.0, mean2spp.forpaper20yr$mat30slopelog)
mean2spp.forpaper20yr$mat45slopelog <- ifelse(mean2spp.forpaper20yr$mat45slopelog>=-0.01 & 
    mean2spp.forpaper20yr$mat45slopelog<=0, 0.0, mean2spp.forpaper20yr$mat45slopelog)
mean2spp.forpaper20yr$mat60slopelog <- ifelse(mean2spp.forpaper20yr$mat60slopelog>=-0.01 & 
    mean2spp.forpaper20yr$mat60slopelog<=0, 0.0, mean2spp.forpaper20yr$mat60slopelog)

names(mean2spp.forpaper20yr) <-  c("years", "species", "31", "45", "60", "mean ST.lo", 
    "31", "45", "60", "var (lo)", "GDD", "31", "45", "60", "31", "45", "60")

mean3spp.forpaper20yr <- mean2spp.forpaper20yr
mean2spp.forpaper20yr <- mean2spp.forpaper20yr[1:6,]



