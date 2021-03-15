## Started 25 January 2021 ##
## By Lizzie ##

## Some code taken from analyseataglance.R 

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/projects/treegarden/decsens/analyses") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("hereliesboomboom")

library(ggplot2)

setwd("~/Documents/git/lab/labgit/projects/decsenspost/winegrape_analyses")

# climate data from paper (received from Ben Cook)
temp <- read.csv("seas_temp_MJJ.onedeg.csv", header=TRUE)
colnames(temp)[1] <- "year"
# phenology data
daux <- read.csv("dauxdata.csv", header=TRUE, skip=2)
daux <- daux[,1:28]
names(daux)[names(daux)=="Abb."] <- "year"
ghd <- as.data.frame(daux)
ghd <- subset(ghd, select=c("year", "Bor", "Bur")) # selecting these two famous region, they are high-quality, complete data

envdata <- merge(ghd, temp, by="year", suffixes=c(".ghd", ".temp"))

## Setting up data to have same years 1980+ and before
subearly <- subset(envdata, year>1950 & year <1980)
sublate <- subset(envdata, year>1979)

# Making sure we have fairly even dataset lengths
nrow(subearly)
nrow(sublate)

burearlylm <- lm(Bur.ghd~Bur.temp, data=subearly)
burlatelm <- lm(Bur.ghd~Bur.temp, data=sublate)

# Note: These data are anomalized to 31 August, I am adding 60 because the GHD can be negative and log-transform does not allow that (effectively I am anomalizing to early July, affects only the intercept)
burearlylog <- lm(log(Bur.ghd+60)~log(Bur.temp), data=subearly) 
burlatelog <-  lm(log(Bur.ghd+60)~log(Bur.temp), data=sublate)

bordearlylm <- lm(Bor.ghd~Bor.temp, data=subearly)
bordlatelm <- lm(Bor.ghd~Bor.temp, data=sublate)

bordearlylog <- lm(log(Bor.ghd+60)~log(Bor.temp), data=subearly) 
bordlatelog <- lm(log(Bor.ghd+60)~log(Bor.temp), data=sublate)

# Ugly code to make up dataframes for plotting
burearly <- data.frame(where="Burgundy",
                    when="1951-1979",
                    meantemp=mean(subearly$Bur.temp, na.rm=TRUE),
                    lmslope=coef(burearlylm)[2],
                    lmslopeconfint11=confint(burearlylm,level=0.89)[2,1],
                    lmslopeconfint89=confint(burearlylm,level=0.89)[2,2],
                    logslope=coef(burearlylog)[2],
                    logslopeconfint11=confint(burearlylog,level=0.89)[2,1],
                    logslopeconfint89=confint(burearlylog,level=0.89)[2,2])

burlate <- data.frame(where="Burgundy",
                    when="1980-2007",
                    meantemp=mean(sublate$Bur.temp, na.rm=TRUE),
                    lmslope=coef(burlatelm)[2],
                    lmslopeconfint11=confint(burlatelm,level=0.89)[2,1],
                    lmslopeconfint89=confint(burlatelm,level=0.89)[2,2],
                    logslope=coef(burlatelog)[2],
                    logslopeconfint11=confint(burlatelog,level=0.89)[2,1],
                    logslopeconfint89=confint(burlatelog,level=0.89)[2,2])

bordearly <- data.frame(where="Bordeaux",
                    when="1951-1979",
                    meantemp=mean(subearly$Bor.temp, na.rm=TRUE),
                    lmslope=coef(bordearlylm)[2],
                    lmslopeconfint11=confint(bordearlylm,level=0.89)[2,1],
                    lmslopeconfint89=confint(bordearlylm,level=0.89)[2,2],
                    logslope=coef(bordearlylog)[2],
                    logslopeconfint11=confint(bordearlylog,level=0.89)[2,1],
                    logslopeconfint89=confint(bordearlylog,level=0.89)[2,2])
bordlate <- data.frame(where="Bordeaux",
                    when="1980-2007",
                    meantemp=mean(sublate$Bor.temp, na.rm=TRUE),
                    lmslope=coef(bordlatelm)[2],
                    lmslopeconfint11=confint(bordlatelm,level=0.89)[2,1],
                    lmslopeconfint89=confint(bordlatelm,level=0.89)[2,2],
                    logslope=coef(bordlatelog)[2],
                    logslopeconfint11=confint(bordlatelog,level=0.89)[2,1],
                    logslopeconfint89=confint(bordlatelog,level=0.89)[2,2])

plotburdat <- rbind(burearly, burlate)
plotborddat <- rbind(bordearly, bordlate)

tempdiffbur <- plotburdat$meantemp[which(plotburdat$when=="1980-2007")]-
    plotburdat$meantemp[which(plotburdat$when=="1951-1979")]
tempdiffburplot <- c(-0.01, tempdiffbur)

tempdiffbord <- plotborddat$meantemp[which(plotborddat$when=="1980-2007")]-
    plotborddat$meantemp[which(plotborddat$when=="1951-1979")]
tempdiffbordplot <- c(0.01, tempdiffbord)


cexhere <- 0.95
cextext <- 0.75
colz <- c("tomato2", "maroon4")
pdf(file.path("..//figures/winedata.pdf"), width = 9, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
par(mfrow=c(1,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.7), ylim=c(-10, -4),
     ylab=expression(paste("Estimated sensitivity"), sep=""), #  (days/", degree, "C)
         xlab=expression(paste("Warming (", degree, "C)")), main="Linear (untransformed)",
     font.main = 1, cex.main = 0.9, cex.lab=1.2, bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
for(i in 1:length(unique(plotburdat$when))){
  pos.x <- tempdiffburplot[i]
  pos.y <- plotburdat$lmslope[i]
  ciherelo <- plotburdat$lmslopeconfint11[i]
  cihereup <- plotburdat$lmslopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col=colz[1])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colz[1])
  text(pos.x + 0.22, pos.y, labels=unique(plotburdat$when)[i], cex=cextext, col=colz[1])
  }
for(i in 1:length(unique(plotborddat$when))){
  pos.x <- tempdiffbordplot[i]
  pos.y <- plotborddat$lmslope[i]
  ciherelo <- plotborddat$lmslopeconfint11[i]
  cihereup <- plotborddat$lmslopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col=colz[2])
  points(pos.x, pos.y, cex=cexhere, pch=17, col=colz[2])
  text(pos.x + 0.24, pos.y, labels=unique(plotborddat$when)[i], cex=cextext, col=colz[2])
  }
legend("bottomright", pch=c(19, 17), col=colz, legend=c("Burgundy", "Bordeaux"),
   cex=1, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.7), ylim=c(-2.5, -0.5),
     ylab=expression(paste("Estimated sensitivity"), sep=""), # (days/", degree, "C)
         xlab=expression(paste("Warming (", degree, "C)")), main="Non-linear (logged)",
     font.main = 1, cex.main = 0.9, cex.lab=1.2, bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
for(i in 1:length(unique(plotburdat$when))){
  pos.x <- tempdiffburplot[i]
  pos.y <- plotburdat$logslope[i]
  ciherelo <- plotburdat$logslopeconfint11[i]
  cihereup <- plotburdat$logslopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col=colz[1])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colz[1])
 text(pos.x + 0.22, pos.y, labels=unique(plotburdat$when)[i], cex=cextext, col=colz[1])
  }
for(i in 1:length(unique(plotborddat$when))){
  pos.x <- tempdiffbordplot[i]
  pos.y <- plotborddat$logslope[i]
  ciherelo <- plotborddat$logslopeconfint11[i]
  cihereup <- plotborddat$logslopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col=colz[2])
  points(pos.x, pos.y, cex=cexhere, pch=17, col=colz[2])
  text(pos.x + 0.24, pos.y-0.1, labels=unique(plotborddat$when)[i], cex=cextext, col=colz[2])
  }
dev.off()
