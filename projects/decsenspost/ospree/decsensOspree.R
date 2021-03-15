## Started 24 February 2021 ##
## By Lizzie ##

## Find and plot a controlled multi-force study to use in decsens ##

# housekeeping
rm(list=ls())
options(stringsAsFactors=FALSE)

# Load libraries
library(ggplot2)

# Setting working directory
setwd("~/Documents/git/lab/labgit/projects/decsenspost/ospree")

# get the data (take from ospree repo)
osp <- read.csv("ospree_clean_withchill_BB.csv", header = TRUE)
dall <- subset(osp, datasetID!="junttila12") # removing junttilla which I noted is a dormancy release study
# do some cleanup of this unwieldy file
source("commoncols.R")
d <- dall[,which(colnames(dall) %in% c(common.cols.wchill, "forcetemp_night", "photoperiod_night"))]
d$datasetstudy <- paste(d$datasetID, d$study)
d$forceday <- as.numeric(d$forcetemp)
d$forcenight <- as.numeric(d$forcetemp_night)
d$photonight <- as.numeric(d$photoperiod_night)

d$chilldaysnum <- as.numeric(d$chilldays)
char2 <- d[which((d$datasetID %in% "charrier11") & (d$study %in% "exp2")),]

library(ggplot2)
library(gridExtra)

##
##

# For manuscript

cex.mainhere <- 0.9
cex.labhere <- 1.1

pdf("..//figures/ospreeforcems.pdf", width=7, height=3.25)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
par(mfrow=c(1,2))
plot(response.time~forceday, data=char2, pch=19,
     ylab="Days to budburst",
     xlab=expression(paste("Temperature (", degree, "C)")),
     main="Linear (untransformed)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=cex.labhere,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01, xlim=c(0.1, 26))
plot(response.time~forceday, data=char2, pch=19, log="xy", 
     ylab="Days to budburst",
     xlab=expression(paste("Temperature (", degree, "C)")),
     main="Non-linear (logged x and y)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=cex.labhere,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
dev.off()





