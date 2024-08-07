## Started 13 Jan 2020
## By Cat - based off code by Lizzie in ospree repo ospree/analyses/bb_analysis/pep_sims/comparetopepsims.R

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

#Libraries
library(dplyr)
library(tidyr)

# Setting working directory. Add in your own path in an if statement for your file structure
setwd("~/Documents/git/decsens/analyses/pep_analyses") # setwd("~/Documents/git/projects/treegarden/decsens/analyses/pep_analyses")

# get some data
# Betula puendula data from PEP (both have has GDD from 1 Jan to leafout)
# bp has mat from March 1st to June 1st and mat.lo is 30 days before leafout (uses tg -- aka mean -- data from E-OBS)

#bp <- read.csv("output/betpen_decsens_1950_2000.csv")  ### 10 yr timeframes (1950-1960 vs 2000-2010)
bp <- read.csv("output/betpen_decsens_1950-2000.csv")  ### 20 yr comparisons (1950-1970 vs 1970-1990 vs 1990-2010)

# loop to extract some model estimates
# this takes mean for each time period then allows comparison acrosgs the two resulting values
bpest <- data.frame(siteslist=numeric(), cc=character(), 
                    meanmat30=numeric(), varmat30=numeric(),  varlogmat30=numeric(), sdmat30=numeric(), 
                    meanmat45=numeric(), varmat45=numeric(),  varlogmat45=numeric(), sdmat45=numeric(), 
                    meanmat60=numeric(), varmat60=numeric(),  varlogmat60=numeric(), sdmat60=numeric(), 
                    meanlo=numeric(), varlo=numeric(), varloglo=numeric(), sdlo=numeric(), meanutah=numeric(), meangdd=numeric(), 
                    mat30slope=numeric(), mat30slopese=numeric(), mat30slopeconfint11=numeric(), mat30slopeconfint89=numeric(),
                    mat45slope=numeric(), mat45slopese=numeric(), mat45slopeconfint11=numeric(), mat45slopeconfint89=numeric(),
                    mat60slope=numeric(), mat60slopese=numeric(), mat60slopeconfint11=numeric(), mat60slopeconfint89=numeric(),
                    meanmatlo=numeric(), 
                    mat30slopelog=numeric(), mat30slopelogse=numeric(), mat30slopelogconfint11=numeric(), mat30slopelogconfint89=numeric(),
                    mat45slopelog=numeric(), mat45slopelogse=numeric(), mat45slopelogconfint11=numeric(), mat45slopelogconfint89=numeric(),
                    mat60slopelog=numeric(), mat60slopelogse=numeric(), mat60slopelogconfint11=numeric(), mat60slopelogconfint89=numeric(),
                    lomatslope=numeric(), lomatslopelog=numeric(), lomatslopelogse=numeric(), lomatslopelogconfint11=numeric(),
                    lomatslopelogconfint89=numeric(),
                    varmatlo=numeric(), sdmatlo=numeric())

sitez <- unique(bp$siteslist)

for(i in c(1:length(sitez))){ # i <- 1
  subby <- subset(bp, siteslist==sitez[i])
  for(ccstate in c(1:3)){ ## ccstate=1
    subbycc <- subset(subby, cc==unique(bp$cc)[ccstate])
    meanmat30 <- mean(subbycc$mat30, na.rm=TRUE)
    varmat30 <- var(subbycc$mat30, na.rm=TRUE)
    varlogmat30 <- var(log(subbycc$mat30), na.rm=TRUE)
    sdmat30 <- sd(subbycc$mat30, na.rm=TRUE)
    meanmat45 <- mean(subbycc$mat45, na.rm=TRUE)
    varmat45 <- var(subbycc$mat45, na.rm=TRUE)
    varlogmat45 <- var(log(subbycc$mat45), na.rm=TRUE)
    sdmat45 <- sd(subbycc$mat45, na.rm=TRUE)
    meanmat60 <- mean(subbycc$mat60, na.rm=TRUE)
    varmat60 <- var(subbycc$mat60, na.rm=TRUE)
    varlogmat60 <- var(log(subbycc$mat60), na.rm=TRUE)
    sdmat60 <- sd(subbycc$mat60, na.rm=TRUE)
    meanmatlo <- mean(subbycc$mat.lo, na.rm=TRUE)
    varmatlo <- var(subbycc$mat.lo, na.rm=TRUE)
    sdmatlo <- sd(subbycc$mat.lo, na.rm=TRUE)
    meanlo <- mean(subbycc$lo, na.rm=TRUE)
    varlo <- var(subbycc$lo, na.rm=TRUE)
    varloglo <- var(log(subbycc$lo), na.rm=TRUE)
    sdlo <- sd(subbycc$lo, na.rm=TRUE)
    meanutah <- mean(subbycc$chillutah, na.rm=TRUE)
    meangdd <- mean(subbycc$gdd, na.rm=TRUE)
    lmmat30 <- lm(lo~mat30, data=subbycc)
    lmmat30se <- summary(lmmat30)$coef[2,2]
    lmmat30confint11 <- confint(lmmat30,level=0.89)[2,1]
    lmmat30confint89 <- confint(lmmat30,level=0.89)[2,2]
    lmmat30log <- lm(log(lo)~log(mat30), data=subbycc)
    lmmat30logse <- summary(lmmat30log)$coef[2,2]
    lmmat30confintlog11 <- confint(lmmat30log,level=0.89)[2,1]
    lmmat30confintlog89 <- confint(lmmat30log,level=0.89)[2,2]
    lmmat45 <- lm(lo~mat45, data=subbycc)
    lmmat45se <- summary(lmmat45)$coef[2,2]
    lmmat45confint11 <- confint(lmmat45,level=0.89)[2,1]
    lmmat45confint89 <- confint(lmmat45,level=0.89)[2,2]
    lmmat45log <- lm(log(lo)~log(mat45), data=subbycc)
    lmmat45logse <- summary(lmmat45log)$coef[2,2]
    lmmat45confintlog11 <- confint(lmmat45log,level=0.89)[2,1]
    lmmat45confintlog89 <- confint(lmmat45log,level=0.89)[2,2]
    lmmat60 <- lm(lo~mat60, data=subbycc)
    lmmat60se <- summary(lmmat60)$coef[2,2]
    lmmat60confint11 <- confint(lmmat60,level=0.89)[2,1]
    lmmat60confint89 <- confint(lmmat60,level=0.89)[2,2]
    lmmat60log <- lm(log(lo)~log(mat60), data=subbycc)
    lmmat60logse <- summary(lmmat60log)$coef[2,2]
    lmmat60confintlog11 <- confint(lmmat60log,level=0.89)[2,1]
    lmmat60confintlog89 <- confint(lmmat60log,level=0.89)[2,2]
    lolmmat <- lm(lo~mat.lo, data=subbycc)
    lolmmatlog <- lm(log(lo)~log(mat.lo), data=subbycc)
    lolmmatlogse <- summary(lolmmatlog)$coef[2,2]
    lolmmatconfintlog11 <- confint(lolmmatlog,level=0.89)[2,1]
    lolmmatconfintlog89 <- confint(lolmmatlog,level=0.89)[2,2]
    bpestadd <- data.frame(siteslist=sitez[i], cc=unique(bp$cc)[ccstate], 
                           meanmat30=meanmat30, varmat30=varmat30, varlogmat30=varlogmat30, sdmat30=sdmat30, 
                           meanmat45=meanmat45, varmat45=varmat45, varlogmat45=varlogmat45, sdmat45=sdmat45, 
                           meanmat60=meanmat60, varmat60=varmat60, varlogmat60=varlogmat60, sdmat60=sdmat60, 
                           meanlo=meanlo, varlo=varlo, varloglo=varloglo, sdlo=sdlo, meanutah=meanutah, 
                           meangdd=meangdd, 
                           mat30slope=coef(lmmat30)["mat30"], mat30slopese=lmmat30se, mat30slopeconfint11=lmmat30confint11, 
                           mat30slopeconfint89=lmmat30confint89,
                           mat45slope=coef(lmmat45)["mat45"], mat45slopese=lmmat45se, mat45slopeconfint11=lmmat45confint11, 
                           mat45slopeconfint89=lmmat45confint89,
                           mat60slope=coef(lmmat60)["mat60"], mat60slopese=lmmat60se, mat60slopeconfint11=lmmat60confint11, 
                           mat60slopeconfint89=lmmat60confint89,
                           mat30slopelog=coef(lmmat30log)["log(mat30)"], mat30slopelogse=lmmat30logse, mat30slopelogconfint11=lmmat30confintlog11, 
                           mat30slopelogconfint89=lmmat30confintlog89,
                           mat45slopelog=coef(lmmat45log)["log(mat45)"], mat45slopelogse=lmmat45logse, mat45slopelogconfint11=lmmat45confintlog11, 
                           mat45slopelogconfint89=lmmat45confintlog89,
                           mat60slopelog=coef(lmmat60log)["log(mat60)"], mat60slopelogse=lmmat60logse, mat60slopelogconfint11=lmmat60confintlog11, 
                           mat60slopelogconfint89=lmmat60confintlog89,
                           lomatslope=coef(lolmmat)["mat.lo"],
                           lomatslopelog=coef(lolmmatlog)["log(mat.lo)"], lomatslopelogse=lolmmatlogse, lomatslopelogconfint11=lolmmatconfintlog11, 
                           lomatslopelogconfint89=lolmmatconfintlog89,
                           meanmatlo=meanmatlo,
                           varmatlo=varmatlo, sdmatlo=sdmatlo)
    bpest <- rbind(bpest, bpestadd)
  }
}    

meanhere <- aggregate(bpest[c("meanmat30", "varmat30", "varlogmat30", "sdmat30", "meanmat45", "varmat45", "varlogmat45", "sdmat45", "meanmat60", "varmat60", "varlogmat60", "sdmat60",
                              "meanmatlo", "varmatlo", "sdmatlo", "meanlo", "varlo", "varloglo", "sdlo", "meanutah", "meangdd",
                              "mat30slope", "mat30slopese", "mat30slopeconfint11", "mat30slopeconfint89",  "mat30slopelog", "mat30slopelogse", "mat30slopelogconfint11", "mat30slopelogconfint89",
                              "mat45slope", "mat45slopese", "mat45slopeconfint11", "mat45slopeconfint89",  "mat45slopelog", "mat45slopelogse", "mat45slopelogconfint11", "mat45slopelogconfint89",
                              "mat60slope", "mat60slopese", "mat60slopeconfint11", "mat60slopeconfint89",  "mat60slopelog", "mat60slopelogse", "mat60slopelogconfint11", "mat60slopelogconfint89",
                             "lomatslope", "lomatslopelog", "lomatslopelogse")], bpest["cc"], FUN=mean)
sdhere <- aggregate(bpest[c("meanmat30", "varmat30", "varlogmat30", "meanmat45", "varmat45", "varlogmat45", "meanmat60", "varmat60", "varlogmat60",
                            "meanmatlo", "varmatlo", "meanlo", "varlo", "varloglo", "meanutah", "meangdd", "mat30slope", "mat45slope", "mat60slope")],
                    bpest["cc"], FUN=sd)


#      cc     meanmat    varmat  varlogmat     sdmat meanmatlo varmatlo   sdmatlo  meanlo    varlo    varloglo
# 1950-1970 7.661285 1.2515107 0.02308514 1.1172600  7.267610 1.076558 0.9952432 114.926 81.65211 0.006449103
# 1990-2010 8.807405 0.8897012 0.01224291 0.9427715  6.600241 1.010194 0.9970972 106.422 42.46895 0.003875612
#     sdlo meanutah  meangdd  matslope matslopese matslopeconfint11 matslopeconfint89 matslopelog matslopelogse
# 8.914409 2038.200 77.86681 -6.042843   1.204039         -8.067031        -4.0186559  -0.3978114    0.08225647
# 6.363890 2287.104 59.19161 -2.291221   1.481032         -4.781080         0.1986369  -0.1839502    0.12221731
#      matslopelogconfint11 matslopelogconfint89
#           -0.5360981          -0.25952477
#           -0.3894176           0.02151724

### Compare logged to unlogged 
bpest$mat30slopelog_exp <- exp(bpest$mat30slopelog)
mean(bpest$mat30slopelog_exp) ## 1.02
bpest$mat45slopelog_exp <- exp(bpest$mat45slopelog)
mean(bpest$mat45slopelog_exp) ## 0.96
bpest$mat60slopelog_exp <- exp(bpest$mat60slopelog)
mean(bpest$mat60slopelog_exp) ## 0.82

mean(bpest$mat30slope) ## -0.36
mean(bpest$mat45slope) ## -0.89
mean(bpest$mat60slope) ## -4.00

write.csv(bpest, file="output/bpenestimates_withlog_1950to2010.csv", row.names = FALSE)
#write.csv(bpest, file="output/bpenestimates_withlog_1950_2000.csv", row.names = FALSE)
write.csv(meanhere, file="output/bpestimates_twentyyrwindows.csv", row.names = FALSE)
#write.csv(meanhere, file="output/bpestimates_tenyrwindows.csv", row.names = FALSE)

## Now do as above, but for 10-year windows  ....
bp10yr <- bp
bp10yr$cc[which(bp10yr$year>1950 & bp10yr$year<1961)] <- "1950-1960"
bp10yr$cc[which(bp10yr$year>1960 & bp10yr$year<1971)] <- "1960-1970"
bp10yr$cc[which(bp10yr$year>1970 & bp10yr$year<1981)] <- "1970-1980"
bp10yr$cc[which(bp10yr$year>1980 & bp10yr$year<1991)] <- "1980-1990"
bp10yr$cc[which(bp10yr$year>1990 & bp10yr$year<2001)] <- "1990-2000"
bp10yr$cc[which(bp10yr$year>2000 & bp10yr$year<2011)] <- "2000-2010"

bpest.10yr <- data.frame(siteslist=numeric(), cc=character(), meanmat=numeric(), varmat=numeric(),  
                    sdmat=numeric(), meanlo=numeric(), varlo=numeric(), sdlo=numeric(), meanutah=numeric(), meangdd=numeric(), 
                    matslope=numeric(), matslopese=numeric(), matslopeconfint11=numeric(), matslopeconfint89=numeric(),
                    meanmatlo=numeric(), 
                    matslopelog=numeric(), matslopelogse=numeric(), matslopelogconfint11=numeric(), matslopelogconfint89=numeric(),
                    varmatlo=numeric(), sdmatlo=numeric())

sitez <- unique(bp10yr$siteslist)

for(i in c(1:length(sitez))){ # i <- 1
  subby <- subset(bp10yr, siteslist==sitez[i])
  for(ccstate in c(1:6)){ ## ccstate=1
    subbycc <- subset(subby, cc==unique(bp10yr$cc)[ccstate])
    meanmat <- mean(subbycc$mat, na.rm=TRUE)
    varmat <- var(subbycc$mat, na.rm=TRUE)
    sdmat <- sd(subbycc$mat, na.rm=TRUE)
    meanmatlo <- mean(subbycc$mat.lo, na.rm=TRUE)
    varmatlo <- var(subbycc$mat.lo, na.rm=TRUE)
    sdmatlo <- sd(subbycc$mat.lo, na.rm=TRUE)
    meanlo <- mean(subbycc$lo, na.rm=TRUE)
    varlo <- var(subbycc$lo, na.rm=TRUE)
    sdlo <- sd(subbycc$lo, na.rm=TRUE)
    meanutah <- mean(subbycc$chillutah, na.rm=TRUE)
    meangdd <- mean(subbycc$gdd, na.rm=TRUE)
    lmmat <- lm(lo~mat, data=subbycc)
    lmmatse <- summary(lmmat)$coef[2,2]
    lmmatconfint11 <- confint(lmmat,level=0.89)[2,1]
    lmmatconfint89 <- confint(lmmat,level=0.89)[2,2]
    lmmatlog <- lm(log(lo)~log(mat), data=subbycc)
    lmmatlogse <- summary(lmmatlog)$coef[2,2]
    lmmatconfintlog11 <- confint(lmmatlog,level=0.89)[2,1]
    lmmatconfintlog89 <- confint(lmmatlog,level=0.89)[2,2]
    bpestadd.10yr <- data.frame(siteslist=sitez[i], cc=unique(bp10yr$cc)[ccstate], meanmat=meanmat, 
                           varmat=varmat, sdmat=sdmat, meanlo=meanlo, varlo=varlo, sdlo=sdlo, meanutah=meanutah, 
                           meangdd=meangdd, matslope=coef(lmmat)["mat"], matslopese=lmmatse, matslopeconfint11=lmmatconfint11, 
                           matslopeconfint89=lmmatconfint89,
                           matslopelog=coef(lmmatlog)["log(mat)"], matslopelogse=lmmatlogse, matslopelogconfint11=lmmatconfintlog11, 
                           matslopelogconfint89=lmmatconfintlog89,
                           meanmatlo=meanmatlo,
                           varmatlo=varmatlo, sdmatlo=sdmatlo)
    bpest.10yr <- rbind(bpest.10yr, bpestadd.10yr)
  }
}    

meanhere.10yr <- aggregate(bpest.10yr[c("meanmat", "varmat", "sdmat", "meanmatlo", "varmatlo", "sdmatlo", "meanlo",
    "varlo", "sdlo", "meanutah", "meangdd","matslope", "matslopese", "matslopelog", "matslopelogse")], bpest.10yr["cc"], FUN=mean)
                              
sdhere.10yr <- aggregate(bpest.10yr[c("meanmat", "varmat", "meanmatlo", "varmatlo", "meanlo", "varlo", "meanutah", "meangdd", "matslope")],
    bpest.10yr["cc"], FUN=sd)



## Also get the difference for each site across two time periods
# This is to compare to sims better

bpest.sitediffs <- data.frame(siteslist=numeric(), matdiff=numeric(), matlodiff=numeric(), diffslope=numeric(),
                              varlodiff=numeric(), varlodiffper=numeric(), varmatdiffper=numeric())

for(i in c(1:length(sitez))){ # i <- 1
  subby <- subset(bpest, siteslist==sitez[i])
  precc <- subset(subby, cc=="1950-1960")
  postcc <- subset(subby, cc=="2000-2010")
  matdiff <- precc$meanmat-postcc$meanmat
  matlodiff <- precc$meanmatlo-postcc$meanmatlo
  diffslope <- precc$matslope-postcc$matslope
  diffslopelog <- precc$matslopelog-postcc$matslopelog
  varlodiff <- precc$varlo-postcc$varlo
  varlodiffper <- postcc$varlo/precc$varlo
  varmatdiffper <- postcc$varmat/precc$varmat
  bpest.sitediffs.add <- data.frame(siteslist=sitez[i], matdiff=matdiff,matlodiff=matlodiff, diffslope=diffslope,
                                    diffslopelog,
                                    varlodiff=varlodiff, varlodiffper=varlodiffper, varmatdiffper=varmatdiffper)
  bpest.sitediffs <- rbind(bpest.sitediffs, bpest.sitediffs.add)
}

bpest.sitediffs$daysperC <- bpest.sitediffs$diffslope/bpest.sitediffs$matdiff
bpest.sitediffs$daysperClog <- bpest.sitediffs$diffslopelog/bpest.sitediffs$matdiff


########################
## Stuff in the supp ##
########################

# estimate change in variance
# take the means, then calculates the percent change between pre and post cc 
1-meanhere$varlo[2]/meanhere$varlo[1]
1-mean.pepsims$var.lo.postcc[1]/mean.pepsims$var.lo.precc[1]
1-median.pepsims$var.lo.postcc[1]/median.pepsims$var.lo.precc[1]

# in supp 'we estimated a decline in sensitivity'
mean(bpest.sitediffs$daysperC)
sd(bpest.sitediffs$daysperC)/sqrt(45) # days per C mean SE bp

# in supp 'given X warming'
mean(bpest.sitediffs$matdiff)
mean(bpest.sitediffs$matdiff)/sqrt(45)

# compare to pep sims (in figure, not text in supp currently)
mean.pepsims$diffbefore.after[1]
sd.pepsims$diffbefore.after[1]/(sqrt(45)) # days per C mean SE pepsims

# Utah units
meanhere$meanutah
sdhere$meanutah/sqrt(45)

# gdd and matlo
meanhere$meangdd
sdhere$meangdd/sqrt(45)
meanhere$meanmatlo
sdhere$meanmatlo/sqrt(45)

##############################
## End of stuff in the supp ##
##############################

# variance! here's what you get if you calculate the % change site-by-site then average
1 - mean(bpest.sitediffs$varlodiffper)
sd(bpest.sitediffs$varlodiffper)/sqrt(45)
1 - mean.pepsims$varlodiffper[1]
sd.pepsims$varlodiffper[1]/(sqrt(45))

# closer look at variance differences due to when we take the mean ...
# when we take the diff at each site, we exacerbate outliers, this explains (I think) the difference between the two ways we calculate variance
pepsims.1d <- subset(pepsims, degwarm==1)
hist(pepsims.1d$varlodiff)
hist(pepsims.1d$var.lo.precc)
hist(pepsims.1d$var.lo.postcc)

# other stuff ... (not using currently)
mean(bpest.sitediffs$diffslope)
sd(bpest.sitediffs$diffslope)/sqrt(45)
mean(bpest.sitediffs$matlodiff)
sd(bpest.sitediffs$matdiff)
sd(bpest.sitediffs$matlodiff)
sd(bpest.sitediffs$daysperC)

##############
## Plotting ##
##############

cexhere <- 0.95
pdf(file.path("figures/peprealandsims.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(0.5,4.5), ylim=c(-3.1, -0.1),
     ylab=expression(paste("Change in estimated sensitivity (days/", degree, "C)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:4){
  pos.x <- mean.pepsims$degwarm[i]
  pos.y <- mean.pepsims$diffbefore.after[i]
  sehere <- sd.pepsims$diffbefore.after[i]/(sqrt(45))
  lines(x=rep(pos.x, 2), y=c(pos.y-sehere, pos.y+sehere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
}
realdat.diff <- mean(bpest.sitediffs$diffslope) 
points(abs(mean(bpest.sitediffs$matdiff)), realdat.diff, cex=cexhere, pch=17, col="salmon")
realdatse <- sd(bpest.sitediffs$diffslope)/sqrt(45)
lines(x=rep(abs(mean(bpest.sitediffs$matdiff)), 2), y=c(realdat.diff-realdatse, realdat.diff+realdatse),
      col="salmon")
# par(xpd=TRUE) # so I can plot legend outside
legend("topright", pch=c(17, 19), col=c("salmon", "darkblue"), legend=c("European data (PEP725)", "Simulations with constant cues"),
       cex=1, bty="n")
dev.off()
