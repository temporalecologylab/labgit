## Started 15 May 2020 ##
## By Lizzie  ##

## Plotting daily climate data ##


# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
require(plyr)
require(dplyr)
require(ggplot2)

setwd("~/Documents/git/projects/treegarden/decsens/analyses/pep_analyses/")

##
## Fall 2020, checking alpha_1 linearity assumption
##
# From Cat: daily temps from Jan 1 to Apr 30. Let me know if you'd like me to make a figure! I used the 17 sites from the BETPEN 20-year window
bp <- read.csv("output/dailytemps_jantoapr.csv")
bp$date <- as.Date(bp$Date, format="%Y-%m-%d")
bp$doy <- format(bp$date, "%j")
bp$mon <- format(bp$date, "%m")
bp$lat.round <- round(bp$lat, 2)
bp$lon.round <- round(bp$lon, 2)
bp$latlon.round <- paste(bp$lat.round, "N", bp$lon.round, "E", sep=" ")
length(unique(bp$lat.long))
length(unique(bp$latlon.round))

bpsm <- subset(bp, as.numeric(doy)>45)

bpsm.select <- subset(bpsm, year==2007)
ggplot(bpsm.select, aes(x=as.numeric(doy), y=Tavg, group=as.factor(year), colour=as.factor(year))) +
    geom_point() +
   # geom_smooth(method="lm") + 
    facet_wrap(.~as.factor(latlon.round))

bpsm$decade <- NA
bpsm$decade[which(bpsm$year<1971)] <- "1951-1970"
bpsm$decade[which(bpsm$year>1990)] <- "1991-2010"
bpsm$decade[which(is.na(bpsm$decade)==TRUE)] <- "1971-1990"

bpsummdec <-
      ddply(bpsm, c("latlon.round", "decade", "doy"), summarise,
      temp = mean(Tavg))

ggplot(bpsummdec, aes(x=as.numeric(doy), y=temp, group=as.factor(decade), color=as.factor(decade))) +
    geom_line() +
    # geom_smooth(method="lm") +
    facet_wrap(.~as.factor(latlon.round))

bpsumm <-
      ddply(bpsm, c("latlon.round", "doy"), summarise,
      temp = mean(Tavg))

getrsq <- function(df, x, y, group){
    dfnew <- data.frame(group=character(), r2=numeric(), slope=numeric())
    getgroups <- unlist(unique(df[group]), use.names=FALSE)
    for(agroup in getgroups){
    subby <- df[which(df[group]==agroup),]
    m <- lm(subby[[y]] ~ as.numeric(subby[[x]]))
    r2 = format(summary(m)$r.squared, digits = 3)
    dfhere <-  data.frame(group=agroup, r2=r2, slope=coef(m)[2])
    dfnew <- rbind(dfnew, dfhere)
}
return(dfnew)
}

bpsummr2 <- getrsq(bpsumm, "doy", "temp", "latlon.round")
names(bpsummr2)[names(bpsummr2)=="group"] <- "latlon.round"

# Jonathan's suggested visualization: https://statmodeling.stat.columbia.edu/2014/04/10/small-multiples-lineplots-maps-ok-always-yes-case/
ggplot(bpsumm, aes(x=as.numeric(doy), y=temp)) +
    geom_line(color="dodgerblue") +
    geom_smooth(method = "lm", linetype = 2, lwd=0.5, color="darkgray", se = FALSE) +
    facet_wrap(.~as.factor(latlon.round)) +
    geom_text(color="dodgerblue", size=3, data=bpsummr2, aes(x = 57, y = 12, label = r2)) +
    xlab("Day of year") +
    ylab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()



## tried to see if moving average helped with visualization... no much
if(FALSE){
require(zoo)
bpma <- bpsumm %>% 
  group_by(lat.long) %>% 
  mutate(ma7day = rollmean(temp, 2, na.pad = TRUE))
ggplot(bpma, aes(x=as.numeric(doy), y=ma7day, group=as.factor(decade), colour=as.factor(decade))) +
    geom_line() +
    facet_wrap(.~as.factor(lat.long))
}


##
## Work from Spring 2020 
##

d <- read.csv("output/betpen_dailytempsandlo_1950to2010.csv")
d$date <- as.Date(d$Date, format="%Y-%m-%d")
d$mon <- format(d$date, "%m")
d$doy <- format(d$date, "%j")
d$cc <- ifelse(d$year<1981, "pre", "post")


ggplot(d, aes(x=doy, y=tmean, group=year, colour=year)) +
    geom_point() +
    facet_wrap(.~lat.lon)

spring <- subset(d, as.numeric(doy)>60)

someyears <- subset(spring, year==1953|year==2003)
ggplot(someyears, aes(x=doy, y=tmean, group=lat.lon, colour=lat.lon)) +
    geom_point() +
    facet_wrap(.~cc)

someyears <- subset(spring, year<1970 | year>1990)
summary(lm(tmean~as.numeric(doy)*cc, data=someyears))

dagg <- aggregate(someyears["tmean"], someyears[c("cc", "year", "doy")], FUN=mean)

ggplot(dagg, aes(x=doy, y=tmean, group=year, colour=year)) +
    geom_point() +
    geom_smooth(method="lm") + 
    facet_wrap(.~cc)


if(FALSE){ # slow
library(rstanarm)
mod <- stan_lmer(tmean~(cc*as.numeric(doy))|lat.lon, data=spring)
}


## Set up the data to zero our <0 temps ...
mstmonths<-c(3:4)
dout <- d[which(as.numeric(d$mon) %in% mstmonths),]
dout$tempadj <- ifelse(dout$tmean<0, 0, dout$tmean)
doutagg <- aggregate(dout["tempadj"], dout[c("lo", "year", "lat.lon")], FUN=mean)
gdd <- aggregate(dout["tempadj"], dout[c("year", "lat.lon")], FUN=sum)
names(gdd)[names(gdd)=="tempadj"] <- "gdd"

datout <- merge(doutagg, gdd, by=c("year", "lat.lon"))
datout$site <- as.numeric(as.factor(datout$lat.lon))

write.csv(datout, "output/betpen_decsens_1950-2000base0.csv", row.names=FALSE)
