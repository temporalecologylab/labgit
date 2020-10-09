## Started 17 April 2020 ##
## By Lizzie, but not going well ##

## Following betpen_chillandgdd_tntx_forsims_odyssey.R some ... ##
## and my climate notes ... ##

## Get leafout date and daily climate (say, Jan 1 to end of Apr for now) for full year range of
# good BETPEN PEP sites ##

# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
require(dplyr)
require(tidyr)
require(ggplot2)
require(lubridate)
require(raster)

#setwd("~/Documents/git/projects/treegarden/decsens/analyses/pep_analyses/")
setwd("~/Documents/git/decsens/analyses/pep_analyses/")

d <- read.csv("input/pep_betpen_all.csv")

df<-d%>%
  filter(BBCH==11)%>%
  filter(YEAR>=1950 & YEAR<=2010)%>%
  dplyr::select(YEAR, DAY, BBCH, PEP_ID, LAT, LON, species)%>%
  rename(year=YEAR)%>%
  rename(lo=DAY)%>%
  rename(lat=LAT)%>%
  rename(lon=LON)
df <- dplyr::select(df, year, PEP_ID, lat, lon, lo)

df <- df[!duplicated(df),]

x <- paste(df$year, df$lo)
df$date <-as.Date(strptime(x, format="%Y %j"))
df$Date <- as.character(df$date)
df$lat.lon <- paste(df$lat, df$lon)
allpeps <- df[(df$year>=1951 & df$year<=2010),]
allpeps$num.years<-ave(allpeps$year, allpeps$PEP_ID, FUN=length)
mostdata <- allpeps[(allpeps$num.years>=60),]

lositeyear <- subset(mostdata, select=c("lo", "lat", "lon", "lat.lon", "year"))
getsitez <- lositeyear[!duplicated(lositeyear$lat.lon),] # bad, should add lon, but works for these sites - changed column to lat.lon and then address below
getsitez$x <- getsitez$lon
getsitez$y  <- getsitez$lat
Coords <- subset(getsitez, select=c(x, y))
Coords <- na.omit(Coords)
nsites <- length(getsitez$lat.lon)

#rg <- brick("/Volumes/OrangeFiend/climate/eobs/tg_ens_mean_0.25deg_reg_v20.0e.nc", sep="")
rg <- brick("~/Desktop/Big Data Files/tg_0.25deg_reg_v16.0.nc")

period <- c(1951:2010)

points <- SpatialPoints(Coords, proj4string = rg@crs)

yearsinclim <- as.numeric(format(as.Date(names(rg),format="X%Y.%m.%d"),"%Y"))
yearsinperiod <- which(yearsinclim%in%period)
climsub <- subset(rg,yearsinperiod) # takes a couple of minutes

## subset climate days
monthsinclim <- as.numeric(format(as.Date(names(climsub),format="X%Y.%m.%d"),"%m"))
dailytempmonths <- c(1:4)
monthsindailytemps <- which(monthsinclim%in%dailytempmonths)
dailytempssub <- subset(climsub,monthsindailytemps)


values <- raster::extract(dailytempssub,points) ## takes many minutes!

dclim <- cbind.data.frame(coordinates(points),values)

require(reshape2)
extclimdata <- melt(dclim, id.vars=c("x","y"))

extclimdata <- extclimdata %>%
  rename(long=x) %>%
  rename(lat=y) %>%
  rename(date=variable) %>%
  rename(Tmean=value)

dailyclimdata <- data.frame(lat=extclimdata$lat, lon=extclimdata$long, date=extclimdata$date, tmean=extclimdata$Tmean)

dailyclimdata$date <- substr(dailyclimdata$date, 2,11)
dailyclimdata$Date <- gsub("[.]", "-", dailyclimdata$date)

dailyclimdata$date <- NULL

dailyclimdata$year <- as.numeric(substr(dailyclimdata$Date, 0, 4))
dailyclimdata$lat.lon <- paste(dailyclimdata$lat, dailyclimdata$lon)

dailyclimdata <- dailyclimdata[!duplicated(dailyclimdata),]

dailyandlo <- left_join(dailyclimdata, lositeyear)

write.csv(dailyandlo, file="~/Documents/git/decsens/analyses/pep_analyses/output/betpen_dailytempsandlo_1950to2010.csv", row.names = FALSE)

