### Look at many species in PEP to choose best options
### 7 Sept 2020 - Cat
## based off of betpen_chillandgdd_tg.R but using Tmin and Tmax now to find Tmean

# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
require(dplyr)
require(tidyr)
require(lubridate)

use.10yr = TRUE

setwd("~/Documents/git/decsens/analyses/pep_analyses")
#d<-read.csv("input/pep_betpen_all.csv", header=TRUE)
#d<-read.csv("input/pep_fagsyl_all.csv", header=TRUE)
#d<-read.csv("input/pep_tilcor_all.csv", header=TRUE) ## 0 sites
#d<-read.csv("input/pep_tilpla_all.csv", header=TRUE) ## 0 sites
#d<-read.csv("input/pep_betpub_all.csv", header=TRUE)  ## 0 sites
#d<-read.csv("input/pep_alnglu_all.csv", header=TRUE)
#d<-read.csv("input/pep_querob_all.csv", header=TRUE)
#d<-read.csv("input/pep_fraexc_all.csv", header=TRUE)
#d<-read.csv("input/pep_samnig_all.csv", header=TRUE)
#d<-read.csv("input/pep_robpse_all.csv", header=TRUE)
#d<-read.csv("input/pep_poptre_all.csv", header=TRUE)
#d<-read.csv("input/pep_cormas_all.csv", header=TRUE)
#d<-read.csv("input/pep_lardec_all.csv", header=TRUE)

df<-d%>%
  filter(BBCH==11)%>%
  filter(YEAR>=1950 & YEAR<=2010)%>%
  dplyr::select(YEAR, DAY, BBCH, PEP_ID, LAT, LON, species)%>%
  rename(year=YEAR)%>%
  rename(lo=DAY)%>%
  rename(lat=LAT)%>%
  rename(long=LON)
## Hmm... can we sequence from budburst to leafout to find the number of freezes between?
df<-dplyr::select(df, year, PEP_ID, lat, long, lo)

df<-df[!duplicated(df),]

x<-paste(df$year, df$lo)
df$date<-as.Date(strptime(x, format="%Y %j"))
df$Date<- as.character(df$date)
df$lat.long <- paste(df$lat, df$long)
allpeps <- df[(df$year>=1951 & df$year<=2010),]
if(use.10yr==TRUE){
  allpeps <- allpeps[!(allpeps$year<=2000 & allpeps$year>=1961),] ### prep for 10 year windows
}

allpeps$cc<-NA
if(use.10yr==TRUE){
allpeps$cc<-ifelse(allpeps$year>1950 & allpeps$year<=1970, "1950-1960", allpeps$cc)
allpeps$cc<-ifelse(allpeps$year>1990 & allpeps$year<=2010, "2000-2010", allpeps$cc)
} else {

allpeps$cc<-ifelse(allpeps$year>1970 & allpeps$year<=1990, "1970-1990", allpeps$cc)
allpeps$cc<-ifelse(allpeps$year>1950 & allpeps$year<=1970, "1950-1970", allpeps$cc)
allpeps$cc<-ifelse(allpeps$year>1990 & allpeps$year<=2010, "1990-2010", allpeps$cc)
}


allpeps$num.years<-ave(allpeps$year, allpeps$PEP_ID, FUN=length)
mostdata <- allpeps[(allpeps$num.years>= if(use.10yr==TRUE){20} else {60}),]

tt<-as.data.frame(table(mostdata$cc, mostdata$lat.long))
tt<-tt[!(tt$Freq==0),]
bestsites<-as.data.frame(table(tt$Var2))
bestsites<-bestsites[(bestsites$Freq>1),]
bestsites <- bestsites$Var1

allpeps.subset<-mostdata[(mostdata$lat.long %in% bestsites),]

sites<-subset(allpeps.subset, select=c(lat, long, lat.long))
sites<-sites[!duplicated(sites$lat.long),]
nsites<-length(sites$lat.long)


## betpen: 20yr = 17  ; 10yr = 45  ### Betula pendula
## fagsyl: 20yr = 24 ; 10yr = 47  ### Fagus sylvatica
## alnglu: 20yr = 5 ; 10yr =  19     ### Alnus glutinosa
## querob: 20yr = 20 ; 10yr =  43     ### Quercus robur
## fraexc: 20yr = 4 ; 10yr =  30     ### Fraxinus excelsior
## cormas: 20yr = 0 ; 10yr =  0   ### Cornus mas
## tilcor: 20yr = 0  ; 10yr =  0   ### Tilia cordata


## poptre: 20yr = 0 ; 10yr = 0   ### Populus tremuloides
## lardec: 20yr = 0 ; 10yr =  0  ### Larix decidua
## tilpla: 20yr = 0  ; 10yr =  0  ### Tilia platyphyllos
## betpub: 20yr = 0 ; 10yr =  0  ### Betula pubescens
## robpse: 20yr = 0 ; 10yr =  0  ### Robinia pseudoacacia
## samnig: 20yr = 0 ; 10yr = 0   ### Sambucus nigra
