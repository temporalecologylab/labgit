### Prepare Quercus robur for comparisons
### 7 June 2019 - Cat
## based off of betpen_chillandgdd_tntx_forsims.R but using Tmin and Tmax now to find Tmean

# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
require(dplyr)
require(tidyr)
require(ggplot2)
require(lubridate)
require(chillR)
require(raster)

use.10yr = TRUE

setwd("~/Documents/git/decsens/analyses/pep_analyses")
d<-read.csv("input/pep_querob_all.csv", header=TRUE)

rn<-brick("~/Desktop/Big Data Files/tn_0.25deg_reg_v16.0.nc", sep="")
rx<-brick("~/Desktop/Big Data Files/tx_0.25deg_reg_v16.0.nc", sep="")

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
badsites<-c("54.5 11.1", "49.7667 11.55", "47.8 11.0167") 
sites<-sites[!(sites$lat.long%in%badsites),]
nsites<-length(sites$lat.long)
sites$x<-sites$long
sites$y<-sites$lat
sites$siteslist<-1:nsites
tmin<-rn
tmax<-rx

lositeyear <- subset(allpeps.subset, select=c("lo", "lat", "long", "lat.long", "year"))
lositeyear <- lositeyear[!duplicated(lositeyear),]
lositeyear <- left_join(lositeyear, sites)
lositeyear<-na.omit(lositeyear)

leaps<-seq(from=1952, to=2010, by=4)
if(use.10yr==TRUE){
period<-1951:1960
} else {period <- 1951:1970}
## set function - depending on the period you are using
extractclimpre<-function(tmin,period){
  #extractclimpost<-function(tmin,period){
  #extractclim90s<-function(tmin,period){
  
  ## define array to store results
  nyears<-length(period)
  climateyears<-array(NA,dim=c(nyears, 5, nsites))
  row.names(climateyears)<-period
  colnames(climateyears)<-c("Mean.Utah", "Mean.Port", "Mean.GDD", "Spring.Temp", "Site Num.")
  
  ## subset climate years
  yearsinclimmin<-as.numeric(format(as.Date(names(tmin),format="X%Y.%m.%d"),"%Y"))
  yearsinclimmax<-as.numeric(format(as.Date(names(tmax),format="X%Y.%m.%d"),"%Y"))
  yearsinperiodmin<-which(yearsinclimmin%in%period)
  yearsinperiodmax<-which(yearsinclimmax%in%period)
  climsubmin<-subset(tmin,yearsinperiodmin)
  climsubmax<-subset(tmax,yearsinperiodmax)
  
  ## subset climate days
  monthsinclimmin<-as.numeric(format(as.Date(names(climsubmin),format="X%Y.%m.%d"),"%m"))
  monthsinclimmax<-as.numeric(format(as.Date(names(climsubmax),format="X%Y.%m.%d"),"%m"))
  chillmonths<-c(9:12,1:3)
  monthsinchillmin<-which(monthsinclimmin%in%chillmonths)
  monthsinchillmax<-which(monthsinclimmax%in%chillmonths)
  chillsubmin<-subset(climsubmin,monthsinchillmin)
  chillsubmax<-subset(climsubmax,monthsinchillmax)
  
  warmmonths<-c(1:5)
  monthsinwarmmin<-which(monthsinclimmin%in%warmmonths)
  monthsinwarmmax<-which(monthsinclimmax%in%warmmonths)
  warmsubmin<-subset(climsubmin,monthsinwarmmin)
  warmsubmax<-subset(climsubmax,monthsinwarmmax)
  
  ## commence loop  
  for (i in 1:nsites){#i=2
    print(i)
    sitesi<-sites$siteslist[i]
    
    ## load shape
    if(sitesi==sites$siteslist[i])
      Coords<-data.frame(sites$x, sites$y)
    points.min <- SpatialPoints(Coords, proj4string = rn@crs)
    points.max <- SpatialPoints(Coords, proj4string = rx@crs)
    
    ## loop across years to extract each years averages
    # save first an array to store results
    yearlyresults<-array(NA,dim=c(length(period),5))
    
    for(j in period){#j=2001
      print(paste(i,j))
      
      # select year's layer
      chillyearsmin<-which(as.numeric(format(as.Date(
        names(chillsubmin),format="X%Y.%m.%d"),"%Y"))==j)
      chillyearsmax<-which(as.numeric(format(as.Date(
        names(chillsubmax),format="X%Y.%m.%d"),"%Y"))==j)
      
      yearschillmin<-subset(chillsubmin,chillyearsmin)
      yearschillmax<-subset(chillsubmax,chillyearsmax)
      
      # extract values and format to compute means
      tempschillsmin<-raster::extract(yearschillmin,points.min)
      tempschillsmax<-raster::extract(yearschillmax,points.max)
      
      #turn into data frame and remove NAs
      chmin<-as.data.frame(tempschillsmin)
      chmin<-subset(chmin,!is.na(rowSums(chmin)))
      chmax<-as.data.frame(tempschillsmax)
      chmax<-subset(chmax,!is.na(rowSums(chmax)))
      
      ## calculate chilling
      chillunitseachcelleachdaymin<-apply(chmin,2,function(x){
        Tmin<-x
        return(Tmin)})
      tminchill<-chillunitseachcelleachdaymin[(as.numeric(rownames(chillunitseachcelleachdaymin))==i)]
      
      chillunitseachcelleachdaymax<-apply(chmax,2,function(x){
        Tmax<-x
        return(Tmax)})
      tmaxchill<-chillunitseachcelleachdaymax[(as.numeric(rownames(chillunitseachcelleachdaymax))==i)]
      
      meandaily <- (tminchill + tmaxchill)/2
      
      x <- as.Date(substr(colnames(chillunitseachcelleachdaymin), 2, 11),format="%Y.%m.%d")
      
      hrly.temp=
        data.frame(
          Temp = c(rep(meandaily, each = 24)),
          Year = c(rep(as.numeric(substr(colnames(chillunitseachcelleachdaymin), 2, 5)), times=24)),
          #JDay = sort(c(rep(seq(1:length(colnames(meandaily))), times = 24)))
          JDay = sort(c(rep(yday(x), times=24)))
        )
      
      
      # select year's layer
      warmyearsmin<-which(as.numeric(format(as.Date(
        names(warmsubmin),format="X%Y.%m.%d"),"%Y"))==j)
      warmyearsmax<-which(as.numeric(format(as.Date(
        names(warmsubmax),format="X%Y.%m.%d"),"%Y"))==j)
      
      yearswarmmin<-subset(warmsubmin,warmyearsmin)
      yearswarmmax<-subset(warmsubmax,warmyearsmax)
      
      # extract values and format to compute means and sdevs
      tempswarmsmin<-raster::extract(yearswarmmin,points.min)
      tempswarmsmax<-raster::extract(yearswarmmax,points.max)
      
      #turn into data frame and remove NAs
      wamin<-as.data.frame(tempswarmsmin)
      wamin<-subset(wamin,!is.na(rowSums(wamin)))
      wamax<-as.data.frame(tempswarmsmax)
      wamax<-subset(wamax,!is.na(rowSums(wamax)))
      
      
      ## calculate forcing (GDD)
      warmunitseachcelleachdaymin<-apply(wamin,2,function(x){
        Tmin.warm<-x
        return(Tmin.warm)})
      tminwarm<-warmunitseachcelleachdaymin[(as.numeric(rownames(warmunitseachcelleachdaymin))==i)]
      
      warmunitseachcelleachdaymax<-apply(wamax,2,function(x){
        Tmax.warm<-x
        return(Tmax.warm)})
      tmaxwarm<-warmunitseachcelleachdaymax[(as.numeric(rownames(warmunitseachcelleachdaymax))==i)]
      
      meandaily.warm <- (tminwarm + tmaxwarm)/2
      
      
      x.warm <- as.Date(substr(colnames(warmunitseachcelleachdaymin), 2, 11),format="%Y.%m.%d")
      
      hrly.temp.warm =
        data.frame(
          Temp = c(rep(meandaily.warm, each = 24)),
          Year = c(rep(as.numeric(substr(colnames(warmunitseachcelleachdaymin), 2, 5)), times=24)),
          JDay = sort(c(rep(yday(x.warm), times = 24)))
        )
      
      lopersite <- lositeyear[(lositeyear$siteslist==i & lositeyear$year==j),]
      lo <- as.numeric(lopersite$lo)
      hrly.temp.warm <- hrly.temp.warm[(hrly.temp.warm$JDay<=lo & hrly.temp.warm$JDay>=lo-30),]
      
      chillcalc.mn<-chilling(hrly.temp, hrly.temp$JDay[1], hrly.temp$JDay[nrow(hrly.temp[1])])
      warmcalc.mn<-chilling(hrly.temp.warm, hrly.temp.warm$JDay[1], hrly.temp.warm$JDay[nrow(hrly.temp.warm[1])])
      
      
      yearlyresults[which(period==j),1]<-chillcalc.mn$Utah_Model[which(chillcalc.mn$End_year==j)]
      yearlyresults[which(period==j),2]<-chillcalc.mn$Chill_portions[which(chillcalc.mn$End_year==j)]
      yearlyresults[which(period==j),3]<-(warmcalc.mn$GDH[which(warmcalc.mn$End_year==j)])/24
      yearlyresults[which(period==j),4]<-mean(hrly.temp.warm$Temp, na.rm=TRUE)
      
      yearlyresults[which(period==j),5]<-sites$siteslist[i]
      
    }
    
    climateyears[,,i]<-yearlyresults
    
  } 
  
  return(climateyears)
  
}


## apply function
clim_pre<-extractclimpre(tmin,period)

pre<-as.data.frame(clim_pre)
if(use.10yr==TRUE){
write.csv(pre, file="~/Desktop/Misc/Ospree Misc/prequerobten.csv", row.names=FALSE)
} else {
  write.csv(pre, file="~/Desktop/Misc/Ospree Misc/prequerob.csv", row.names=FALSE)
}


if(use.10yr==TRUE){
  period<-2001:2010
} else {period<-1991:2010}
extractclimpost<-function(tmin,period){
  #extractclim90s<-function(tmin,period){
  
  ## define array to store results
  nyears<-length(period)
  climateyears<-array(NA,dim=c(nyears, 5, nsites))
  row.names(climateyears)<-period
  colnames(climateyears)<-c("Mean.Utah", "Mean.Port", "Mean.GDD", "Spring.Temp", "Site Num.")
  
  ## subset climate years
  yearsinclimmin<-as.numeric(format(as.Date(names(tmin),format="X%Y.%m.%d"),"%Y"))
  yearsinclimmax<-as.numeric(format(as.Date(names(tmax),format="X%Y.%m.%d"),"%Y"))
  yearsinperiodmin<-which(yearsinclimmin%in%period)
  yearsinperiodmax<-which(yearsinclimmax%in%period)
  climsubmin<-subset(tmin,yearsinperiodmin)
  climsubmax<-subset(tmax,yearsinperiodmax)
  
  ## subset climate days
  monthsinclimmin<-as.numeric(format(as.Date(names(climsubmin),format="X%Y.%m.%d"),"%m"))
  monthsinclimmax<-as.numeric(format(as.Date(names(climsubmax),format="X%Y.%m.%d"),"%m"))
  chillmonths<-c(9:12,1:3)
  monthsinchillmin<-which(monthsinclimmin%in%chillmonths)
  monthsinchillmax<-which(monthsinclimmax%in%chillmonths)
  chillsubmin<-subset(climsubmin,monthsinchillmin)
  chillsubmax<-subset(climsubmax,monthsinchillmax)
  
  warmmonths<-c(1:5)
  monthsinwarmmin<-which(monthsinclimmin%in%warmmonths)
  monthsinwarmmax<-which(monthsinclimmax%in%warmmonths)
  warmsubmin<-subset(climsubmin,monthsinwarmmin)
  warmsubmax<-subset(climsubmax,monthsinwarmmax)
  
  ## commence loop  
  for (i in 1:nsites){#i=2
    print(i)
    sitesi<-sites$siteslist[i]
    
    ## load shape
    if(sitesi==sites$siteslist[i])
      Coords<-data.frame(sites$x, sites$y)
    points.min <- SpatialPoints(Coords, proj4string = rn@crs)
    points.max <- SpatialPoints(Coords, proj4string = rx@crs)
    
    ## loop across years to extract each years averages
    # save first an array to store results
    yearlyresults<-array(NA,dim=c(length(period),5))
    
    for(j in period){#j=2001
      print(paste(i,j))
      
      # select year's layer
      chillyearsmin<-which(as.numeric(format(as.Date(
        names(chillsubmin),format="X%Y.%m.%d"),"%Y"))==j)
      chillyearsmax<-which(as.numeric(format(as.Date(
        names(chillsubmax),format="X%Y.%m.%d"),"%Y"))==j)
      
      yearschillmin<-subset(chillsubmin,chillyearsmin)
      yearschillmax<-subset(chillsubmax,chillyearsmax)
      
      # extract values and format to compute means
      tempschillsmin<-raster::extract(yearschillmin,points.min)
      tempschillsmax<-raster::extract(yearschillmax,points.max)
      
      #turn into data frame and remove NAs
      chmin<-as.data.frame(tempschillsmin)
      chmin<-subset(chmin,!is.na(rowSums(chmin)))
      chmax<-as.data.frame(tempschillsmax)
      chmax<-subset(chmax,!is.na(rowSums(chmax)))
      
      ## calculate chilling
      chillunitseachcelleachdaymin<-apply(chmin,2,function(x){
        Tmin<-x
        return(Tmin)})
      tminchill<-chillunitseachcelleachdaymin[(as.numeric(rownames(chillunitseachcelleachdaymin))==i)]
      
      chillunitseachcelleachdaymax<-apply(chmax,2,function(x){
        Tmax<-x
        return(Tmax)})
      tmaxchill<-chillunitseachcelleachdaymax[(as.numeric(rownames(chillunitseachcelleachdaymax))==i)]
      
      meandaily <- (tminchill + tmaxchill)/2
      
      x <- as.Date(substr(colnames(chillunitseachcelleachdaymin), 2, 11),format="%Y.%m.%d")
      
      hrly.temp=
        data.frame(
          Temp = c(rep(meandaily, each = 24)),
          Year = c(rep(as.numeric(substr(colnames(chillunitseachcelleachdaymin), 2, 5)), times=24)),
          #JDay = sort(c(rep(seq(1:length(colnames(meandaily))), times = 24)))
          JDay = sort(c(rep(yday(x), times=24)))
        )
      
      
      # select year's layer
      warmyearsmin<-which(as.numeric(format(as.Date(
        names(warmsubmin),format="X%Y.%m.%d"),"%Y"))==j)
      warmyearsmax<-which(as.numeric(format(as.Date(
        names(warmsubmax),format="X%Y.%m.%d"),"%Y"))==j)
      
      yearswarmmin<-subset(warmsubmin,warmyearsmin)
      yearswarmmax<-subset(warmsubmax,warmyearsmax)
      
      # extract values and format to compute means and sdevs
      tempswarmsmin<-raster::extract(yearswarmmin,points.min)
      tempswarmsmax<-raster::extract(yearswarmmax,points.max)
      
      #turn into data frame and remove NAs
      wamin<-as.data.frame(tempswarmsmin)
      wamin<-subset(wamin,!is.na(rowSums(wamin)))
      wamax<-as.data.frame(tempswarmsmax)
      wamax<-subset(wamax,!is.na(rowSums(wamax)))
      
      
      ## calculate forcing (GDD)
      warmunitseachcelleachdaymin<-apply(wamin,2,function(x){
        Tmin.warm<-x
        return(Tmin.warm)})
      tminwarm<-warmunitseachcelleachdaymin[(as.numeric(rownames(warmunitseachcelleachdaymin))==i)]
      
      warmunitseachcelleachdaymax<-apply(wamax,2,function(x){
        Tmax.warm<-x
        return(Tmax.warm)})
      tmaxwarm<-warmunitseachcelleachdaymax[(as.numeric(rownames(warmunitseachcelleachdaymax))==i)]
      
      meandaily.warm <- (tminwarm + tmaxwarm)/2
      
      
      x.warm <- as.Date(substr(colnames(warmunitseachcelleachdaymin), 2, 11),format="%Y.%m.%d")
      
      hrly.temp.warm =
        data.frame(
          Temp = c(rep(meandaily.warm, each = 24)),
          Year = c(rep(as.numeric(substr(colnames(warmunitseachcelleachdaymin), 2, 5)), times=24)),
          JDay = sort(c(rep(yday(x.warm), times = 24)))
        )
      
      lopersite <- lositeyear[(lositeyear$siteslist==i & lositeyear$year==j),]
      lo <- as.numeric(lopersite$lo)
      hrly.temp.warm <- hrly.temp.warm[(hrly.temp.warm$JDay<=lo & hrly.temp.warm$JDay>=lo-30),]
      
      chillcalc.mn<-chilling(hrly.temp, hrly.temp$JDay[1], hrly.temp$JDay[nrow(hrly.temp[1])])
      warmcalc.mn<-chilling(hrly.temp.warm, hrly.temp.warm$JDay[1], hrly.temp.warm$JDay[nrow(hrly.temp.warm[1])])
      
      
      yearlyresults[which(period==j),1]<-chillcalc.mn$Utah_Model[which(chillcalc.mn$End_year==j)]
      yearlyresults[which(period==j),2]<-chillcalc.mn$Chill_portions[which(chillcalc.mn$End_year==j)]
      yearlyresults[which(period==j),3]<-(warmcalc.mn$GDH[which(warmcalc.mn$End_year==j)])/24
      yearlyresults[which(period==j),4]<-mean(hrly.temp.warm$Temp, na.rm=TRUE)
      
      yearlyresults[which(period==j),5]<-sites$siteslist[i]
      
    }
    
    climateyears[,,i]<-yearlyresults
    
  } 
  
  return(climateyears)
  
}

clim_post<-extractclimpost(tmin,period) 
post<-as.data.frame(clim_post)
if(use.10yr==TRUE){
  write.csv(post, file="~/Desktop/Misc/Ospree Misc/postquerobten.csv", row.names=FALSE)
} else {
  write.csv(post, file="~/Desktop/Misc/Ospree Misc/postquerob.csv", row.names=FALSE)
}

if(use.10yr==FALSE){
  period<-1971:1990
  extractclimmid<-function(tmin,period){
    #extractclim90s<-function(tmin,period){
    
    ## define array to store results
    nyears<-length(period)
    climateyears<-array(NA,dim=c(nyears, 5, nsites))
    row.names(climateyears)<-period
    colnames(climateyears)<-c("Mean.Utah", "Mean.Port", "Mean.GDD", "Spring.Temp", "Site Num.")
    
    ## subset climate years
    yearsinclimmin<-as.numeric(format(as.Date(names(tmin),format="X%Y.%m.%d"),"%Y"))
    yearsinclimmax<-as.numeric(format(as.Date(names(tmax),format="X%Y.%m.%d"),"%Y"))
    yearsinperiodmin<-which(yearsinclimmin%in%period)
    yearsinperiodmax<-which(yearsinclimmax%in%period)
    climsubmin<-subset(tmin,yearsinperiodmin)
    climsubmax<-subset(tmax,yearsinperiodmax)
    
    ## subset climate days
    monthsinclimmin<-as.numeric(format(as.Date(names(climsubmin),format="X%Y.%m.%d"),"%m"))
    monthsinclimmax<-as.numeric(format(as.Date(names(climsubmax),format="X%Y.%m.%d"),"%m"))
    chillmonths<-c(9:12,1:3)
    monthsinchillmin<-which(monthsinclimmin%in%chillmonths)
    monthsinchillmax<-which(monthsinclimmax%in%chillmonths)
    chillsubmin<-subset(climsubmin,monthsinchillmin)
    chillsubmax<-subset(climsubmax,monthsinchillmax)
    
    warmmonths<-c(1:5)
    monthsinwarmmin<-which(monthsinclimmin%in%warmmonths)
    monthsinwarmmax<-which(monthsinclimmax%in%warmmonths)
    warmsubmin<-subset(climsubmin,monthsinwarmmin)
    warmsubmax<-subset(climsubmax,monthsinwarmmax)
    
    ## commence loop  
    for (i in 1:nsites){#i=1
      print(i)
      sitesi<-sites$siteslist[i]
      
      ## load shape
      if(sitesi==sites$siteslist[i])
        Coords<-data.frame(sites$x, sites$y)
      points.min <- SpatialPoints(Coords, proj4string = rn@crs)
      points.max <- SpatialPoints(Coords, proj4string = rx@crs)
      
      ## loop across years to extract each years averages
      # save first an array to store results
      yearlyresults<-array(NA,dim=c(length(period),5))
      
      for(j in period){#j=1951
        print(paste(i,j))
        
        # select year's layer
        chillyearsmin<-which(as.numeric(format(as.Date(
          names(chillsubmin),format="X%Y.%m.%d"),"%Y"))==j)
        chillyearsmax<-which(as.numeric(format(as.Date(
          names(chillsubmax),format="X%Y.%m.%d"),"%Y"))==j)
        
        yearschillmin<-subset(chillsubmin,chillyearsmin)
        yearschillmax<-subset(chillsubmax,chillyearsmax)
        
        # extract values and format to compute means
        tempschillsmin<-raster::extract(yearschillmin,points.min)
        tempschillsmax<-raster::extract(yearschillmax,points.max)
        
        #turn into data frame and remove NAs
        chmin<-as.data.frame(tempschillsmin)
        chmin<-subset(chmin,!is.na(rowSums(chmin)))
        chmax<-as.data.frame(tempschillsmax)
        chmax<-subset(chmax,!is.na(rowSums(chmax)))
        
        ## calculate chilling
        chillunitseachcelleachdaymin<-apply(chmin,2,function(x){
          Tmin<-x
          return(Tmin)})
        tminchill<-chillunitseachcelleachdaymin[(as.numeric(rownames(chillunitseachcelleachdaymin))==i)]
        
        chillunitseachcelleachdaymax<-apply(chmax,2,function(x){
          Tmax<-x
          return(Tmax)})
        tmaxchill<-chillunitseachcelleachdaymax[(as.numeric(rownames(chillunitseachcelleachdaymax))==i)]
        
        meandaily <- (tminchill + tmaxchill)/2
        
        x <- as.Date(substr(colnames(chillunitseachcelleachdaymin), 2, 11),format="%Y.%m.%d")
        
        hrly.temp=
          data.frame(
            Temp = c(rep(meandaily, each = 24)),
            Year = c(rep(as.numeric(substr(colnames(chillunitseachcelleachdaymin), 2, 5)), times=24)),
            #JDay = sort(c(rep(seq(1:length(colnames(meandaily))), times = 24)))
            JDay = sort(c(rep(yday(x), times=24)))
          )
        
        
        # select year's layer
        warmyearsmin<-which(as.numeric(format(as.Date(
          names(warmsubmin),format="X%Y.%m.%d"),"%Y"))==j)
        warmyearsmax<-which(as.numeric(format(as.Date(
          names(warmsubmax),format="X%Y.%m.%d"),"%Y"))==j)
        
        yearswarmmin<-subset(warmsubmin,warmyearsmin)
        yearswarmmax<-subset(warmsubmax,warmyearsmax)
        
        # extract values and format to compute means and sdevs
        tempswarmsmin<-raster::extract(yearswarmmin,points.min)
        tempswarmsmax<-raster::extract(yearswarmmax,points.max)
        
        #turn into data frame and remove NAs
        wamin<-as.data.frame(tempswarmsmin)
        wamin<-subset(wamin,!is.na(rowSums(wamin)))
        wamax<-as.data.frame(tempswarmsmax)
        wamax<-subset(wamax,!is.na(rowSums(wamax)))
        
        
        ## calculate forcing (GDD)
        warmunitseachcelleachdaymin<-apply(wamin,2,function(x){
          Tmin.warm<-x
          return(Tmin.warm)})
        tminwarm<-warmunitseachcelleachdaymin[(as.numeric(rownames(warmunitseachcelleachdaymin))==i)]
        
        warmunitseachcelleachdaymax<-apply(wamax,2,function(x){
          Tmax.warm<-x
          return(Tmax.warm)})
        tmaxwarm<-warmunitseachcelleachdaymax[(as.numeric(rownames(warmunitseachcelleachdaymax))==i)]
        
        meandaily.warm <- (tminwarm + tmaxwarm)/2
        
        
        x.warm <- as.Date(substr(colnames(warmunitseachcelleachdaymin), 2, 11),format="%Y.%m.%d")
        
        hrly.temp.warm =
          data.frame(
            Temp = c(rep(meandaily.warm, each = 24)),
            Year = c(rep(as.numeric(substr(colnames(warmunitseachcelleachdaymin), 2, 5)), times=24)),
            JDay = sort(c(rep(yday(x.warm), times = 24)))
          )
        
        lopersite <- lositeyear[(lositeyear$siteslist==i & lositeyear$year==j),]
        lo <- as.numeric(lopersite$lo)
        hrly.temp.warm <- hrly.temp.warm[(hrly.temp.warm$JDay<=lo & hrly.temp.warm$JDay>=lo-30),]
        
        chillcalc.mn<-chilling(hrly.temp, hrly.temp$JDay[1], hrly.temp$JDay[nrow(hrly.temp[1])])
        warmcalc.mn<-chilling(hrly.temp.warm, hrly.temp.warm$JDay[1], hrly.temp.warm$JDay[nrow(hrly.temp.warm[1])])
        
        
        yearlyresults[which(period==j),1]<-chillcalc.mn$Utah_Model[which(chillcalc.mn$End_year==j)]
        yearlyresults[which(period==j),2]<-chillcalc.mn$Chill_portions[which(chillcalc.mn$End_year==j)]
        yearlyresults[which(period==j),3]<-(warmcalc.mn$GDH[which(warmcalc.mn$End_year==j)])/24
        yearlyresults[which(period==j),4]<-mean(hrly.temp.warm$Temp, na.rm=TRUE)
        
        yearlyresults[which(period==j),5]<-sites$siteslist[i]
        
      }
      
      climateyears[,,i]<-yearlyresults
      
    } 
    
    return(climateyears)
    
  }
  
  clim_mid<-extractclimmid(tmin,period) 
  mid<-as.data.frame(clim_mid)
  write.csv(mid, file="~/Desktop/Misc/Ospree misc/midquerob.csv", row.names=FALSE)
}


#setwd("~/Documents/git/decsens/analyses/pep_analyses/output/zarchive")

savepre <- pre
savepost <- post
savemid <- mid

if(use.10yr==TRUE){
pre <- read.csv("~/Desktop/Misc/Ospree Misc/prequerobten.csv")
post <- read.csv("~/Desktop/Misc/Ospree Misc/postquerobten.csv")
} else{
  pre <- read.csv("~/Desktop/Misc/Ospree Misc/prequerob.csv")
  post <- read.csv("~/Desktop/Misc/Ospree Misc/postquerob.csv")
  mid <- read.csv("~/Desktop/Misc/Ospree Misc/midquerob.csv")
}

predata<-data.frame(chillutah = c(pre$Mean.Utah.1, pre$Mean.Utah.2,
                                  pre$Mean.Utah.3, pre$Mean.Utah.4,
                                  pre$Mean.Utah.5, pre$Mean.Utah.6,
                                  pre$Mean.Utah.7, pre$Mean.Utah.8,
                                  pre$Mean.Utah.9, pre$Mean.Utah.10,
                                  pre$Mean.Utah.11, pre$Mean.Utah.12,
                                  pre$Mean.Utah.13, pre$Mean.Utah.14,
                                  pre$Mean.Utah.15, pre$Mean.Utah.16,
                                  pre$Mean.Utah.17, pre$Mean.Utah.18,
                                  pre$Mean.Utah.19, pre$Mean.Utah.20,
                                  pre$Mean.Utah.21, pre$Mean.Utah.22,
                                  pre$Mean.Utah.23, pre$Mean.Utah.24,
                                  pre$Mean.Utah.25, pre$Mean.Utah.26,
                                  pre$Mean.Utah.27, pre$Mean.Utah.28,
                                  pre$Mean.Utah.29, pre$Mean.Utah.30,
                                  pre$Mean.Utah.31, pre$Mean.Utah.32,
                                  pre$Mean.Utah.33, pre$Mean.Utah.34,
                                  pre$Mean.Utah.35, pre$Mean.Utah.36,
                                  pre$Mean.Utah.37, pre$Mean.Utah.38,
                                  pre$Mean.Utah.39, pre$Mean.Utah.40,
                                  pre$Mean.Utah.41, pre$Mean.Utah.42,
                                  pre$Mean.Utah.43),
                    
                    chillports = c(pre$Mean.Port.1, pre$Mean.Port.2,
                                   pre$Mean.Port.3, pre$Mean.Port.4,
                                   pre$Mean.Port.5, pre$Mean.Port.6,
                                   pre$Mean.Port.7, pre$Mean.Port.8,
                                   pre$Mean.Port.9, pre$Mean.Port.10,
                                   pre$Mean.Port.11, pre$Mean.Port.12,
                                   pre$Mean.Port.13, pre$Mean.Port.14,
                                   pre$Mean.Port.15, pre$Mean.Port.16,
                                   pre$Mean.Port.17, pre$Mean.Port.18,
                                   pre$Mean.Port.19, pre$Mean.Port.20,
                                   pre$Mean.Port.21, pre$Mean.Port.22,
                                   pre$Mean.Port.23, pre$Mean.Port.24,
                                   pre$Mean.Port.25, pre$Mean.Port.26,
                                   pre$Mean.Port.27, pre$Mean.Port.28,
                                   pre$Mean.Port.29, pre$Mean.Port.30,
                                   pre$Mean.Port.31, pre$Mean.Port.32,
                                   pre$Mean.Port.33, pre$Mean.Port.34,
                                   pre$Mean.Port.35, pre$Mean.Port.36,
                                   pre$Mean.Port.37, pre$Mean.Port.38,
                                   pre$Mean.Port.39, pre$Mean.Port.40,
                                   pre$Mean.Port.41, pre$Mean.Port.42,
                                   pre$Mean.Port.43),
                    
                    
                    gdd = c(pre$Mean.GDD.1, pre$Mean.GDD.2,
                            pre$Mean.GDD.3, pre$Mean.GDD.4,
                            pre$Mean.GDD.5, pre$Mean.GDD.6,
                            pre$Mean.GDD.7, pre$Mean.GDD.8,
                            pre$Mean.GDD.9, pre$Mean.GDD.10,
                            pre$Mean.GDD.11, pre$Mean.GDD.12,
                            pre$Mean.GDD.13, pre$Mean.GDD.14,
                            pre$Mean.GDD.15, pre$Mean.GDD.16,
                            pre$Mean.GDD.17, pre$Mean.GDD.18,
                            pre$Mean.GDD.19, pre$Mean.GDD.20,
                            pre$Mean.GDD.21, pre$Mean.GDD.22,
                            pre$Mean.GDD.23, pre$Mean.GDD.24,
                            pre$Mean.GDD.25, pre$Mean.GDD.26,
                            pre$Mean.GDD.27, pre$Mean.GDD.28,
                            pre$Mean.GDD.29, pre$Mean.GDD.30,
                            pre$Mean.GDD.31, pre$Mean.GDD.32,
                            pre$Mean.GDD.33, pre$Mean.GDD.34,
                            pre$Mean.GDD.35, pre$Mean.GDD.36,
                            pre$Mean.GDD.37, pre$Mean.GDD.38,
                            pre$Mean.GDD.39, pre$Mean.GDD.40,
                            pre$Mean.GDD.41, pre$Mean.GDD.42,
                            pre$Mean.GDD.43),
                    
                    mat.lo = c(pre$Spring.Temp.1, pre$Spring.Temp.2,
                               pre$Spring.Temp.3, pre$Spring.Temp.4,
                               pre$Spring.Temp.5, pre$Spring.Temp.6,
                               pre$Spring.Temp.7, pre$Spring.Temp.8,
                               pre$Spring.Temp.9, pre$Spring.Temp.10,
                               pre$Spring.Temp.11, pre$Spring.Temp.12,
                               pre$Spring.Temp.13, pre$Spring.Temp.14,
                               pre$Spring.Temp.15, pre$Spring.Temp.16,
                               pre$Spring.Temp.17, pre$Spring.Temp.18,
                               pre$Spring.Temp.19, pre$Spring.Temp.20,
                               pre$Spring.Temp.21, pre$Spring.Temp.22,
                               pre$Spring.Temp.23, pre$Spring.Temp.24,
                               pre$Spring.Temp.25, pre$Spring.Temp.26,
                               pre$Spring.Temp.27, pre$Spring.Temp.28,
                               pre$Spring.Temp.29, pre$Spring.Temp.30,
                               pre$Spring.Temp.31, pre$Spring.Temp.32,
                               pre$Spring.Temp.33, pre$Spring.Temp.34,
                               pre$Spring.Temp.35, pre$Spring.Temp.36,
                               pre$Spring.Temp.37, pre$Spring.Temp.38,
                               pre$Spring.Temp.39, pre$Spring.Temp.40,
                               pre$Spring.Temp.41, pre$Spring.Temp.42,
                               pre$Spring.Temp.43),
                    
                    siteslist = c(pre$Site.Num..1, pre$Site.Num..2,
                                  pre$Site.Num..3, pre$Site.Num..4,
                                  pre$Site.Num..5, pre$Site.Num..6,
                                  pre$Site.Num..7, pre$Site.Num..8,
                                  pre$Site.Num..9, pre$Site.Num..10,
                                  pre$Site.Num..11, pre$Site.Num..12,
                                  pre$Site.Num..13, pre$Site.Num..14,
                                  pre$Site.Num..15, pre$Site.Num..16,
                                  pre$Site.Num..17, pre$Site.Num..18,
                                  pre$Site.Num..19, pre$Site.Num..20,
                                  pre$Site.Num..21, pre$Site.Num..22,
                                  pre$Site.Num..23, pre$Site.Num..24,
                                  pre$Site.Num..25, pre$Site.Num..26,
                                  pre$Site.Num..27, pre$Site.Num..28,
                                  pre$Site.Num..29, pre$Site.Num..30,
                                  pre$Site.Num..31, pre$Site.Num..32,
                                  pre$Site.Num..33, pre$Site.Num..34,
                                  pre$Site.Num..35, pre$Site.Num..36,
                                  pre$Site.Num..37, pre$Site.Num..38,
                                  pre$Site.Num..39, pre$Site.Num..40,
                                  pre$Site.Num..41, pre$Site.Num..42,
                                  pre$Site.Num..43),
                    year = (as.numeric(rownames(pre))+1950))

site<-full_join(predata, sites)
site$x<-NULL
site$y<-NULL  

postdata<-data.frame(chillutah = c(post$Mean.Utah.1, post$Mean.Utah.2,
                                   post$Mean.Utah.3, post$Mean.Utah.4,
                                   post$Mean.Utah.5, post$Mean.Utah.6,
                                   post$Mean.Utah.7, post$Mean.Utah.8,
                                   post$Mean.Utah.9, post$Mean.Utah.10,
                                   post$Mean.Utah.11, post$Mean.Utah.12,
                                   post$Mean.Utah.13, post$Mean.Utah.14,
                                   post$Mean.Utah.15, post$Mean.Utah.16,
                                   post$Mean.Utah.17, post$Mean.Utah.18,
                                   post$Mean.Utah.19, post$Mean.Utah.20,
                                   post$Mean.Utah.21, post$Mean.Utah.22,
                                   post$Mean.Utah.23, post$Mean.Utah.24,
                                   post$Mean.Utah.25, post$Mean.Utah.26,
                                   post$Mean.Utah.27, post$Mean.Utah.28,
                                   post$Mean.Utah.29, post$Mean.Utah.30,
                                   post$Mean.Utah.31, post$Mean.Utah.32,
                                   post$Mean.Utah.33, post$Mean.Utah.34,
                                   post$Mean.Utah.35, post$Mean.Utah.36,
                                   post$Mean.Utah.37, post$Mean.Utah.38,
                                   post$Mean.Utah.39, post$Mean.Utah.40,
                                   post$Mean.Utah.41, post$Mean.Utah.42,
                                   post$Mean.Utah.43, post$Mean.Utah.44,
                                   post$Mean.Utah.45, post$Mean.Utah.46,
                                   post$Mean.Utah.47),
                     
                     chillports = c(post$Mean.Port.1, post$Mean.Port.2,
                                    post$Mean.Port.3, post$Mean.Port.4,
                                    post$Mean.Port.5, post$Mean.Port.6,
                                    post$Mean.Port.7, post$Mean.Port.8,
                                    post$Mean.Port.9, post$Mean.Port.10,
                                    post$Mean.Port.11, post$Mean.Port.12,
                                    post$Mean.Port.13, post$Mean.Port.14,
                                    post$Mean.Port.15, post$Mean.Port.16,
                                    post$Mean.Port.17, post$Mean.Port.18,
                                    post$Mean.Port.19, post$Mean.Port.20,
                                    post$Mean.Port.21, post$Mean.Port.22,
                                    post$Mean.Port.23, post$Mean.Port.24,
                                    post$Mean.Port.25, post$Mean.Port.26,
                                    post$Mean.Port.27, post$Mean.Port.28,
                                    post$Mean.Port.29, post$Mean.Port.30,
                                    post$Mean.Port.31, post$Mean.Port.32,
                                    post$Mean.Port.33, post$Mean.Port.34,
                                    post$Mean.Port.35, post$Mean.Port.36,
                                    post$Mean.Port.37, post$Mean.Port.38,
                                    post$Mean.Port.39, post$Mean.Port.40,
                                    post$Mean.Port.41, post$Mean.Port.42,
                                    post$Mean.Port.43, post$Mean.Port.44,
                                    post$Mean.Port.45, post$Mean.Port.46,
                                    post$Mean.Port.47),
                     
                     
                     gdd = c(post$Mean.GDD.1, post$Mean.GDD.2,
                             post$Mean.GDD.3, post$Mean.GDD.4,
                             post$Mean.GDD.5, post$Mean.GDD.6,
                             post$Mean.GDD.7, post$Mean.GDD.8,
                             post$Mean.GDD.9, post$Mean.GDD.10,
                             post$Mean.GDD.11, post$Mean.GDD.12,
                             post$Mean.GDD.13, post$Mean.GDD.14,
                             post$Mean.GDD.15, post$Mean.GDD.16,
                             post$Mean.GDD.17, post$Mean.GDD.18,
                             post$Mean.GDD.19, post$Mean.GDD.20,
                             post$Mean.GDD.21, post$Mean.GDD.22,
                             post$Mean.GDD.23, post$Mean.GDD.24,
                             post$Mean.GDD.25, post$Mean.GDD.26,
                             post$Mean.GDD.27, post$Mean.GDD.28,
                             post$Mean.GDD.29, post$Mean.GDD.30,
                             post$Mean.GDD.31, post$Mean.GDD.32,
                             post$Mean.GDD.33, post$Mean.GDD.34,
                             post$Mean.GDD.35, post$Mean.GDD.36,
                             post$Mean.GDD.37, post$Mean.GDD.38,
                             post$Mean.GDD.39, post$Mean.GDD.40,
                             post$Mean.GDD.41, post$Mean.GDD.42,
                             post$Mean.GDD.43, post$Mean.GDD.44,
                             post$Mean.GDD.45, post$Mean.GDD.46,
                             post$Mean.GDD.47),
                     
                     mat.lo = c(post$Spring.Temp.1, post$Spring.Temp.2,
                                post$Spring.Temp.3, post$Spring.Temp.4,
                                post$Spring.Temp.5, post$Spring.Temp.6,
                                post$Spring.Temp.7, post$Spring.Temp.8,
                                post$Spring.Temp.9, post$Spring.Temp.10,
                                post$Spring.Temp.11, post$Spring.Temp.12,
                                post$Spring.Temp.13, post$Spring.Temp.14,
                                post$Spring.Temp.15, post$Spring.Temp.16,
                                post$Spring.Temp.17, post$Spring.Temp.18,
                                post$Spring.Temp.19, post$Spring.Temp.20,
                                post$Spring.Temp.21, post$Spring.Temp.22,
                                post$Spring.Temp.23, post$Spring.Temp.24,
                                post$Spring.Temp.25, post$Spring.Temp.26,
                                post$Spring.Temp.27, post$Spring.Temp.28,
                                post$Spring.Temp.29, post$Spring.Temp.30,
                                post$Spring.Temp.31, post$Spring.Temp.32,
                                post$Spring.Temp.33, post$Spring.Temp.34,
                                post$Spring.Temp.35, post$Spring.Temp.36,
                                post$Spring.Temp.37, post$Spring.Temp.38,
                                post$Spring.Temp.39, post$Spring.Temp.40,
                                post$Spring.Temp.41, post$Spring.Temp.42,
                                post$Spring.Temp.43, post$Spring.Temp.44,
                                post$Spring.Temp.45, post$Spring.Temp.46, 
                                post$Spring.Temp.47),
                     
                     siteslist = c(post$Site.Num..1, post$Site.Num..2,
                                   post$Site.Num..3, post$Site.Num..4,
                                   post$Site.Num..5, post$Site.Num..6,
                                   post$Site.Num..7, post$Site.Num..8,
                                   post$Site.Num..9, post$Site.Num..10,
                                   post$Site.Num..11, post$Site.Num..12,
                                   post$Site.Num..13, post$Site.Num..14,
                                   post$Site.Num..15, post$Site.Num..16,
                                   post$Site.Num..17, post$Site.Num..18,
                                   post$Site.Num..19, post$Site.Num..20,
                                   post$Site.Num..21, post$Site.Num..22,
                                   post$Site.Num..23, post$Site.Num..24,
                                   post$Site.Num..25, post$Site.Num..26,
                                   post$Site.Num..27, post$Site.Num..28,
                                   post$Site.Num..29, post$Site.Num..30,
                                   post$Site.Num..31, post$Site.Num..32,
                                   post$Site.Num..33, post$Site.Num..34,
                                   post$Site.Num..35, post$Site.Num..36,
                                   post$Site.Num..37, post$Site.Num..38,
                                   post$Site.Num..39, post$Site.Num..40,
                                   post$Site.Num..41, post$Site.Num..42,
                                   post$Site.Num..43, post$Site.Num..44,
                                   post$Site.Num..45, post$Site.Num..46,
                                   post$Site.Num..47),
                     year = (as.numeric(rownames(post))+if(use.10yr==TRUE){2000}else{1990}))

site.post<-full_join(postdata, sites)
site.post$x<-NULL
site.post$y<-NULL

if(use.10yr==FALSE){
  middata<-data.frame(chillutah = c(mid$Mean.Utah.1, mid$Mean.Utah.2,
                                    mid$Mean.Utah.3, mid$Mean.Utah.4,
                                    mid$Mean.Utah.5, mid$Mean.Utah.6,
                                    mid$Mean.Utah.7, mid$Mean.Utah.8,
                                    mid$Mean.Utah.9, mid$Mean.Utah.10,
                                    mid$Mean.Utah.11, mid$Mean.Utah.12,
                                    mid$Mean.Utah.13, mid$Mean.Utah.14,
                                    mid$Mean.Utah.15, mid$Mean.Utah.16,
                                    mid$Mean.Utah.17, mid$Mean.Utah.18,
                                    mid$Mean.Utah.19, mid$Mean.Utah.20,
                                    mid$Mean.Utah.21, mid$Mean.Utah.22,
                                    mid$Mean.Utah.23, mid$Mean.Utah.24,
                                    mid$Mean.Utah.25, mid$Mean.Utah.26,
                                    mid$Mean.Utah.27, mid$Mean.Utah.28,
                                    mid$Mean.Utah.29, mid$Mean.Utah.30,
                                    mid$Mean.Utah.31, mid$Mean.Utah.32,
                                    mid$Mean.Utah.33, mid$Mean.Utah.34,
                                    mid$Mean.Utah.35, mid$Mean.Utah.36,
                                    mid$Mean.Utah.37, mid$Mean.Utah.38),
                      
                      chillports = c(mid$Mean.Port.1, mid$Mean.Port.2,
                                     mid$Mean.Port.3, mid$Mean.Port.4,
                                     mid$Mean.Port.5, mid$Mean.Port.6,
                                     mid$Mean.Port.7, mid$Mean.Port.8,
                                     mid$Mean.Port.9, mid$Mean.Port.10,
                                     mid$Mean.Port.11, mid$Mean.Port.12,
                                     mid$Mean.Port.13, mid$Mean.Port.14,
                                     mid$Mean.Port.15, mid$Mean.Port.16,
                                     mid$Mean.Port.17, mid$Mean.Port.18,
                                     mid$Mean.Port.19, mid$Mean.Port.20,
                                     mid$Mean.Port.21, mid$Mean.Port.22,
                                     mid$Mean.Port.23, mid$Mean.Port.24,
                                     mid$Mean.Port.25, mid$Mean.Port.26,
                                     mid$Mean.Port.27, mid$Mean.Port.28,
                                     mid$Mean.Port.29, mid$Mean.Port.30,
                                     mid$Mean.Port.31, mid$Mean.Port.32,
                                     mid$Mean.Port.33, mid$Mean.Port.34,
                                     mid$Mean.Port.35, mid$Mean.Port.36,
                                     mid$Mean.Port.37, mid$Mean.Port.38),
                      
                      
                      gdd = c(mid$Mean.GDD.1, mid$Mean.GDD.2,
                              mid$Mean.GDD.3, mid$Mean.GDD.4,
                              mid$Mean.GDD.5, mid$Mean.GDD.6,
                              mid$Mean.GDD.7, mid$Mean.GDD.8,
                              mid$Mean.GDD.9, mid$Mean.GDD.10,
                              mid$Mean.GDD.11, mid$Mean.GDD.12,
                              mid$Mean.GDD.13, mid$Mean.GDD.14,
                              mid$Mean.GDD.15, mid$Mean.GDD.16,
                              mid$Mean.GDD.17, mid$Mean.GDD.18,
                              mid$Mean.GDD.19, mid$Mean.GDD.20,
                              mid$Mean.GDD.21, mid$Mean.GDD.22,
                              mid$Mean.GDD.23, mid$Mean.GDD.24,
                              mid$Mean.GDD.25, mid$Mean.GDD.26,
                              mid$Mean.GDD.27, mid$Mean.GDD.28,
                              mid$Mean.GDD.29, mid$Mean.GDD.30,
                              mid$Mean.GDD.31, mid$Mean.GDD.32,
                              mid$Mean.GDD.33, mid$Mean.GDD.34,
                              mid$Mean.GDD.35, mid$Mean.GDD.36,
                              mid$Mean.GDD.37, mid$Mean.GDD.38),
                      
                      mat.lo = c(mid$Spring.Temp.1, mid$Spring.Temp.2,
                                 mid$Spring.Temp.3, mid$Spring.Temp.4,
                                 mid$Spring.Temp.5, mid$Spring.Temp.6,
                                 mid$Spring.Temp.7, mid$Spring.Temp.8,
                                 mid$Spring.Temp.9, mid$Spring.Temp.10,
                                 mid$Spring.Temp.11, mid$Spring.Temp.12,
                                 mid$Spring.Temp.13, mid$Spring.Temp.14,
                                 mid$Spring.Temp.15, mid$Spring.Temp.16,
                                 mid$Spring.Temp.17, mid$Spring.Temp.18,
                                 mid$Spring.Temp.19, mid$Spring.Temp.20,
                                 mid$Spring.Temp.21, mid$Spring.Temp.22,
                                 mid$Spring.Temp.23, mid$Spring.Temp.24,
                                 mid$Spring.Temp.25, mid$Spring.Temp.26,
                                 mid$Spring.Temp.27, mid$Spring.Temp.28,
                                 mid$Spring.Temp.29, mid$Spring.Temp.30,
                                 mid$Spring.Temp.31, mid$Spring.Temp.32,
                                 mid$Spring.Temp.33, mid$Spring.Temp.34,
                                 mid$Spring.Temp.35, mid$Spring.Temp.36,
                                 mid$Spring.Temp.37, mid$Spring.Temp.38),
                      
                      siteslist = c(mid$Site.Num..1, mid$Site.Num..2,
                                    mid$Site.Num..3, mid$Site.Num..4,
                                    mid$Site.Num..5, mid$Site.Num..6,
                                    mid$Site.Num..7, mid$Site.Num..8,
                                    mid$Site.Num..9, mid$Site.Num..10,
                                    mid$Site.Num..11, mid$Site.Num..12,
                                    mid$Site.Num..13, mid$Site.Num..14,
                                    mid$Site.Num..15, mid$Site.Num..16,
                                    mid$Site.Num..17, mid$Site.Num..18,
                                    mid$Site.Num..19, mid$Site.Num..20,
                                    mid$Site.Num..21, mid$Site.Num..22,
                                    mid$Site.Num..23, mid$Site.Num..24,
                                    mid$Site.Num..25, mid$Site.Num..26,
                                    mid$Site.Num..27, mid$Site.Num..28,
                                    mid$Site.Num..29, mid$Site.Num..30,
                                    mid$Site.Num..31, mid$Site.Num..32,
                                    mid$Site.Num..33, mid$Site.Num..34,
                                    mid$Site.Num..35, mid$Site.Num..36,
                                    mid$Site.Num..37, mid$Site.Num..38),
                      year = (as.numeric(rownames(mid))+1970))
  
  site.mid<-full_join(middata, sites)
  site.mid$x<-NULL
  site.mid$y<-NULL
}

full.site<-full_join(site, site.post)
if(use.10yr==FALSE){full.site<-full_join(full.site, site.mid)}
full.site$year<-as.numeric(full.site$year)
full.site$cc <- NA

if(use.10yr==TRUE){
full.site$cc <- ifelse(full.site$year<=1960, "1950-1960", full.site$cc)
full.site$cc <- ifelse(full.site$year>2000 & full.site$year<=2010, "2000-2010", full.site$cc)
} else{
full.site$cc <- ifelse(full.site$year<=1970, "1950-1970", full.site$cc)
full.site$cc <- ifelse(full.site$year>1970 & full.site$year<=1990, "1970-1990", full.site$cc)
full.site$cc <- ifelse(full.site$year>1990 & full.site$year<=2010, "1990-2010", full.site$cc)
}
  
lodata <- subset(allpeps.subset, select=c("year", "lat", "long", "lo"))
full.site <- left_join(full.site, lodata)
full.site.nonas <- full.site[!is.na(full.site$lo),]


if(use.10yr==TRUE){
  write.csv(full.site.nonas, file="querob_allchillsandgdds_nomat_tenyr.csv", row.names = FALSE)
} else {
  write.csv(full.site.nonas, file="querob_allchillsandgdds_nomat_twentyyr.csv", row.names = FALSE)
}

##################################################################################################
############################### MEAN TEMP instead of GDD #########################################
##################################################################################################
#full.site.nonas <- read.csv("querob_allchillsandgdds_nomat_twentyyr.csv", header = TRUE)

if(use.10yr==TRUE){
period <- c(1951:1960, 2001:2010)
} else {
period <- c(1951:1970, 1971:1990, 1991:2010)
}

sites<-subset(full.site.nonas, select=c(lat, long, lat.long))
sites<-sites[!duplicated(sites$lat.long),]
sites$x<-sites$long
sites$y<-sites$lat
Coords<-subset(sites, select=c(x, y))
Coords <- na.omit(Coords)
nsites<-length(sites$lat.long)
tmin <- rn
tmax <- rx

points.min <- SpatialPoints(Coords, proj4string = rn@crs)
points.max <- SpatialPoints(Coords, proj4string = rx@crs)

yearsinclim<-as.numeric(format(as.Date(names(tmin),format="X%Y.%m.%d"),"%Y"))
yearsinperiod<-which(yearsinclim%in%period)
climsubmin<-subset(tmin,yearsinperiod)
climsubmax<-subset(tmax,yearsinperiod)

## subset climate days
monthsinclim<-as.numeric(format(as.Date(names(climsubmin),format="X%Y.%m.%d"),"%m"))
mstmonths<-c(3:4)
monthsinmst<-which(monthsinclim%in%mstmonths)
mstsubmin<-subset(climsubmin,monthsinmst)
mstsubmax<-subset(climsubmax,monthsinmst)

valuesmin <- raster::extract(mstsubmin,points.min)
valuesmax <- raster::extract(mstsubmax,points.max)

dclimmin <- cbind.data.frame(coordinates(points.min),valuesmin)
dclimmax <- cbind.data.frame(coordinates(points.max),valuesmax)

require(reshape2)
dxmin<-melt(dclimmin, id.vars=c("x","y"))
dxmax<-melt(dclimmax, id.vars=c("x","y"))

dxmin<-dxmin%>%
  rename(long=x)%>%
  rename(lat=y)%>%
  rename(date=variable)%>%
  rename(Tmin=value)

dxmax<-dxmax%>%
  rename(long=x)%>%
  rename(lat=y)%>%
  rename(date=variable)%>%
  rename(Tmax=value)

dx <- data.frame(lat=dxmin$lat, long=dxmin$long, date=dxmin$date, tmin=dxmin$Tmin, tmax=dxmax$Tmax)
dx$Tavg <- (dx$tmin+dx$tmax)/2

dx$date<-substr(dx$date, 2,11)
dx$Date<- gsub("[.]", "-", dx$date)

dx$tmin <- NULL
dx$tmax <- NULL
dx$date<-NULL

dx$year<-as.numeric(substr(dx$Date, 0, 4))
dx$lat.long<-paste(dx$lat, dx$long)
dx$month <- substr(dx$Date, 6,7)
dx$doy <- as.numeric(strftime(dx$Date, format = "%j"))

### Now, let's vary pre-season length. We'll add 30, 45 and 60 days
dx$mat60<-ave(dx$Tavg, dx$year, dx$lat.long)
dx$mat30 <- ifelse(dx$doy>=74 & dx$doy<=105, 
                   ave(dx$Tavg[dx$doy>=74 & dx$doy<=105], dx$year[dx$doy>=74 & dx$doy<=105], 
                       dx$lat.long[dx$doy>=74 & dx$doy<=105]), NA)

dx$mat45 <- ifelse(dx$doy>=60 & dx$doy<=105, 
                   ave(dx$Tavg[dx$doy>=60 & dx$doy<=105], dx$year[dx$doy>=60 & dx$doy<=105], 
                       dx$lat.long[dx$doy>=60 & dx$doy<=105]), NA)

dx <- na.omit(dx)

mst<-dx%>%dplyr::select(-Tavg, -Date, -doy, -month)
mst$id <- paste(mst$year, mst$lat.long)
mst<-mst[!duplicated(c(mst$id)),]
mst$id <- NULL

fullsites45 <- left_join(full.site, mst)
fullsites45 <- fullsites45[!is.na(fullsites45$lat.long),]

if(use.10yr==TRUE){
write.csv(fullsites45, file="~/Documents/git/decsens/analyses/pep_analyses/output/querob_decsens_1950_2000.csv", row.names = FALSE)
} else {
  write.csv(fullsites45, file="~/Documents/git/decsens/analyses/pep_analyses/output/querob_decsens_1950-2000.csv", row.names = FALSE)
  }

