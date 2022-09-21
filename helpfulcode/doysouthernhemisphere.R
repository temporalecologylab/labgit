## Started on my brother's birthday, 2022 ##
## Code to translate doy in northern and southern hemispheres ##

# BTW, how get DOY for nothern hemisphere
# df$doy <- as.numeric(format(df$date , "%j"))

# Here we'll use DOY to get a sequence of dates
nhemi <- c(1:365)
nhemidate <- as.Date(paste(nhemi, 2002, sep="-"), format="%j-%Y") # I picked a random non-leap year

# Now, we start counting on 1 July for southern hemisphere
nhemidate[182]
shemi <- c(185:365, 1:184)

# Make a look up table
northsouth <- data.frame(caldate=nhemidate, northdoy=nhemi, southdoy=shemi)

# Write it out to your location
write.csv(northsouth, "/Users/lizzie/Documents/git/projects/vin/adelaideclimate/analyses/output/northsouthlookup.csv", row.names=FALSE)
