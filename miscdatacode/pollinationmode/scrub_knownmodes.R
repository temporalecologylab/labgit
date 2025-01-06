### Started 19 October 2010 ###
### By Lizzie ###

## Updated 13 November 2010 to fix some errors associated with me mis-noting that species in the known modes is the full binomial ##
## Updated 11 April 2011 to update file names (to be safe we're re-writing) them

## Takes Winny's known modes file and runs it through phenology taxa datasets ##

## NOTES ##
## Be sure to update knownmodes.csv! ##

options(stringsAsFactors=FALSE)

library(Hmisc)

setwd("/Volumes/clelandlab/SharedData/PhenologyTraits/R")

# First, get a list of all species we have with family #
namesall <- unique(read.csv("names_all.csv", colClasses="character"))
pheno.species <- paste(capitalize(namesall$genus), tolower(namesall$species))
namesall$pheno.species <- paste(capitalize(namesall$genus),
   tolower(namesall$species))

ipni <- read.csv("TaxonScrubber/TS_AllTaxa.csv")
ipni.species <-paste(ipni$Genus, ipni$Species)
ipni$ipni.species <- paste(ipni$Genus, ipni$Species)

pnames.inIPNI <- namesall[pheno.species %in% ipni.species,]
ipnifam <- subset(ipni, select=c("Family", "ipni.species"))

masterlistsub <- ipnifam[ipni.species %in% pheno.species,]
masterlistsub$boo <- NA
masterlistsub1 <- aggregate(masterlistsub["boo"],
    masterlistsub[c("ipni.species","Family")], FUN=length)
masterlist <- merge(masterlistsub1, namesall, by.x="ipni.species",
    by.y="pheno.species", all.x=TRUE)

# Now get known modes file #
known <- read.csv("knownmodes2011.csv", header=TRUE)
names(known) # yuck

knownf <- subset(known, select=c("Family",
   "Resource", "Pollination.Mode" , "Notes"))
knowng <- subset(known, select=c("Genus",
   "Resource", "Pollination.Mode" , "Notes"))
knownsp <- subset(known, select=c(
   "Species..I.e..Genus...specific.epithet.",
   "Resource", "Pollination.Mode" , "Notes"))
knownsp$specieslower <- knownsp$Species..I.e..Genus...specific.epithet. # cleans up an old mistake

# Grab a site (you could do this for the whole masterlist but let's start here #
bertel <- subset(masterlist, site=="bertel"|site=="fitter")

# first, pull out families that are known #
bertel1 <- merge(bertel, knownf, by="Family", all.x=TRUE, suffixes=c("", "fam"))
bertel2 <- merge(bertel1, knowng, by.x="genus",
   by.y="Genus", all.x=TRUE, suffixes=c("", "gen"))
bertel3 <- merge(bertel2, knownsp, by.x="ipni.species",
   by.y="specieslower", all.x=TRUE, suffixes=c("", "sp"))

# Now, since this is totally hacked, clean up the mess #
bertel3$p.mode <- NA
bertel3$p.mode <- as.character(bertel3$p.mode)

for (fixer in c(1:nrow(bertel3))){
  if ((is.na(bertel3$Pollination.Mode[fixer]))==FALSE)
       bertel3$p.mode[fixer] <- bertel3$Pollination.Mode[fixer]
  else
    {
    if ((is.na(bertel3$Pollination.Modegen[fixer]))==FALSE)
    bertel3$p.mode[fixer] <- bertel3$Pollination.Modegen[fixer]
    else
      {
        if ((is.na(bertel3$Pollination.Modesp[fixer]))==FALSE)
        bertel3$p.mode[fixer] <- bertel3$Pollination.Modesp[fixer]
      }
  }

}

bertelout <- subset(bertel3, select=c("site", "Family", "genus", "species", "p.mode"))


## Read in work Winny has already done and merge ##
dater <- read.csv("/Volumes/clelandlab/SharedData/PhenologyTraits/R/Data_sweep2011.csv", header=TRUE)

dater1 <- subset(dater, pollination..wind.animal.self. !="")

donework <- merge(dater1, bertelout, by=c("site", "genus", "species"), all.y=TRUE)

for (j in c(1:nrow(donework))){
  if (is.na(donework$pollination..wind.animal.self.[j])==TRUE)
  donework$pollination..wind.animal.self.[j] <- donework$p.mode[j]
}

donework1 <- donework[,1:13]

write.csv(donework1, "/Volumes/clelandlab/SharedData/PhenologyTraits/R/outputApr112011.csv", row.names=FALSE)


## Now do it for the whole list ##
# first, pull out families that are known #
masterlist1 <- merge(masterlist, knownf, by="Family", all.x=TRUE, suffixes=c("", "fam"))
masterlist2 <- merge(masterlist1, knowng, by.x="genus",
   by.y="Genus", all.x=TRUE, suffixes=c("", "gen"))
masterlist3 <- merge(masterlist2, knownsp, by.x="ipni.species",
   by.y="specieslower", all.x=TRUE, suffixes=c("", "sp"))

# Now, since this is totally hacked, clean up the mess #
masterlist3$p.mode <- NA
masterlist3$p.mode <- as.character(masterlist3$p.mode)

for (fixer in c(1:nrow(masterlist3))){
  if ((is.na(masterlist3$Pollination.Mode[fixer]))==FALSE)
       masterlist3$p.mode[fixer] <- masterlist3$Pollination.Mode[fixer]
  else
    {
    if ((is.na(masterlist3$Pollination.Modegen[fixer]))==FALSE)
    masterlist3$p.mode[fixer] <- masterlist3$Pollination.Modegen[fixer]
    else
      {
        if ((is.na(masterlist3$Pollination.Modesp[fixer]))==FALSE)
        masterlist3$p.mode[fixer] <- masterlist3$Pollination.Modesp[fixer]
      }
  }

}

masterout <- subset(masterlist3,
    select=c("site", "Family", "genus", "species", "p.mode"))

write.csv(masterout, "/Volumes/clelandlab/SharedData/PhenologyTraits/R/masteroutputApr2011.csv", row.names=FALSE)

dim(subset(masterout, is.na(p.mode)))
notcape <- subset(masterout, site=="cape")
dim(subset(notcape, is.na(p.mode)))
dim(notcape)
