#### Mu plot source code for ungulate example
spnum <- length(unique(distall$species))
quartz()
par(xpd=FALSE)
par(mar=c(5,10,3,10))
plot(x=NULL,y=NULL, xlim=c(-100,100), yaxt='n', ylim=c(0,6),
     xlab="Model estimate change in migration distance (km)", ylab="")
axis(2, at=1:6, labels=rev(c("Canadian", "Herbivore", "Canadian Herbivores",
                             "Sigma Canadian", "Sigma \nHerbivore", 
                             "Sigma Interaction")), las=1)
abline(v=0, lty=2, col="darkgrey")
rownameshere <- c("mu_b_canadian_sp", "mu_b_herb_sp", "mu_b_ch_sp", "sigma_b_canadian_sp",
                  "sigma_b_herb_sp", "sigma_b_ch_sp")
for(i in 1:6){
  pos.y<-(6:1)[i]
  pos.x<-noncps[rownameshere[i],"mean"]
  lines(noncps[rownameshere[i],c("25%","75%")],rep(pos.y,2),col="darkgrey")
  points(pos.x,pos.y,cex=1.5,pch=19,col="darkblue")
  for(spsi in 1:spnum){
    pos.sps.i<-which(grepl(paste0("[",spsi,"]"),rownames(noncps),fixed=TRUE))[3:5]
    jitt<-(spsi/40) + 0.08
    pos.y.sps.i<-pos.y-jitt
    pos.x.sps.i<-noncps[pos.sps.i[i],"mean"]
    lines(noncps[pos.sps.i[i],c("25%","75%")],rep(pos.y.sps.i,2),
          col=alpha(my.pal[spsi], alphahere))
    points(pos.x.sps.i,pos.y.sps.i,cex=0.8, pch=my.pch[spsi], col=alpha(my.pal[spsi], alphahere))
    
  }
}
par(xpd=TRUE) # so I can plot legend outside
legend(120, 6, sort(unique(gsub("_", " ", distall$species))), pch=my.pch[1:spnum],
       col=alpha(my.pal[1:spnum], alphahere),
       cex=1, bty="n", text.font=3)

