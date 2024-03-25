install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel","ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
library("ggplot2")
theme_set(theme_bw())
library("sf")
world <- ne_countries(scale = "medium", returnclass = "sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("cowplot")
#Load in data
xJ<-read.csv("Japan.csv",header=TRUE)
xM<-read.csv("MAT7623.csv",header=TRUE)
xJ<-xJ[which(xJ$mag>=4.75),]
xM<-xM[which(xM$mag>=4),]
xM<-xM[1259:5393,]
xJ<-xJ[2600:4100,]

xEtas=c()
xJ$time<-gsub("[Z]","",xJ$time)
xJ$time<-gsub("[T]"," ",xJ$time)
xEtas$time<-as.numeric(difftime(xJ$time,xJ$time[1],units="days"))
xEtas$magnitude<-xJ$mag+0.1*runif(nrow(xJ))-0.05 ##Jitter magnitude for plotting
xEtas<-as.data.frame(xEtas)



xEtasM=c()
xM$time<-gsub("[Z]","",xM$time)
xM$time<-gsub("[T]"," ",xM$time)
xEtasM$time<-as.numeric(difftime(xM$time,xM$time[1],units="days"))
xEtasM$magnitude<-xM$mag+0.1*runif(nrow(xM))-0.05##Relative magbitude
xEtasM<-as.data.frame(xEtasM)


Mlong<-c(floor(min(xM$longitude)),ceiling(max(xM$longitude)))
Mlat<-c(floor(min(xM$latitude)),ceiling(max(xM$latitude)))

#Middle America Trench map
MAT<-ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = Mlong, ylim = Mlat, expand = FALSE)+
  ggtitle("B")+
  geom_point(data = xM, aes(x = xM$longitude, y = xM$latitude),col=1,cex=0.33)               

Jlong<-c(floor(min(xJ$longitude))-1,ceiling(max(xJ$longitude))+1)
Jlat<-c(floor(min(xJ$latitude))-1,ceiling(max(xJ$latitude))+1)

#Japan map
JAP<-ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = Jlong, ylim = Jlat, expand = FALSE)+
  ggtitle("A")+
  geom_point(data = xJ, aes(x = xJ$longitude, y = xJ$latitude),col=1,cex=0.33)               

#Japan magnitude time plot
JAPMT<-ggplot(xEtas,aes(x=time, y=magnitude))+geom_point(cex=0.3)+
  labs(x="Days since 07/12/93",y="Magnitude")+ggtitle("C")

#Middle America Trench magnitude time plot
MATmt<-ggplot(xEtasM,aes(x=time, y=magnitude))+geom_point(cex=0.3)+
  labs(x="Days since 01/13/98",y="Magnitude")+ggtitle("D")
#Into one plot
plot_grid(JAP, MAT,JAPMT,MATmt, labels = NULL)



###### Residual plot in main body of text
par(mfrow=c(3,2))
par(mar=c(4.1,4.1,1.6,2.1))
residSP1J=read.csv("Japan55SP1.csv",header=FALSE)
residSP1J<-residSP1J[,]
plot(1:length(residSP1J),residSP1J-(1:length(residSP1J)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-30,30))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP1J)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP1J)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP1J)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP1J)),col="gray",lty=3)
mtext("A",adj=0)


residSP1M<-read.table("MAT435SP1.csv",header=FALSE)
residSP1M<-residSP1M[,]
plot(1:length(residSP1M),residSP1M-(1:length(residSP1M)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-50,50))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP1M)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP1M)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP1M)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP1M)),col="gray",lty=3)
mtext("B",adj=0)



residSP2J=read.csv("Japan55SP2.csv",header=FALSE)
residSP2J<-residSP2J[,]
plot(1:length(residSP2J),residSP2J-(1:length(residSP2J)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-60,60))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP2J)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP2J)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP2J)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP2J)),col="gray",lty=3)
mtext("C",adj=0)

residSP2M<-read.table("MAT435SP2.csv",header=FALSE)
residSP2M<-residSP2M[,]
plot(1:length(residSP2M),residSP2M-(1:length(residSP2M)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-100,100))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP2M)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP2M)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP2M)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP2M)),col="gray",lty=3)
mtext("D",adj=0)

load("JapanETAStimefit.RData")
residJE<-residuals(x0)
plot(1:length(residJE),residJE-(1:length(residJE)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-120,65))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residJE)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residJE)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residJE)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residJE)),col="gray",lty=3)
mtext("E",adj=0)

load("MATETAStimefit.RData")
residJM<-residuals(x0)
plot(1:length(residJM),residJM-(1:length(residJM)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-105,105))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residJM)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residJM)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residJM)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residJM)),col="gray",lty=3)
mtext("F",adj=0)






######## Residual Plot in Appendix
par(mfrow=c(4,2))
par(mar=c(4.1,4.1,1.6,2.1))
residSP1J=read.csv("~/GitHub/CensoredLikelihood/MAT455SP1.csv",header=FALSE)
residSP1J<-residSP1J[,]
plot(1:length(residSP1J),residSP1J-(1:length(residSP1J)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-40,40))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP1J)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP1J)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP1J)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP1J)),col="gray",lty=3)
mtext("A",adj=0)

residSP1J=read.csv("Japan575SP1.csv",header=FALSE)
residSP1J<-residSP1J[,]
plot(1:length(residSP1J),residSP1J-(1:length(residSP1J)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-22,22))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP1J)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP1J)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP1J)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP1J)),col="gray",lty=3)
mtext("B",adj=0)

residSP1M<-read.table("MAT455SP2.csv",header=FALSE)
residSP1M<-residSP1M[,]
plot(1:length(residSP1M),residSP1M-(1:length(residSP1M)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-100,100))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP1M)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP1M)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP1M)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP1M)),col="gray",lty=3)
mtext("C",adj=0)

residSP2J=read.csv("Japan575SP2.csv",header=FALSE)
residSP2J<-residSP2J[,]
plot(1:length(residSP2J),residSP2J-(1:length(residSP2J)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-60,60))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP2J)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP2J)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP2J)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP2J)),col="gray",lty=3)
mtext("D",adj=0)

residSP2J=read.csv("MAT475SP1.csv",header=FALSE)
residSP2J<-residSP2J[,]
plot(1:length(residSP2J),residSP2J-(1:length(residSP2J)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-27,27))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP2J)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP2J)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP2J)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP2J)),col="gray",lty=3)
mtext("E",adj=0)

residSP2M<-read.table("Japan6SP1.csv",header=FALSE)
residSP2M<-residSP2M[,]
plot(1:length(residSP2M),residSP2M-(1:length(residSP2M)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-18,18))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP2M)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP2M)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP2M)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP2M)),col="gray",lty=3)
mtext("F",adj=0)

residSP2M<-read.table("MAT475SP2.csv",header=FALSE)
residSP2M<-residSP2M[,]
plot(1:length(residSP2M),residSP2M-(1:length(residSP2M)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-100,100))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP2M)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP2M)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP2M)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP2M)),col="gray",lty=3)
mtext("G",adj=0)

residSP2M<-read.table("Japan6SP2.csv",header=FALSE)
residSP2M<-residSP2M[,]
plot(1:length(residSP2M),residSP2M-(1:length(residSP2M)),
     xlab="CNE",ylab="MRTT",type="l",ylim=c(-60,60))
abline(h=0,col="gray")
abline(h=1.36*sqrt(length(residSP2M)),col="gray",lty=2)
abline(h=-1.36*sqrt(length(residSP2M)),col="gray",lty=2)
abline(h=1.63*sqrt(length(residSP2M)),col="gray",lty=3)
abline(h=-1.63*sqrt(length(residSP2M)),col="gray",lty=3)
mtext("H",adj=0)

## Compute ETAS intensity
#load appropriate output file e.g. "JapanETAStimefit.RData"
library("PtProcess")
ETASint<-etas_gif(xEtas,xEtas$time,x0$params)
mean(ETASint)
