#AREE figures
library(ggplot2)
library(cowplot)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(gridExtra)
library(grid)
library(patchwork)

#source("analysis\\TempcyclesAnalysis.R")

#ROBOMUSSEL ANALYSIS

#SITES
# WaOr Tatoosh Island, WA 1660.1 48.39 124.74
# WaOr Boiler Bay, OR 1260.7 44.83 124.05
# WaOr Strawberry Hill, OR 1196 44.25 124.12
# CenCal Hopkins, CA 327.1 36.62 121.90
# CenCal Piedras Blancas, CA 208.11 35.66 121.28
# CenCal Cambria, CA 185.66 35.54 121.10
# SoCal Lompoc, CA 84.175 34.72 120.61
# SoCal Jalama, CA 57.722 34.50 120.50
# SoCal Alegria, CA 37.284 34.47 120.28
# SoCal Coal Oil Point (COP), CA 0 34.41 119.88

#-----------------
#Site data
#setwd( "C:\\Users\\Buckley\\Documents\\ClimateBiology\\")
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/ICBClimateBiology/data/")

site.dat= read.csv("musselREADME.csv")

#Load robomussel data
te.max <- readRDS("tedat.rds")
te.max= ungroup(te.max)
te.max= as.data.frame(te.max)

#Fix duplicate CP in WA
te.max$lat= as.character(te.max$lat)
te.max$site= as.character(te.max$site)
te.max[which(te.max$lat==48.45135),"site"]<-"CPWA"
te.max$site= as.factor(te.max$site)

#drop 1 of two close sites
te.max<- te.max[-which(te.max$site=="LB"),]

#----------------------
#PLOTS

#Time series
#clim2 = clim2 %>% group_by(Year,Site) %>% summarise(Min= mean(Min, na.rm=TRUE),Max= mean(Max, na.rm=TRUE),Mean= mean(Mean, na.rm=TRUE) )

#count by site, subsite, year
te.count = te.max %>% group_by(year,site) %>% summarise( count=length(MaxTemp_C)  )
te.count= as.data.frame(te.count)

#subset sites
te.max2= subset(te.max, te.max$site %in% c("BB","LL") ) # "HS","SD","AG"
#te.max2= subset(te.max, te.max$lat %in% c(48.39137,44.83064,35.66582,34.46717) )

# USWACC	48.5494	-123.0059667	Colins Cove
# USWACP	48.45135	-122.9617833	Cattle Point
#* USWASD	48.39136667	-124.7383667	Strawberry Point
#* USORBB	44.83064	-124.06005	Boiler Bay
#* USCAPD	35.66581667	-121.2867167	Piedras
# USCAAG	34.46716667	-120.2770333	Alegria

#time series
te.max1= te.max2
#te.max1= subset(te.max2, te.max2$year==2002)
#2002

#restrict to summer
#May 1 through September: 121:273 
te.max1= subset(te.max1, te.max1$doy>120 & te.max1$doy<274)

#ggplot(data=te.max1, aes(x=doy, y = MaxTemp_C, color=subsite ))+geom_line() +theme_bw()+facet_wrap(~site)
#by tidal height
#ggplot(data=te.max1, aes(x=doy, y = MaxTemp_C, color=height ))+geom_line() +theme_bw()+facet_wrap(~site)

#update labels
te.max1$labs= as.character(te.max1$site)
te.max1$labs[which(te.max1$labs=="BB")]<- "Boiler Bay, OR, 44.8°N"
te.max1$labs[which(te.max1$labs=="LL")]<- "Lompoc, CA, 34.7°N"

#restrict zones
te.max1$zone2= factor(te.max1$zone, levels=c("Low","Mid","High") )
#order sites
te.max1$labs= factor(te.max1$labs, levels=c("Boiler Bay, OR, 44.8°N", "Lompoc, CA, 34.7°N", ordered="TRUE")  )
  
te.max1= subset(te.max1, !is.na(te.max1$zone2==2002))

te.max.yr= subset(te.max1, te.max1$year==2008)

#FIG 1A
#by lat
sens.h <- data.frame(doy = 230,MaxTemp_C = 41,subsite=10, zone2="Mid",labs = "Boiler Bay, OR, 44.8°N")
sens.l <- data.frame(doy = 238,MaxTemp_C = 35,subsite=10, zone2="Mid",labs = "Boiler Bay, OR, 44.8°N")

fig.t<- ggplot(data=te.max.yr, aes(x=doy, y = MaxTemp_C, color=zone2, group=subsite))+geom_line(alpha=0.8) +
  theme_classic() +
  labs(x = "Day of year",y="Maximum daily temperature (°C)")+
  scale_color_manual(values=c("purple","darkgreen","darkorange"), name="intidal height")+
  geom_abline(aes(intercept=36.2, slope=0), col="purple")+ #low, fast
  geom_abline(aes(intercept=39.4, slope=0), col="darkorange")+ #high, fast
  geom_abline(aes(intercept=36.6, slope=0), col="darkorange", lty="dotted")+ #high slow
  geom_abline(aes(intercept=39.4-1.07, slope=0), col="darkorange", lty="dashed")+ #high, fast, acclimated
  theme(legend.position="bottom")+
  ggtitle('b')+
  geom_text(data = sens.h,label = "sensitivity: Tcrit high intertidal", color="darkorange")+
  geom_text(data = sens.l,label = "Tcrit low intertidal", color="purple")+
  facet_wrap(~labs, ncol=1)
#solid: fast heating, dashed: fast heating with acclimation, dotted: slow heating

#Just sensitivity for lompoc
#geom_abline(data = subset(te.max.yr, labs == "34.7° Lompoc, CA"), aes(intercept=39.4-1.07, slope=0), col="darkorange", lty="dashed")+ #high, fast, acclimated
  
#===================================================
#Quilt plot

#round lat
te.max$lat= as.numeric(as.character(te.max$lat))
te.max$lat.lab= round(te.max$lat,2)

#mean daily maximum by month
te.month = te.max %>% group_by(lat, month, lat.lab) %>% summarise( max=max(MaxTemp_C), mean.max=mean(MaxTemp_C), q75= quantile(MaxTemp_C, 0.75), q95= quantile(MaxTemp_C, 0.95) ) 

fig.quilt<- ggplot(te.month) + 
  aes(x = month, y = as.factor(lat.lab) ) + 
  geom_tile(aes(fill = mean.max)) + 
  coord_equal()+
  scale_fill_gradientn(colours = viridis(10), name="temperature (°C)" )+
  #scale_fill_distiller(palette="Spectral", na.value="white", name="max temperature (°C)") + 
  theme_classic()+xlab("month")+ylab("latitude (°)")+ theme(legend.position="bottom")+ #+ coord_fixed(ratio = 4)
  geom_hline(yintercept = 7.5, color="white", lwd=2) +
  scale_x_continuous(breaks=seq(1,12,2))+
  ggtitle('c')

#==================================================
#annual means
te.ann = te.max %>% group_by(year, lat, lat.lab) %>% summarise( max=max(MaxTemp_C), mean.max=mean(MaxTemp_C), q75= quantile(MaxTemp_C, 0.75), q95= quantile(MaxTemp_C, 0.95) ) 

#line
ggplot(te.ann[te.ann$lat.lab %in% c("44.83","34.72"),]) + 
  aes(x = year, y = mean.max, color=lat.lab, group=lat.lab) +geom_line()
  #just boiler bay and lompoc?

#=========================
#sensitivity

#Moyen: https://jeb.biologists.org/content/222/17/jeb203166
#(1) the critical temperature (Tcrit) at which heart rate (HR) precipitously declines
#zone
#https://jeb.biologists.org/content/223/13/jeb222893.full

#As heating rate significantly impacted high- but not low-zone mussels' cardiac thermal tolerance
#(1) the critical temperature (Tcrit) at which heart rate (HR) precipitously declines

#heating rates vary wave exposed / protected, height on shore
#hypothesized that, because high-zone mussels have a higher cardiac thermal tolerance and experience faster mean daily heating rates than low-zone mussels

#geo var: https://www.int-res.com/articles/meps2012/450/m450p093.pdf

#---------------
#combine
setwd("/Volumes/GoogleDrive/Shared drives/TrEnCh/Projects/AREE_TPCevolution/figures/")

#image
library(ggpubr)
library(magick)
ir <- image_read('VisIR_SandsIntertidal.JPG')
ir <- ggplot() +
  background_image(ir) + #coord_fixed()+
  ggtitle('a')


pdf("Fig0.pdf", height = 10, width = 8)
 (ir | fig.quilt)/fig.t 
dev.off()
