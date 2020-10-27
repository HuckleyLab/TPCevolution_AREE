#AREE figures
library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)

#FUNCTIONS
#Deutsch et al. TPC
#Performance Curve Function from Deutsch et al. 2008
tpc.plot= function(T,Topt,CTmin, CTmax){
  F=T
  F[]=NA
  sigma= (Topt-CTmin)/4
  F[T<=Topt & !is.na(T)]= exp(-((T[T<=Topt & !is.na(T)]-Topt)/(2*sigma))^2) 
  F[T>Topt & !is.na(T)]= 1- ((T[T>Topt & !is.na(T)]-Topt)/(Topt-CTmax))^2
  #set negetative to zero
  F[F<0]<-0
  
  return(F)
}

#Logan TPC
#Tb is body temperature
#Tbr is TPC breadth
#CTmax is critical thermal maxima
#p
ltpc=function(Tb,Tbr,CTmax,p) 4*(exp(p*(Tb-(CTmax-Tb)))-exp(p*Tbr-1,2*p*(CTmax-Tb)))

TPC.beta= function(T, shift=-1, breadth=0.1, aran=0, tolerance= 43, skew=0.7){ 
  T = T + 273.15 #Convert temperature in degrees Celsius to Kelvin
  shift= shift + 273.15 #Convert temperature in degrees Celsius to Kelvin         
  z=rep(0.01, length(T))
  z[which(is.na(T))]=NA  #account for NAs
  sel= which(T-shift>=0 & T-shift<=tolerance)
  z[sel]= ((T[sel]-shift)/tolerance)^(skew/breadth-1)*(1-(T[sel]-shift)/tolerance)^((1-skew)/breadth-1) / beta(skew/breadth,(1-skew)/breadth) 
  if(aran==1) z[sel]=z[sel]*exp(-0.6/(T[sel]*8.61734*10^(-5)))*10^10 #add scaling factor
  return(z)
}

#temp vector
temps=1:60

#--------
#Figure 1. 
#A. TPC, specialist generalist, hotter is better

breadth= c(0.1,0.2,0.2) 
shift= c(-10,-5,-7)
arans=c(0,0,1)
scen=c("cool","warm: specialist-generalist","warm: thermodynamics")

vars= cbind(breadth, shift, arans)

for(k in 1:3){ 
  perf = TPC.beta(temps, shift=vars[k,2], breadth=vars[k,1], aran=vars[k,3], tolerance=50)  
  if(k==1) ps= perf
  if(k>1) ps= c(ps, perf)
}

#other vars
ts= rep(temps,3)
ps[is.nan(ps)]=0

pdat=cbind(ts,ps, wvs, avs)
pdat= as.data.frame(pdat)
pdat$scen= rep(scen, each=60)

#plot
fig1a=ggplot(pdat,aes(x=ts, y=ps, color=scen))+geom_line(size=1.1)+
  theme_bw(base_size=14)+xlab("")+ylab("")+
  theme(legend.position = c(0.4, 0.9))+
  theme(legend.background = element_rect(fill=NA))+scale_color_viridis_d()+ 
  labs(color='') 

fig1a= fig1a +geom_text(data=pdat, x=5, y=0.15, label="CTmin", show.legend=FALSE, color="black",size=5)+
  geom_text(data=pdat, x=25, y=3, label="Topt", show.legend=FALSE, color="black",size=5)+
  geom_text(data=pdat, x=35, y=0.15, label="CTmax", show.legend=FALSE, color="black",size=5)

#--------
#Gilchrist 1995, envi var
#DON'T USE

#Parameters from Table 1
#AG and WG variation
# Tbr, p, CTmax, Topt
T_10_0.5= c(11.29, 0.26, 31.23, 27.72)
T_10_4.5= c(1.01, 5.80, 26.65, 26.49)
T_10_16.5= c(1.03, 5.66, 26.66, 26.49)

T_20_0.5= c(22.91, 0.09, 37.94, 28.04)
T_20_4.5= c(17.38, 0.14, 36.24, 29.70)
T_20_16.5= c(1.03, 5.66, 28.00, 27.84)

t1= cbind(10, 0.05, temps, tpc.plot(temps,T_10_0.5[4], T_10_0.5[3]-T_10_0.5[1], T_10_0.5[3]))
t2= cbind(10, 4.5, temps, tpc.plot(temps,T_10_4.5[4], T_10_4.5[3]-T_10_4.5[1], T_10_4.5[3]))
t3= cbind(10, 16.5, temps, tpc.plot(temps,T_10_16.5[4], T_10_16.5[3]-T_10_16.5[1], T_10_16.5[3]))
t4= cbind(20, 0.05, temps, tpc.plot(temps,T_20_0.5[4], T_20_0.5[3]-T_20_0.5[1], T_20_0.5[3]))
t5= cbind(20, 4.5, temps, tpc.plot(temps,T_20_4.5[4], T_20_4.5[3]-T_20_4.5[1], T_20_0.5[3]))
t6= cbind(20, 16.5, temps, tpc.plot(temps,T_20_16.5[4], T_20_16.5[3]-T_20_16.5[1], T_20_16.5[3]))

tpcs= rbind(t1,t2,t3,t4,t5,t6)
colnames(tpcs)=c("AG","WG","temperature","performance")
tpcs= as.data.frame(tpcs)
tpcs$AG= as.factor(tpcs$AG)
tpcs$WG= as.factor(tpcs$WG)

#plot
p<-ggplot(tpcs)+aes(x=temperature, y = performance)+geom_line(aes(color=WG, lty=AG)) +xlim(10,40)
#p + geom_area(aes(fill=WG))
  
#--------
#B. Ashbury and Angilletta 2010, hotter is better

#parameters
#just fecundity, fig 2
#breadth= c(0.002, 0.00035, 0.00033, 0.0116, 0.00014, 0.00013, 0.0217, 0.02, 0.0008) 
#shift= c(-20, -19, -16.5, -18.67, -18.0, -15.1, -17, -16.2, -14)

#fecund and survival, fig 3
breadth= c(0.005, 0.012, 0.03, 0.01, 0.02, 0.057, 0.016, 0.052, 0.127) 
shift= c(-20,-19,-20,-18.5,-20,-21,-17,-19,-20)

wv= c(0.5, 3.5, 6.5, 0.5, 3.5, 6.5, 0.5, 3.5, 6.5)
av= c(0,0,0,10,10,10,20,20,20)
vars= cbind(breadth, shift, wv, av)
#drop wv 0.5
inds=which(wv>0.5)
wv=wv[inds]
av=av[inds]
vars=vars[inds,]

for(k in 1:nrow(vars)){ 
 perf = TPC.beta(temps, shift=vars[k,2], breadth=vars[k,1], aran=1, tolerance= 80, skew=0.5)  
  if(k==1) ps= perf
  if(k>1) ps= c(ps, perf)
 }

#other vars
ts= rep(temps,nrow(vars))
wvs= rep(wv, each=60)
avs= rep(av, each=60)
ps[is.nan(ps)]=0

pdat=cbind(ts,ps, wvs, avs)
pdat= as.data.frame(pdat)
pdat$wvs= factor(pdat$wvs)
pdat$avs= factor(pdat$avs)
pdat$groups= paste(pdat$wvs,pdat$avs, sep="_")

#plot
fig1b= ggplot(pdat,aes(x=ts, y=ps))+geom_line(aes(color=wvs, lty=avs),size=1.1)+
theme_bw(base_size=14)+xlab("")+ylab("")+
  theme(legend.position = c(0.75, 0.75), legend.background = element_rect(fill="transparent"))+
  scale_color_viridis_d()+ 
  labs(color='within generation',lty='among generations')

#--------
#C. Williams et al. Carry over

#parameters from CarryoverICB
shifts= c(-8.416198, -9.609429, -7.507768, -12.56133, -0.8443629, -14.53758, -13.06035)
breadths= c(0.08008858, 0.06062592, 0.09135031, 0.1375593, 0.1351597, 0.15, 0.05)
  
#carryover= c("no carryover","beneficial acclimation", "cumulative damage","no carryover","beneficial acclimation", "cumulative damage","no carryover")
carryover= c("no carryover","acclimation", "damage","no carryover","acclimation", "damage","no carryover")
component=c("injury","injury","injury","mortality","mortality","mortality","performance")
dat= data.frame(shifts, breadths, carryover, component)

for(i in 1:nrow(vars) ){ 
  
  p= TPC.beta(temps, shift=dat[i,"shifts"], breadth=dat[i,"breadths"])
  p1=as.data.frame(cbind(temps,p))
  p1$carryover=dat$carryover[i]
  p1$component=dat$component[i]
  
  if(i==1) p1.all= p1
  if(i>1) p1.all=rbind(p1, p1.all)
}
  
#plot
p1.all$carryover=factor(p1.all$carryover, levels=c("no carryover","acclimation","damage"))
p1.all$component=factor(p1.all$component, levels=c("performance","injury","mortality"))

fig1c= ggplot()+
  geom_density(data=clim.dat, aes(x=V2,y=..scaled..),color="grey", fill="grey")+
  geom_line(data=p1.all, aes(x=temps, y=p, color=component, lty=carryover),size=1.1)+
  theme_bw(base_size=14)+xlab("")+ylab("")+
  theme(legend.position = c(0.8, 0.75),legend.background = element_rect(fill="transparent"))+
  theme(legend.background = element_rect(fill=NA))+scale_color_viridis_d()

#========================
#combine
setwd("/Volumes/GoogleDrive/Shared drives/TrEnCh/Projects/AREE_TPCevolution/figures/")

pdf("Fig1.pdf", height = 6, width = 12)

plot<-plot_grid(fig1a, fig1b, fig1c, align='vh', ncol=3, vjust=1, scale = 1,labels=c("a","b","c"))

#create common x and y labels
y.grob <- textGrob("performance", 
                   gp=gpar(fontsize=20), rot=90)
x.grob <- textGrob("temperature (°C)", 
                   gp=gpar(fontsize=20))

#add to plot
grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))

dev.off()

