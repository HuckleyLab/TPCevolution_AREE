#AREE figures
library(ggplot2)

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

breadth= c(0.05,0.1,0.1) 
shift= c(-10,-10,-10)
arans=c(0,0,1)
scen=c("","spec-gen","thermodynamics")

vars= cbind(breadth, shift, arans)

for(k in 1:3){ 
  perf = TPC.beta(temps, shift=vars[k,2], breadth=vars[k,1], aran=vars[k,3])  
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
ggplot(pdat,aes(x=ts, y=ps, color=scen))+geom_line()



#--------
#B. Gilchrist 1995, envi var

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
p<-ggplot(tpcs)+aes(x=temperature, y = performance)+geom_line(aes(color=WG, lty=AG)) 
#p + geom_area(aes(fill=WG))
  
#--------
#C. Ashbury and Angilletta 2010, hotter is better

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

for(k in 1:9){ 
 perf = TPC.beta(temps, shift=vars[k,2], breadth=vars[k,1], aran=1, tolerance= 80, skew=0.5)  
  if(k==1) ps= perf
  if(k>1) ps= c(ps, perf)
 }

#other vars
ts= rep(temps,9)
wvs= rep(wv, each=60)
avs= rep(av, each=60)
ps[is.nan(ps)]=0

pdat=cbind(ts,ps, wvs, avs)
pdat= as.data.frame(pdat)
pdat$wvs= factor(pdat$wvs)
pdat$avs= factor(pdat$avs)
pdat$groups= paste(pdat$wvs,pdat$avs, sep="_")

#plot
ggplot(pdat,aes(x=ts, y=ps))+geom_line(aes(color=wvs, lty=avs))

#--------
#D. Williams et al. Carry over

#parameters from CarryoverICB
shifts= c(-8.416198, -9.609429, -7.507768, -12.56133, -0.8443629, -14.53758, -13.06035)
breadths= c(0.08008858, 0.06062592, 0.09135031, 0.1375593, 0.1351597, 0.15, 0.05)
  
#---
ltys= c("solid","dashed","dotted","solid","dashed","dotted", "solid")
cols= c("black","black","black","darkgrey","darkgrey","darkgrey", "grey")
lwds= c(3,3,3,3,3,3,5)

d <- density(clim.dat$V2)  
plot(d, xlim=range(5,45), main="", ylab="Performance", xlab="Temperature (?C)", col="gray", lty="solid",yaxt="n", cex.lab=1.2)
polygon(d, col="lightgray", border="lightgray")

for(z.count in 1:7){ 
  
  z1= TPC.beta(T, shift=shifts[z.count], breadth=breadths[z.count]) 
  par(new=TRUE)
  if(z.count==1){
    zscale= TPC.beta(temps, shift=shifts[z.count], breadth=min(breadths, na.rm=TRUE))
    plot(temps,z1, type="l", ylim=range(0,max(zscale)), xlim=range(5,45),  ylab="", xlab="", lty=ltys[z.count], lwd=lwds[z.count], col=cols[z.count], axes=FALSE)
  }
  if(z.count>1) points(temps,z1, type="l",lty=ltys[z.count], lwd=lwds[z.count], col=cols[z.count])
  
  #axis(side=4)
  
} #end loop plots

legend("topright", legend=c("noncumulative","hardening", "stress response"), lty=c("solid","dashed","dotted"), col=c("black","black","black") , bty="n")
#legend("bottomright", legend=c("injury","mortality"), lty=c("solid"), col=c("black","grey") , bty="n")

