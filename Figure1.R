#AREE figures

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

#--------
#Figure 1. 
#A. TPC, specialist generalist, hotter is better

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

#--------
#D. Williams et al. Carry over

#--------
#E. Kingsolver and woods

