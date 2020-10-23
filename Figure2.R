#AREE figures
library(ggplot2)
library(viridis)

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

#=============================
#Fig 2.
temps=0:50

#Colias butterfly larval feeding
#Higgins
tpc= function(T, Fmax,To, row, sigma) Fmax*exp(-exp(row*(T-To)-6)-sigma*(T-To)^2)

setwd("/Volumes/GoogleDrive/Shared drives/TrEnCh/Projects/AREE_TPCevolution/data/")
dat= read.csv("HigginsTPC.csv")
#keep comparison
dat= dat[c(1,3,5,6),]

for(i in 1:nrow(dat)){
p= tpc(T=temps, Fmax=dat$Fmax[i],To=dat$Topt[i], row=dat$row[i], sigma=dat$sigma[i])
p1=as.data.frame(cbind(temps,p))
p1$taxa="butterfly"
p1$species=dat$species[i]
p1$population=dat$population[i]
p1$year= dat$year[i]

if(i==1) p1.all= p1
if(i>1) p1.all=rbind(p1, p1.all)
}

#plot
p1.all$year= as.factor(p1.all$year)
fig2.butterfly= ggplot(p1.all,aes(x=temps, y=p))+geom_line(aes(color=population, lty=year),size=2)+
  theme_bw(base_size=14)+xlab("temperature (°C)")+ylab("performance")

#------
# Ant urban adaptation https://academic.oup.com/biolinnean/article/121/2/248/3038290
#Evol App https://datadryad.org/stash/dataset/doi:10.5061/dryad.0r75421
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.0r75421
# data in dryad

setwd("/Volumes/GoogleDrive/Shared drives/TrEnCh/Projects/AREE_TPCevolution/data/ants/")

ctmin= read.csv("Field_versus_F2_CTmin.csv")
ctmax= read.csv("Field_versus_F2_CTmax.csv")

ctmin= aggregate(ctmin$ctmin, by=list(colonyID=ctmin$colonyID,generation=ctmin$generation,population=ctmin$population), FUN="mean")
ctmax= aggregate(ctmax$ctmax, by=list(colonyID=ctmax$colonyID,generation=ctmax$generation,population=ctmax$population), FUN="mean")
names(ctmin)[4]<-"ctmin"
names(ctmax)[4]<-"ctmax"
dat=cbind(ctmin,ctmax$ctmax)
#match(ctmin$colonyID, ctmax$colonyID)

for(i in 1:nrow(dat)){
  topt= dat$ctmin+(dat$ctmax-dat$ctmin)*0.7
  p= tpc.plot(T=temps, Topt=topt[i], CTmin=dat$ctmin, CTmax=dat$ctmax)
  p1=as.data.frame(cbind(temps,p))
  p1$taxa="butterfly"
  p1$colonyID=dat$colonyID[i]
  p1$generation=dat$generation[i]
  p1$population= dat$population[i]
  
  if(i==1) p1.all= p1
  if(i>1) p1.all=rbind(p1, p1.all)
}

#plot
ggplot(p1.all,aes(x=temps, y=p))+geom_line(aes(color=colonyID, lty=population))

#------
# Lizard TPCs https://onlinelibrary.wiley.com/doi/full/10.1111/1749-4877.12309c
# data in supplement

tpc= function(T,Topt,CTmin, CTmax, Pmax){
  F=T
  F[]=NA
  sigma= (Topt-CTmin)/4
  F[T<=Topt & !is.na(T)]= exp(-((T[T<=Topt & !is.na(T)]-Topt)/(2*sigma))^2)*Pmax 
  F[T>Topt & !is.na(T)]= Pmax- ((T[T>Topt & !is.na(T)]-Topt)/(Topt-CTmax))^2
  #set negetative to zero
  F[F<0]<-0
  
  return(F)
}

setwd("/Volumes/GoogleDrive/Shared drives/TrEnCh/Projects/AREE_TPCevolution/data/")
dat= read.csv("lizardTPC.csv")

for(i in 1:nrow(dat)){
  p= tpc(T=temps, Topt=dat$Topt[i],CTmin=dat$CTmin[i], CTmax=dat$CTmax[i], Pmax=dat$Pmax[i] )
  p1=as.data.frame(cbind(temps,p))
  p1$taxa="lizard"
  p1$pop=dat$pop[i]
  p1$gen= dat$gen[i]
  
  if(i==1) p1.all= p1
  if(i>1) p1.all=rbind(p1, p1.all)
}

#plot
ggplot(p1.all,aes(x=temps, y=p))+geom_line(aes(color=gen, lty=pop))+
  theme_bw(base_size=14)+xlab("temperature (°C)")+ylab("performance")

#------
#Plants
#Native invasive mimulus: https://www.biorxiv.org/content/10.1101/2020.09.10.291252v1.full [FIG?]

#Mimulus temp resurrection: https://onlinelibrary.wiley.com/doi/10.1111/evo.14041
#https://github.com/rwoolive/Cardinalis_TPC_evolution

dat= read.csv("MimulusTPc.csv")
#pick one population from each
dat= dat[dat$Pop %in% c("C1","N1","S1"),]
dat$Pop[dat$Pop=="C1"]<-"central"
dat$Pop[dat$Pop=="N1"]<-"north"
dat$Pop[dat$Pop=="S1"]<-"south"

for(i in 1:nrow(dat)){
  p= tpc(T=temps, Topt=dat$Topt[i], CTmin=dat$LTL[i], CTmax=dat$UTL[i],Pmax=dat$PerformanceMax[i])
  p1=as.data.frame(cbind(temps,p))
  p1$taxa="mimulus"
  p1$population=dat$Pop[i]
  p1$year= dat$year[i]
  
  if(i==1) p1.all= p1
  if(i>1) p1.all=rbind(p1, p1.all)
}

#plot
p1.all$year= as.factor(p1.all$year)
fig2.mimulus=ggplot(p1.all,aes(x=temps, y=p))+geom_line(aes(color=population, lty=year))+
  theme_bw(base_size=14)+xlab("temperature (°C)")+ylab("performance")

#------
#Microbe
#Phytoplankton fig https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12545 Fig 3

#Diatom fig https://www.nature.com/articles/s41467-018-03906-5 Fig 1

tpc= function(T, Ea, k=8.62*10^{-5}, mu, Tc=291.15, Eh, Th) exp(Ea*(1/(k*Tc)-1/(k*T))+mu-log(1+exp(Eh*(1/(k*Th)-1/(k*T)))) )

dat= read.csv("DiatomTPC.csv")
#change ancestor label
dat$environment[dat$time=="ancestor"]<-"ancestor"

temps=273:(273+50)

for(i in 1:nrow(dat)){
  p= tpc(T=temps, Ea=dat$Ea[i], mu=dat$mu[i], Eh=dat$Eh[i], Th=dat$Th[i])
  p1=as.data.frame(cbind(temps,p))
  p1$taxa="diatom"
  p1$environment=dat$environment[i]
  p1$time=dat$time[i]
  
  if(i==1) p1.all= p1
  if(i>1) p1.all=rbind(p1, p1.all)
}

#plot
p1.all$environment=factor(p1.all$environment, levels=c("ancestor","22","26","32","FS"))
cols= viridis(5)
cols[1]<-"darkgray"

fig2.diatom=ggplot(p1.all,aes(x=temps-273, y=p))+geom_line(aes(color=environment),size=2)+
  theme_bw(base_size=14)+xlab("temperature (°C)")+ylab("performance")+
  scale_colour_manual(values = cols)


