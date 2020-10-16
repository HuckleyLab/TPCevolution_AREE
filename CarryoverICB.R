 #rm(list=ls(all=TRUE))
library(msm) #for rtnorm
library(EnvStats) #for truncated lognormal

#source functions
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Extremes/Analysis/EvolModel/")
source("EvolModel_functions_10Aug2015.R")

#read MEBOURNE DATA
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Extremes/ExtremesSynched2/ExtremesEvol/analysis/")
clim.dat= read.table("acorn.sat.maxT.086071.MELBOURNE.daily.txt", header=FALSE, na.strings="99999.9", skip=1)
clim.dat= na.omit(clim.dat)

library(ks)
mel.d <- kde(clim.dat$V2)  
plot(mel.d, xlim=range(0,60), main="", ylab="Performance", xlab="Temperature (?C)", col="gray", lty="solid",yaxt="n", cex.lab=1.2)
#temp.dist=rkde(fhat=d, n=1000)

#-------------------------------------------------
#SURVIVAL FUNCTIONS
#SURVIVAL COST ABOVE CTmax
#Survival
surv.mat<- function(t.mat1, td=8.68){
  CTmin=t.mat1[1]; CTmax=t.mat1[2]; T=t.mat1[3:length(t.mat1)]   
  
  s1= ifelse(T<CTmax, s<-1, s<- exp(-(T-CTmax)/td) ) ## exponential
  #s1= ifelse(T<CTmax, s<-1, s<- 1+1/(CTmax-60)*(T-CTmax) ) #linear
  
  s1[which(s1<=0)]=0.01 #correct negative survivals
  s1[which(s1>1)]=1
  
  return( prod(s1) )
}

#t.mat1= t.mat[,1]
surv.mat.stress<- function(t.mat1, td=8.68){
  CTmin=t.mat1[1]; CTmax=t.mat1[2]; T=t.mat1[3:length(t.mat1)]  
  
  ext= ifelse(T<CTmax, s<-0, s<- 1 )
  exts= cumsum(ext)*ext
  exts.fact= 1*ext+ exts*(-0.01) #percent less survival each extreme
  exts.fact[which(exts.fact<0)]=0
  
  s1= ifelse(T<CTmax, s<-1, s<- exp(-(T-CTmax)/td) ) ## exponential
  #s1= ifelse(T<CTmax, s<-1, s<- 1+1/(CTmax-60)*(T-CTmax) ) #linear
  
  ext.ind= which(s1<1)
  s1[ext.ind]= s1[ext.ind]*exts.fact[ext.ind] #multiply survival of extremes by factor
  s1[which(s1<=0)]=0.01 #correct negative survivals
  s1[which(s1>1)]=1
  
  return( prod(s1) )
}

surv.mat.hard<- function(t.mat1, td=8.68){
  CTmin=t.mat1[1]; CTmax=t.mat1[2]; T=t.mat1[3:length(t.mat1)]  
  
  ext= ifelse(T<CTmax, s<-0, s<- 1 )
  exts= cumsum(ext)*ext
  exts.fact= 1*ext+ exts*(0.01) #percent more survival each extreme
  exts.fact[which(exts.fact<0)]=0
  
  s1= ifelse(T<CTmax, s<-1, s<- exp(-(T-CTmax)/td) ) ## exponential
  #s1= ifelse(T<CTmax, s<-1, s<- 1+1/(CTmax-60)*(T-CTmax) ) #linear
  
  ext.ind= which(s1<1)
  s1[ext.ind]= s1[ext.ind]*exts.fact[ext.ind] #multiply survival of extremes by factor
  s1[which(s1<=0)]=0.01 #correct negative survivals
  s1[which(s1>1)]=1
  
  return( prod(s1) )
}

#----------------------------------
## DAMAGE
#t.matz1= t.matz[,1]

damage.mat<- function(t.matz1, td=8.68, pdet ){
  CTmin=t.matz1[1]; CTmax=t.matz1[2]; T=t.matz1[3:302] ; z=t.matz1[303:length(t.matz1)]  
  
  ext= ifelse(T<CTmax, s<-0, s<- 1 )
  exts= cumsum(ext)
  #LOSS EACH EXTREME
  ploss.fact= 1-pdet*exts #performance cost of extreme
  ploss.fact[which(ploss.fact<0)]=0
  
  return( sum(z*ploss.fact) )
}

damage.mat.hard<- function(t.matz1, td=8.68, pdet ){
  CTmin=t.matz1[1]; CTmax=t.matz1[2]; T=t.matz1[3:302] ; z=t.matz1[303:length(t.matz1)]  
  
  ext= ifelse(T<CTmax, s<-0, s<- 1 )
  exts= cumsum(ext)
  #LOSS EACH EXTREME
  ploss.fact= 1-(1-0.02*exts)*pdet*exts #performance cost of extreme
  ploss.fact[which(ploss.fact>1)]=1
  ploss.fact[which(ploss.fact<0)]=0
  
  return( sum(z*ploss.fact) )
}

damage.mat.stress<- function(t.matz1, td=8.68, pdet ){
  CTmin=t.matz1[1]; CTmax=t.matz1[2]; T=t.matz1[3:302] ; z=t.matz1[303:length(t.matz1)]  
  
  ext= ifelse(T<CTmax, s<-0, s<- 1 )
  exts= cumsum(ext)
  #LOSS EACH EXTREME
  ploss.fact= 1-(1+0.02*exts)*pdet*exts #performance cost of extreme
  
  ploss.fact[which(ploss.fact>1)]=1
  ploss.fact[which(ploss.fact<0)]=0
  
  return( sum(z*ploss.fact) )
}


#==========================================================
### ANALYSIS

#SET up DATA STORAGE
shift.means= matrix(NA, ncol=7, nrow=1)
breadth.means= matrix(NA, ncol=7, nrow=1)

#set up possible performance curves
#define combinations
shift=seq(-15,4,1) #change from negative 10
breadth=seq(0.05,0.15,0.02)
comb.mb= expand.grid(shift, breadth)
comb.mb= as.matrix( cbind(as.numeric(comb.mb[,1]), as.numeric(comb.mb[,2])) )

#----------------------------------------------
#MODEL EVOLUTION

#Specify population size
Nind=200

#SHIFT BREADTH G MATRIX
#g= matrix(data = c(0.3, 0.1, 0.1, 0.3), nrow = 2) #Roughly based on Kingsolver et. al 2004
g= matrix(data = c(0.7, -0.1, -0.1, 0.7), nrow = 2) 

#Specify initial mean and variance shift and breadth, popn sample size
shift.sd=1
breadth.sd=0.02

#MODEL 4 SCENARIOS
# WITH / WITHOUT EXTREMES
#TROPICAL TEMPERATE

for(k in 1:1){  #LOOP SCENARIOS
print(k)

mean.temp= 20
sd.temp= 5

n=1000 #Number times steps for environmental data
sd.therm.hetero= 2

#loop through seven survival scenarios
zs=1:7

for(z.count in zs){

for(gen in 1:200){ 

#simulate climate
#various distributions perviously tried 
#temp.dist=rlnorm(n, meanlog = 3, sdlog = 0.4)
#temp.dist= rnorm(n, mean = mean.temp, sd = sd.temp)
#temp.dist=rlnormTrunc(n, meanlog = 3, sdlog = 0.3, min = 0, max = 60)
#temp.dist=rcauchy(n, location = mean.temp, scale = 0.5) 
#temp.dist= 15+ rexp(n, rate = 1/7) #var 1/(rate^2)

#SAMPLE FROM LEMBOUNE TEMP DISTRIBUTION
set.seed(1)
temp.dist=rkde(fhat=mel.d, n=n)
#d <- density(temp.dist)  
#plot(d, xlim=range(0,60), main="", ylab="Performance", xlab="Temperature (?C)", col="gray", lty="solid",yaxt="n", cex.lab=1.2)

#itexp <- function(u, m, t) { -log(1-u*(1-exp(-t*m)))/m }
#rtexp <- function(n, m, t) { itexp(runif(n), m, t) }

#set initial parameters 
if(gen==1 & k %in% c(1,4)){
  optim= comb.mb[which.max(apply(comb.mb, MARGIN=1, FUN=TPC.beta.sum, temps=temp.dist)),]
  shift.mean.init=optim[1]
  breadth.mean.init=optim[2]
  }

if(gen==1){shift.mean=shift.mean.init; breadth.mean=breadth.mean.init}

#### SIMULATE INDIVIDUALS
shift.sample <- rtnorm(Nind, mean = shift.mean, sd = shift.sd, lower=-10, upper=4)
breadth.sample <- rtnorm(Nind, mean = breadth.mean, sd = breadth.sd, lower=0.05, upper=0.15)
shift.breadth = as.matrix(cbind(shift.sample, breadth.sample))

#calculate Topt
ts= seq(0,60, 0.2)
topts= apply(shift.breadth, FUN= Topt.beta.mat,T=ts, MARGIN=1)
t.opt= matrix( rep(topts, length(temp.dist)), nrow= length(temp.dist), ncol=Nind, byrow = TRUE) 
tcts= apply(shift.breadth, FUN= Tcts.beta.mat,T=ts, MARGIN=1)

#****************************************
## CHANGE TO MATRIX CALCULATION

t.het= matrix( rep(temp.dist, Nind), nrow= length(temp.dist), ncol=Nind) 

#add random thermal heterogeneity for each time step
therm.var= rnorm( length(temp.dist)*Nind, mean=0, sd= sd.therm.hetero)
therm.var= matrix(therm.var,  nrow= length(temp.dist), ncol=Nind, byrow = FALSE)

#bound potential environments
t.low= t.het-abs(therm.var)
t.high= t.het+abs(therm.var)

#find optimal temp with behavior
t.sel<- t.het
t.sel[]=t.opt
inds= which(t.high<t.opt, arr.ind=T)
t.sel[inds]=t.high[inds]
inds= which(t.low>t.opt, arr.ind=T)
t.sel[inds]=t.low[inds] 

######NO BEHAVIOR
#t.sel<-t.het

#calculate performance
tsb= rbind(shift.sample, breadth.sample, t.sel)
z= apply(tsb, FUN=TPC.beta.mat, MARGIN=2)

#sum or products
z.sum= colSums(z)
z.prod= colSums( log(z))
#-----------------

### VARY COST OF EXCEEDING CTmax
t.mat= rbind(tcts, t.sel)
t.matz= rbind(tcts, t.sel,z)

#No cost
#if(z.count==1) perf.mat= z.sum

#Performance cost with hardening /stress
if(z.count==1){  perf.mat= apply(t.matz, FUN=damage.mat, MARGIN=2, pdet=0.02)}
 #   ploss= apply(t.mat, FUN=perfloss.mat, MARGIN=2, pdet=max(z))
 #   perf.mat= z.sum+ploss
#    perf.mat[which(perf.mat<0)]=0}
 

if(z.count==2){ perf.mat= apply(t.matz, FUN=damage.mat.hard, MARGIN=2, pdet=0.02)}
 # ploss= apply(t.mat, FUN=perfloss.mat.hard, MARGIN=2, pdet=max(z))
#  perf.mat= z.sum+ploss
#  perf.mat[which(perf.mat<0)]=0}
  
if(z.count==3){ perf.mat= apply(t.matz, FUN=damage.mat.stress, MARGIN=2, pdet=0.02)}
 # ploss= apply(t.mat, FUN=perfloss.mat.stress, MARGIN=2, pdet=max(z))
#  perf.mat= z.sum+ploss
#  perf.mat[which(perf.mat<0)]=0}

#Survival cost
if(z.count==4){
  surv= apply(t.mat, FUN=surv.mat, MARGIN=2)
  surv[which(surv<=0)]=0.001
  perf.mat= z.sum*surv}

if(z.count==5){
  surv= apply(t.mat, FUN=surv.mat.hard, MARGIN=2)
  surv[which(surv<=0)]=0.001
  perf.mat= z.sum*surv
}

if(z.count==6){
  surv= apply(t.mat, FUN=surv.mat.stress, MARGIN=2)
  surv[which(surv<=0)]=0.001
  perf.mat= z.sum*surv
}

if(z.count==7){
  perf.mat= z.sum
}

#***************************************
#standardize to relative fitness and centered on trait mean
shift.sample.mean= mean(shift.sample)
breadth.sample.mean= mean(breadth.sample)
fit.mean= mean(perf.mat)

rel.fit= perf.mat/fit.mean
st.shift= shift.sample- shift.sample.mean
st.breadth= breadth.sample- breadth.sample.mean

#UNIVARIATE REGRESSION
##selection analysis
Seln.shift.mod<- lm(rel.fit~st.shift+I(st.shift^2))
Seln.breadth.mod<- lm(rel.fit~st.breadth+I(st.breadth^2))
##  summary(Seln.mod)
Beta.shift <-Seln.shift.mod$coeff[2]
Gamma.shift <- 2*Seln.shift.mod$coeff[3]
Beta.breadth <-Seln.breadth.mod$coeff[2]
Gamma.breadth <- 2*Seln.breadth.mod$coeff[3]

#Response to selection
R2seln.shift <- g[1,1]*(shift.sd^2)*Beta.shift + g[1,2]*(breadth.sd^2)*Beta.breadth
R2seln.breadth <- g[1,2]*(shift.sd^2)*Beta.shift + g[2,2]*(breadth.sd^2)*Beta.breadth
R2seln.shift <- g[1,1]*(shift.sd^2)*Beta.shift 
R2seln.breadth <- g[2,2]*(breadth.sd^2)*Beta.breadth

#update Abs value
shift.mean <- shift.mean + R2seln.shift
breadth.mean <- breadth.mean + R2seln.breadth
#Constain breadth
if(breadth.mean>0.15) breadth.mean<-0.15
if(breadth.mean<0.05) breadth.mean<-0.05
#Constain shift
if(shift.mean>4) shift.mean<-4
if(shift.mean< (-15)) shift.mean<- (-15)

} #end generation loop

  shift.means[k, z.count]= shift.mean
  breadth.means[k, z.count]= breadth.mean  
  
} #end survival loop

} #end loop scenarios 

#-----------
#write out shifts and breadths
shifts= shift.means
breadths= breadth.means

#========================================================================
#PLOTS

## FIGURE
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Extremes/Figs/EvolFigs/")
file<-paste("Q1Impacts_2020.pdf" ,sep="", collapse=NULL)
pdf(file,height = 5, width = 5)

par(mfrow=c(1,1), mar=c(4,4,2,0), oma=c(0,0,0,0), bty="l", lty="solid")

ltys= c("solid","dashed","dotted","solid","dashed","dotted", "solid")
cols= c("black","black","black","darkgrey","darkgrey","darkgrey", "grey")
lwds= c(3,3,3,3,3,3,5)

d <- density(clim.dat$V2)  
plot(d, xlim=range(5,45), main="", ylab="Performance", xlab="Temperature (?C)", col="gray", lty="solid",yaxt="n", cex.lab=1.2)
polygon(d, col="lightgray", border="lightgray")

for(z.count in 1:7){ ##CHANGE TO 3 TO CORRESPOND TO z.count
  
  T=1:60
  z1= TPC.beta(T, shift=shift.means[1,z.count], breadth=breadth.means[1,z.count]) 
  par(new=TRUE)
  if(z.count==1){
    zscale= TPC.beta(T, shift=shift.means[1,z.count], breadth=min(breadth.means, na.rm=TRUE))
    plot(T,z1, type="l", ylim=range(0,max(zscale)), xlim=range(5,45),  ylab="", xlab="", lty=ltys[z.count], lwd=lwds[z.count], col=cols[z.count], axes=FALSE)
  }
  if(z.count>1) points(T,z1, type="l",lty=ltys[z.count], lwd=lwds[z.count], col=cols[z.count])
  
  #axis(side=4)
  
} #end loop plots

legend("topright", legend=c("noncumulative","hardening", "stress response"), lty=c("solid","dashed","dotted"), col=c("black","black","black") , bty="n")
#legend("bottomright", legend=c("injury","mortality"), lty=c("solid"), col=c("black","grey") , bty="n")

dev.off()

