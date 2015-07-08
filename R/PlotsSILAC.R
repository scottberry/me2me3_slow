# Remove all variables
rm(list=ls())
library(ggplot2)
library(plyr)
library(reshape2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
four_color <- c(cbPalette[7],cbPalette[2],cbPalette[4],cbPalette[6])

# Set the working directory
setwd("~/local/Modelling/me2me3_slow/")

s <- 60
ctrl <- 60
cc <- 50
a <- 0.0
b <- 1.0
f <- 0.4
tau <- 4.0
st <- 10

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
bstr <- paste('b',gsub("\\.", "_",sprintf("%0.2f",b)),sep="")
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
taustr <- paste('tau',gsub("\\.", "_",sprintf("%0.2f",tau)),sep="")

time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me2_file <- paste("me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me2_LIGHT_file <- paste("LIGHT_me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me2_HEAVY_file <- paste("HEAVY_me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me2_UNLABELLED_file <- paste("UNLABELLED_me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me3_LIGHT_file <- paste("LIGHT_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me3_HEAVY_file <- paste("HEAVY_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me3_UNLABELLED_file <- paste("UNLABELLED_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")

time <- t(read.table(time_file))
me2 <- t(read.table(tDep_me2_file))
me3 <- t(read.table(tDep_me3_file))
Firing <- t(read.table(tDep_firing_file))
me2_LIGHT <- t(read.table(tDep_me2_LIGHT_file))
me2_HEAVY <- t(read.table(tDep_me2_HEAVY_file))
me2_UNLABELLED <- t(read.table(tDep_me2_UNLABELLED_file))
me3_LIGHT <- t(read.table(tDep_me3_LIGHT_file))
me3_HEAVY <- t(read.table(tDep_me3_HEAVY_file))
me3_UNLABELLED <- t(read.table(tDep_me3_UNLABELLED_file))

par(mfrow=c(3,1),mar=c(2,4,0,0)+0.5,oma=c(3,3,0,0))
plot(time/3600,me2,type="l",ylim=c(0,1),col="green3")
lines(time/3600,me2_LIGHT,type="l",ylim=c(0,1),col="orange")
lines(time/3600,me2_HEAVY,type="l",ylim=c(0,1),col="blue3")
lines(time/3600,me2_UNLABELLED,type="l",ylim=c(0,1),col="grey60")
plot(time/3600,me3,type="l",ylim=c(0,1),col="green3")
lines(time/3600,me3_LIGHT,type="l",ylim=c(0,1),col="orange")
lines(time/3600,me3_HEAVY,type="l",ylim=c(0,1),col="blue3")
lines(time/3600,me3_UNLABELLED,type="l",ylim=c(0,1),col="grey60")
hist(Firing/3600,breaks=seq(0,max(Firing/3600+1),by=0.5),freq=TRUE,main="",ylab="Firing events")
mtext("time (hours)",side=1,line=0,outer=TRUE)
