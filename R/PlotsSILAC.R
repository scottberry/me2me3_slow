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
cc <- 6
a <- 1.0
st <- 1

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")

time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me0_file <- paste("me0_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me1_file <- paste("me1_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me2_file <- paste("me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")

tDep_me0_LIGHT_file <- paste("LIGHT_me0_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me1_LIGHT_file <- paste("LIGHT_me1_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me2_LIGHT_file <- paste("LIGHT_me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me3_LIGHT_file <- paste("LIGHT_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")

tDep_me0_HEAVY_file <- paste("HEAVY_me0_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me1_HEAVY_file <- paste("HEAVY_me1_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me2_HEAVY_file <- paste("HEAVY_me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me3_HEAVY_file <- paste("HEAVY_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")

time <- t(read.table(time_file))
me0 <- t(read.table(tDep_me0_file))
me1 <- t(read.table(tDep_me1_file))
me2 <- t(read.table(tDep_me2_file))
me3 <- t(read.table(tDep_me3_file))
Firing <- t(read.table(tDep_firing_file))

me0_OLD <- t(read.table(tDep_me0_LIGHT_file))
me1_OLD <- t(read.table(tDep_me1_LIGHT_file))
me2_OLD <- t(read.table(tDep_me2_LIGHT_file))
me3_OLD <- t(read.table(tDep_me3_LIGHT_file))

me0_NEW <- t(read.table(tDep_me0_HEAVY_file))
me1_NEW <- t(read.table(tDep_me1_HEAVY_file))
me2_NEW <- t(read.table(tDep_me2_HEAVY_file))
me3_NEW <- t(read.table(tDep_me3_HEAVY_file))

par(mfrow=c(5,1),mar=c(2,4,0,0)+0.5,oma=c(3,3,0,0))
plot(time/3600,me0,type="l",ylim=c(0,1),col="green3")
lines(time/3600,me0_OLD,type="l",ylim=c(0,1),col="orange")
lines(time/3600,me0_NEW,type="l",ylim=c(0,1),col="blue3")
plot(time/3600,me1,type="l",ylim=c(0,1),col="green3")
lines(time/3600,me1_OLD,type="l",ylim=c(0,1),col="orange")
lines(time/3600,me1_NEW,type="l",ylim=c(0,1),col="blue3")
plot(time/3600,me2,type="l",ylim=c(0,1),col="green3")
lines(time/3600,me2_OLD,type="l",ylim=c(0,1),col="orange")
lines(time/3600,me2_NEW,type="l",ylim=c(0,1),col="blue3")
plot(time/3600,me3,type="l",ylim=c(0,1),col="green3")
lines(time/3600,me3_OLD,type="l",ylim=c(0,1),col="orange")
lines(time/3600,me3_NEW,type="l",ylim=c(0,1),col="blue3")
hist(Firing/3600,breaks=seq(0,max(Firing/3600+1),by=0.5),freq=TRUE,main="",ylab="Firing events")
mtext("time (hours)",side=1,line=0,outer=TRUE)
mtext("modification level",side=2,line=0,outer=TRUE)
# 
# # plot modifications as a fraction of OLD and NEW histones
# tot_OLD = me0_OLD + me1_OLD + me2_OLD + me3_OLD
# tot_OLD[tot_OLD==0] <- NA
# me0_OLD_frac = me0_OLD/tot_OLD
# me1_OLD_frac = me1_OLD/tot_OLD
# me2_OLD_frac = me2_OLD/tot_OLD
# me3_OLD_frac = me3_OLD/tot_OLD
# 
# tot_NEW = me0_NEW + me1_NEW + me2_NEW + me3_NEW
# tot_NEW[tot_NEW==0] <- NA
# me0_NEW_frac = me0_NEW/tot_NEW
# me1_NEW_frac = me1_NEW/tot_NEW
# me2_NEW_frac = me2_NEW/tot_NEW
# me3_NEW_frac = me3_NEW/tot_NEW
# 
# par(mfrow=c(5,1),mar=c(2,4,0,0)+0.5,oma=c(3,3,0,0))
# plot(time/3600,me0_OLD_frac,type="l",ylim=c(0,1),col="orange",xlim=c(0,24*days),ylab="me0")
# lines(time/3600,me0_NEW_frac,type="l",ylim=c(0,1),col="blue3")
# plot(time/3600,me1_OLD_frac,type="l",ylim=c(0,1),col="orange",xlim=c(0,24*days),ylab="me1")
# lines(time/3600,me1_NEW_frac,type="l",ylim=c(0,1),col="blue3")
# plot(time/3600,me2_OLD_frac,type="l",ylim=c(0,1),col="orange",xlim=c(0,24*days),ylab="me2")
# lines(time/3600,me2_NEW_frac,type="l",ylim=c(0,1),col="blue3")
# plot(time/3600,me3_OLD_frac,type="l",ylim=c(0,1),col="orange",xlim=c(0,24*days),ylab="me3")
# lines(time/3600,me3_NEW_frac,type="l",ylim=c(0,1),col="blue3")
# hist(Firing/3600,breaks=seq(0,max(Firing/3600+1),by=0.5),freq=TRUE,xlim=c(0,24*days),main="",ylab="Firing events")
# mtext("time (hours)",side=1,line=0,outer=TRUE)
# mtext("modification level (fraction of total)",side=2,line=0,outer=TRUE)

# par(mfrow=c(1,1),mar=c(2,4,0,0)+0.5,oma=c(3,3,0,0))
# plot(time/3600,me3_OLD_frac,type="l",ylim=c(0,1),col="orange",xlim=c(0,24*days),ylab="me3")
# lines(time/3600,me3_NEW_frac,type="l",ylim=c(0,1),col="blue3")