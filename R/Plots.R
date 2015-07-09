# Remove all variables
rm(list=ls())
library(ggplot2)
library(reshape2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
four_color <- c(cbPalette[7],cbPalette[2],cbPalette[4],cbPalette[6])

# Set the working directory
setwd("~/local/Modelling/me2me3_slow/")

s <- 60
ctrl <- 60
cc <- 20
a <- 1.0
b <- 1.0
f <- 1.0
tau <- 0.0
st <- 1

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
bstr <- paste('b',gsub("\\.", "_",sprintf("%0.2f",b)),sep="")
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
taustr <- paste('tau',gsub("\\.", "_",sprintf("%0.2f",tau)),sep="")

time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me0_file <- paste("me0_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me1_file <- paste("me1_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me2_file <- paste("me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")

time <- read.table(time_file) ; colnames(time) <- "time"
me0 <- read.table(tDep_me0_file) ; colnames(me0) <- "level"
me1 <- read.table(tDep_me1_file) ; colnames(me1) <- "level"
me2 <- read.table(tDep_me2_file) ; colnames(me2) <- "level"
me3 <- read.table(tDep_me3_file) ; colnames(me3) <- "level"
Firing <- read.table(tDep_firing_file)

par(mfrow=c(5,1),mar=c(2,4,0,0)+0.5,oma=c(3,3,0,0))
plot(time$time/(3600*22),me0$level,type="l",ylim=c(0,1),col="blue3")
plot(time$time/(3600*22),me1$level,type="l",ylim=c(0,1),col="blue3")
plot(time$time/(3600*22),me2$level,type="l",ylim=c(0,1),col="red3")
plot(time$time/(3600*22),me3$level,type="l",ylim=c(0,1),col="red3")
hist(Firing$V1/(3600*22),seq(0,(max(Firing)/(3600*22) + 1/22),by=1/22),ylab="Firing/hour",main="")
mtext("time (cell cycles)",side=1,line=0,outer=TRUE)
mtext("modification level",side=2,line=0,outer=TRUE)

