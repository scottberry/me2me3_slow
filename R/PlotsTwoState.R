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
cc <- 10
f <- 1.0
rep <- "NoRep"
st <- 1
id <- ""

fstr <- paste('thresh',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")

time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,fstr,rep,"_st",st,id,".txt",sep="")
tDep_me0_file <- paste("me0_t_s",s,"ctrl",ctrl,"cc",cc,fstr,rep,"_st",st,id,".txt",sep="")
tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,fstr,rep,"_st",st,id,".txt",sep="")
tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,fstr,rep,"_st",st,id,".txt",sep="")

time <- read.table(time_file) ; colnames(time) <- "time"
me0 <- read.table(tDep_me0_file) ; colnames(me0) <- "level"
me3 <- read.table(tDep_me3_file) ; colnames(me3) <- "level"
Firing <- read.table(tDep_firing_file)

par(mfrow=c(3,1),mar=c(2,4,0,0)+0.5,oma=c(3,3,0,0))
plot(time$time/(3600*22),me0$level,type="l",ylim=c(0,1),col="blue3")
#plot(time$time/(3600*22),me1$level,type="l",ylim=c(0,1),col="blue3")
#plot(time$time/(3600*22),me2$level,type="l",ylim=c(0,1),col="red3")
plot(time$time/(3600*22),me3$level,type="l",ylim=c(0,1),col="red3")
#plot(time$time/(3600*22),H3_1$level,type="l",ylim=c(0,1),col="green3")
#plot(time$time/(3600*22),H3_3$level,type="l",ylim=c(0,1),col="orange3")
# plot(alpha ~ time, data=stochasticAlpha,ylim=c(0,max(stochasticAlpha$alpha*1.02)),type="l")
hist(Firing$V1/(3600*22),seq(0,(max(Firing)/(3600*22) + 1/22),by=1/22),ylab="Firing/hour",main="")
mtext("time (cell cycles)",side=1,line=0,outer=TRUE)
mtext("modification level",side=2,line=0,outer=TRUE)
