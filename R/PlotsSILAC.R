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
cc <- 15
a <- 1.0
st <- 1

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
 
time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me2_file <- paste("me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me2_LIGHT_file <- paste("LIGHT_me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me2_HEAVY_file <- paste("HEAVY_me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me2_UNLABELLED_file <- paste("UNLABELLED_me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me3_LIGHT_file <- paste("LIGHT_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me3_HEAVY_file <- paste("HEAVY_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me3_UNLABELLED_file <- paste("UNLABELLED_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")

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

silacAbs_file <- paste("SilacAbs_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
silacRel_file <- paste("SilacRel_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")

silacAbs <- read.table(silacAbs_file,header = TRUE)
silacRel <- read.table(silacRel_file,header = TRUE)

m <- aggregate(silacRel$level,
               by=list(time=silacRel$time,
                       mod=silacRel$mod,
                       label=silacRel$label),
               FUN=mean, na.rm=TRUE)
colnames(m) <- c("time","mod","label","mean")
stdev <- aggregate(silacRel$level,
                   by=list(time=silacRel$time,
                           mod=silacRel$mod,
                           label=silacRel$label),
                   FUN=sd, na.rm=TRUE)
colnames(stdev) <- c("time","mod","label","sd")

silacRel_summary <- cbind(m,sd=stdev$sd)

K27me3_endcycle <- silacRel_summary[silacRel_summary$label=="LIGHT" & 
                                      silacRel_summary$mod=="me3" &
                                      silacRel_summary$time==min(silacRel_summary$time),]$mean

silacRel_summary_K27me3 <- subset(silacRel_summary,mod=="me3",select=-c(mod))
silacRel_summary_K27me3$mean <- silacRel_summary_K27me3$mean*0.3010592/K27me3_endcycle

ggplot(subset(silacAbs,mod=="me3"),aes(x=time,y=level,col=label)) + geom_point()
ggplot(subset(silacRel,mod=="me3"),aes(x=time,y=level,col=label)) + geom_point()

ggplot(silacRel_summary_K27me3,aes(x=time,y=mean,col=label)) + geom_line() +
  scale_y_continuous(limits=c(0,0.5))

K27me3_Alabert <- readRDS(file = "R/K27me3_Time_Alabert.rds")
K27me3_Alabert <- subset(K27me3_Alabert,select=-c(Modification))
K27me3_Alabert$Time <- as.numeric(as.character(K27me3_Alabert$Time))

## Adjust summary to fit experimental data format
colnames(silacRel_summary_K27me3) <- c("Time","Label","Mean_Percent","Sd_Percent")
silacRel_summary_K27me3$Label <- revalue(silacRel_summary_K27me3$Label,c("HEAVY"="H","LIGHT"="L"))
silacRel_summary_K27me3$Time <- (silacRel_summary_K27me3$Time - min(silacRel_summary$time))/3600
silacRel_summary_K27me3

ggplot(K27me3_Alabert,aes(x=Time,y=Mean_Percent,col=Label)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=Mean_Percent-Sd_Percent,ymax=Mean_Percent+Sd_Percent),width=0.1) +
  geom_line(data=silacRel_summary_K27me3)
