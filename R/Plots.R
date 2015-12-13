# Remove all variables
rm(list=ls())
library(ggplot2)
library(reshape2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
four_color <- c(cbPalette[7],cbPalette[2],cbPalette[4],cbPalette[6])

# Set the working directory
setwd("~/local/Modelling/me2me3_slow/")

"Firing_t_s60ctrl60cc20a1_00b1_00thresh0_33turn0_00100000Rep_st1_repressed.txt"
s <- 60
ctrl <- 60
cc <- 20
a <- 1.0
b <- 1.0
f <- 0.333
turn <- 0.001
rep <- "Rep"
st <- 1
id <- "_active"

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
bstr <- paste('b',gsub("\\.", "_",sprintf("%0.2f",b)),sep="")
fstr <- paste('thresh',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
turnstr <- paste('turn',gsub("\\.", "_",sprintf("%0.8f",turn)),sep="")

time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,id,".txt",sep="")
tDep_me0_file <- paste("me0_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,id,".txt",sep="")
tDep_me1_file <- paste("me1_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,id,".txt",sep="")
tDep_me2_file <- paste("me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,id,".txt",sep="")
tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,id,".txt",sep="")
#tDep_H3_1_file <- paste("H3_1_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,id,".txt",sep="")
#tDep_H3_3_file <- paste("H3_3_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,id,".txt",sep="")
tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,id,".txt",sep="")

space_file <- paste("spatial_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,id,".txt",sep="")


time <- read.table(time_file) ; colnames(time) <- "time"
me0 <- read.table(tDep_me0_file) ; colnames(me0) <- "level"
me1 <- read.table(tDep_me1_file) ; colnames(me1) <- "level"
me2 <- read.table(tDep_me2_file) ; colnames(me2) <- "level"
me3 <- read.table(tDep_me3_file) ; colnames(me3) <- "level"
#H3_1 <- read.table(tDep_H3_1_file) ; colnames(H3_1) <- "level"
#H3_3 <- read.table(tDep_H3_3_file) ; colnames(H3_3) <- "level"
Firing <- read.table(tDep_firing_file)
space <- read.csv(space_file,header=TRUE)

## Plot alpha
# 
# alpha_file <- paste("alpha_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,id,".txt",sep="")
# stochasticAlpha <- read.table(alpha_file,header=TRUE)

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

ggplot(space,aes(x=i,y=pos)) + 
  geom_tile(aes(fill=as.factor(K27))) + 
  scale_fill_manual("H3K27me",values=c(cbPalette[4],cbPalette[5],cbPalette[2],cbPalette[7])) + 
  scale_x_continuous("time") +
  scale_y_continuous("position")
ggsave("ExampleRepressedLocus.pdf",width=10,height=7,units="cm")
