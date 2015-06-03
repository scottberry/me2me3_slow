# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")

# Set the working directory
setwd("~/Network/group-share/berrys/me2me3_slow/")

s <- 60
ctrl <- 60
cc <- 100 
st <- 1

# Read data
time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")
tDep_me0_file <- paste("me0_t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")
tDep_me1_file <- paste("me1_t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")
tDep_me2_file <- paste("me2_t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")
tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")
tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")

time <- read.table(time_file) ; colnames(time) <- "time"
me0 <- read.table(tDep_me0_file) ; colnames(me0) <- "level"
me1 <- read.table(tDep_me1_file) ; colnames(me1) <- "level"
me2 <- read.table(tDep_me2_file) ; colnames(me2) <- "level"
me3 <- read.table(tDep_me3_file) ; colnames(me3) <- "level"
Firing <- read.table(tDep_firing_file)

# Extract each cell cycle and interpolate "constant" for equal 
# time intervals in order to average

library(scales)
# First cell cycle
par(mfrow=c(1,1),mar=c(2,4,0,0)+0.5,oma=c(3,3,0,0))
t <- time[time$time>=16*3600 & time$time<=16*2*3600,]
sub_me3 <- me3[time$time>=16*3600 & time$time<=16*2*3600,]
plot((t %% (16*3600))/3600,sub_me3,type="l",ylim=c(0,1),col=alpha("red", 0.2),lwd=5)

# Interpolation
t_val <- seq(0.25,15.75,by=0.25)
accum <- as.data.frame(approx(x=(t %% (16*3600))/3600,
       y=sub_me3,
       xout=t_val,
       method="constant",f=0))

# Repeat for remaining cell cycles
for (i in 2:(cc-1)) {
  t <- time[time$time>=16*i*3600 & time$time<=16*(i+1)*3600,]
  sub_me3 <- me3[time$time>=16*i*3600 & time$time<=16*(i+1)*3600,]
  lines((t %% (16*3600))/3600,sub_me3,col=alpha("red", 0.1),lwd=3)
  inter <- approx(x=(t %% (16*3600))/3600,
                   y=sub_me3,
                   xout=t_val,
                   method="constant",f=0)
  accum <- cbind(accum,inter$y)
}

# count non-NA rows
n <- rowSums(!is.na(accum[,2:50]))
# replace NA values
accum[is.na(accum)] <- 0
# average over all non-NA rows
avg <- rowSums(accum[,2:50])/n

# Add average line to plot
lines(t_val,avg,col="black",lwd=3)
