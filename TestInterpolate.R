# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")

# Set the working directory
setwd("~/Network/group-share/berrys/me2me3_slow/")

s <- 60
ctrl <- 60
cc <- 10
st <- 1
id <- 1
a <- 1.0

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")

time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
tDep_me1_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
inter_me1_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,"_int.txt",sep="")

time <- read.table(time_file) ; colnames(time) <- "time"
me3 <- read.table(tDep_me1_file) ; colnames(me3) <- "level"
int <- read.table(inter_me1_file) ; colnames(int) <- c("time","level")

raw <- data.frame(time,me3)
raw$d <- "gillespie"
int$d <- "interpolated"

all <- rbind(raw,int)

ggplot(all,aes(x=time/(3600*16),y=level,col=d)) + geom_point(data=raw) + geom_line(data=int) +
  scale_x_continuous(limits=c(0,4))
