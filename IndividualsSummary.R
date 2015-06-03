# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")

# Set the working directory
setwd("~/Network/group-share/berrys/Results/me2me3_individuals/startM")

s <- 60
ctrl <- 60
cc <- 10
st <- 1
id <- 1
g_act <- c(0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28,2.56,5.12,10.24)

transcriptsLastTwoCycles <- function(a) {
  astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
  tot <- 0
  for (id in seq(1,10)) {
    tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,"_",id,".txt",sep="")
    Firing <- read.table(tDep_firing_file)
    tot <- tot + sum(Firing>8*3600*16)
  }
  return(tot/10)
}

startM <- data.frame(g_act)
startM$transcripts <- apply(startM, 1, function(x) transcriptsLastTwoCycles(x[1]) )
startM$initial <- "me3"

setwd("~/Network/group-share/berrys/Results/me2me3_individuals/startU")

startU <- data.frame(g_act)
startU$transcripts <- apply(startU, 1, function(x) transcriptsLastTwoCycles(x[1]) )
startU$initial <- "me0"

res <- rbind(startM,startU)

ggplot(res,aes(x=g_act,y=transcripts,col=initial)) + 
  geom_line() + 
  geom_point() + 
  scale_colour_manual(values=four_colour) +
  scale_x_log10() + 
  scale_y_log10("transcripts in last 2 cell cycles") +
  annotation_logticks(short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm"), size=0.25) +
  theme_thesis_multiplanel

setwd("~/Network/group-share/berrys/Results/me2me3_individuals/startM")

# Wrapper for compliled c code to interpolate gillespie algorithm output
interpolateGillespie <- function(tFile,
                                 dFile,
                                 steps=50,
                                 tmpfile=tempfile())
{
  exec <- file.path(".","ConstTimeInterpolate")
  command <- paste(exec,"-t",tFile,"-d",dFile,"-s",steps,"-o",">",tmpfile)
  system(command)
  timeSeries <- read.table(tmpfile,header=FALSE)
  colnames(timeSeries) <- c("time","level")
  return(timeSeries)
}

interpolate_me3 <- function(id,a) {
  astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
  time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,"_",id,".txt",sep="")
  tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,"_",id,".txt",sep="")
  dat <- interpolateGillespie(time_file,tDep_me3_file)
  dat$id <- id
  return(dat)
}

g_act <- c(0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28,2.56,5.12,10.24)
g <- interpolate_me3(1,g_act[5])
for (i in seq(2,10)) {
  g <- rbind(g,interpolate_me3(i,g_act[5]))
}

ggplot(g,aes(x=time/(3600*16),y=level,group=id)) + 
  geom_line(colour="red3",size=1,alpha=0.2) + 
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,5))

# par(mfrow=c(1,1),mar=c(2,4,0,0)+0.5,oma=c(3,3,0,0))
# hist(Firing$V1/3600,seq(0,(max(Firing)/3600 + 1),by=1),ylab="Firing/hour",main="")
# mtext("time (hours)",side=1,line=0,outer=TRUE)
# mtext("modification level",side=2,line=0,outer=TRUE)

# time <- read.table(time_file) ; colnames(time) <- "time"
# me0 <- read.table(tDep_me0_file) ; colnames(me0) <- "level"
# me1 <- read.table(tDep_me1_file) ; colnames(me1) <- "level"
# me2 <- read.table(tDep_me2_file) ; colnames(me2) <- "level"
# me3 <- read.table(tDep_me3_file) ; colnames(me3) <- "level"

#time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,"_",id,".txt",sep="")
#tDep_me0_file <- paste("me0_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,"_",id,".txt",sep="")
#tDep_me1_file <- paste("me1_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,"_",id,".txt",sep="")
#tDep_me2_file <- paste("me2_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,"_",id,".txt",sep="")
#tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,"_",id,".txt",sep="")

0.01*2^seq(1,20)
