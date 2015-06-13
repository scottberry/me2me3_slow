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
  dat$act <- a
  return(dat)
}

g <- interpolate_me3(1,g_act[1])
for (act in g_act) {
  print(act)
   for (i in seq(2,10)) {
     g <- rbind(g,interpolate_me3(i,act))
   }
}

ggplot(g,aes(x=time/(3600*16),y=level,group=id)) + 
  geom_line(colour="red3",size=1,alpha=0.2) + 
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,5)) + facet_grid(act ~ .) + theme_thesis_multiplanel

ggsave(filename = "startM.pdf",width=10,height=10,units="cm")

setwd("~/Network/group-share/berrys/Results/me2me3_individuals/startU")

g <- interpolate_me3(1,g_act[1])
for (act in g_act) {
  print(act)
  for (i in seq(2,10)) {
    g <- rbind(g,interpolate_me3(i,act))
  }
}

ggplot(g,aes(x=time/(3600*16),y=level,group=id)) + 
  geom_line(colour="red3",size=1,alpha=0.2) + 
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,5)) + facet_grid(act ~ .) + theme_thesis_multiplanel

ggsave(filename = "startU.pdf",width=10,height=10,units="cm")
