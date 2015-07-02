# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")

# Set the working directory
setwd("~/local/Modelling/me2me3_slow/")

s <- 60
ctrl <- 60
cc <- 20
a <- 1.0
f <- 1.0
tau <- 4.0
st <- 1

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
taustr <- paste('tau',gsub("\\.", "_",sprintf("%0.2f",tau)),sep="")

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

interpolate_me3 <- function(id,steps=50) {
  astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
  time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,"_",id,".txt",sep="")
  tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,"_",id,".txt",sep="")
  tDep_me3_LIGHT_file <- paste("LIGHT_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,"_",id,".txt",sep="")
  tDep_me3_HEAVY_file <- paste("HEAVY_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,"_",id,".txt",sep="")
  tDep_me3_UNLABELLED_file <- paste("UNLABELLED_me3_t_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,"_",id,".txt",sep="")
  total <- interpolateGillespie(time_file,tDep_me3_file,steps=steps)
  total$species <- "Total"
  light <- interpolateGillespie(time_file,tDep_me3_LIGHT_file,steps=steps)
  light$species <- "Old"
  heavy <- interpolateGillespie(time_file,tDep_me3_HEAVY_file,steps=steps)
  heavy$species <- "New"
  unlabelled <- interpolateGillespie(time_file,tDep_me3_UNLABELLED_file,steps=steps)
  unlabelled$species <- "Unlabelled"
  dat <- rbind(total,light,heavy,unlabelled)
  dat$id <- id
  dat$species <- factor(dat$species)
  return(dat)
}

g <- interpolate_me3(1,100)
for (i in seq(2,20)) {
  g <- rbind(g,interpolate_me3(i,100))
}

g$species <- factor(g$species,levels=c("Old","New","Unlabelled","Total"),ordered = TRUE)

ggplot(data=g,aes(x=time/(3600*22),y=level,col=species)) + 
  geom_line(aes(group=id),size=0.4,alpha=0.2) + 
  facet_grid(species ~ .) + 
  scale_color_manual(values=c("orange","blue2","grey60","darkgreen")) +
  scale_y_continuous(name="K27me3 (proportion of total histone)",limits=c(0,1.02),breaks=c(0,1)) +
  scale_x_continuous(name="Time (cell cycles)",limits=c(0,13),breaks=c(0,6,12)) + 
  theme_thesis_multiplanel + theme(legend.position="none")

ggsave("Threshold1_0_SILAC_timecourses.png",width=8,height=6,units="cm",dpi=600)

## Replot as a function of labelled histone using same method as earlier.
## Also look at average firing rate through the cell cycle to determine why there is a "dip"?
g_sub <- subset(g,species %in% c("Old","New") & time > 22*6*3600 & time < 22*9*3600)
g_sub$time <- round(g_sub$time/3600)-(6*22)
agg <- aggregate(g_sub$level,by=list(g_sub$species,
                                     g_sub$time),
          FUN=mean)
colnames(agg) <- c("species","time","mean")
ggplot(data=agg,aes(x=time,y=mean,col=species)) + geom_line() + theme_thesis_multiplanel


