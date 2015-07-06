# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")

library(gridExtra)

# Set the working directory
setwd("~/local/Modelling/me2me3_slow/")

s <- 60
ctrl <- 60
cc <- 20
a <- 1.0
f <- 0.3
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

g <- interpolate_me3(1,200)
for (i in seq(2,20)) {
  g <- rbind(g,interpolate_me3(i,200))
}

g$species <- factor(g$species,levels=c("Old","New","Unlabelled","Total"),ordered = TRUE)

ggplot(data=g,aes(x=time/(3600*22),y=level,col=species)) + 
  geom_line(aes(group=id),size=0.4,alpha=0.2) + 
  facet_grid(species ~ .) + 
  scale_color_manual(values=c("orange","blue2","grey60","darkgreen")) +
  scale_y_continuous(name="K27me3 (proportion of total histone)",limits=c(0,1.02),breaks=c(0,1)) +
  scale_x_continuous(name="Time (cell cycles)",limits=c(0,13),breaks=c(0,6,12)) + 
  theme_thesis_multiplanel + theme(legend.position="none")

ggsave("Threshold0_3_SILAC_timecourses.png",width=8,height=6,units="cm",dpi=600)

## Average time-courses around SILAC experiment with experimental data

g_sub <- subset(g,species %in% c("Old","New") & time > 22*6*3600 & time < 22*9*3600)
g_sub$time <- round(g_sub$time/3600)-(6*22)
agg <- aggregate(g_sub$level,by=list(g_sub$species,
                                     g_sub$time),
                 FUN=mean)
colnames(agg) <- c("species","time","mean")

## Extract cell-cycle end mean from experimental data
ends <- c(22*3600*5:10)
sim_end_val <- mean(subset(g,time %in% ends & species=="Total")$level)

## Read experimental data file and adjust names
K27me3_Alabert <- readRDS(file = "R/K27me3_SILAC_expt.rds")
K27me3_Alabert$Time <- as.numeric(as.character(K27me3_Alabert$Time))
colnames(K27me3_Alabert) <- c("mod","time","label","level")
K27me3_Alabert$mod <- revalue(K27me3_Alabert$mod,c("K27me3"="me3"))
K27me3_Alabert$label <- revalue(K27me3_Alabert$label,c("H"="New","L"="Old"))

## Extract the cell-cyle-end K27me3 from the experimental data and use this to 
## normalise to the experimental data
me3_light_0 <- mean(K27me3_Alabert[K27me3_Alabert$label=="Old" & K27me3_Alabert$time==0,]$level)
K27me3_Alabert$adjusted_level <- K27me3_Alabert$level/me3_light_0

## Summarise the experimental data by mean and sd
K27me3_Alabert_mean <- aggregate(K27me3_Alabert$adjusted_level,
                                 by=list(K27me3_Alabert$mod,
                                         K27me3_Alabert$time,
                                         K27me3_Alabert$label),
                                 FUN=mean)
colnames(K27me3_Alabert_mean) <- c("mod","time","label","adjusted_level")

K27me3_Alabert_sd <- aggregate(K27me3_Alabert$adjusted_level,
                               by=list(K27me3_Alabert$mod,
                                       K27me3_Alabert$time,
                                       K27me3_Alabert$label),
                               FUN=sd)
colnames(K27me3_Alabert_sd) <- c("mod","time","label","adjusted_sd")

K27me3_Alabert_summary <- cbind(K27me3_Alabert_mean,
                                adjusted_sd=K27me3_Alabert_sd$adjusted_sd)

# Re-adjust experimental data to simulation cell-cycle end value and also to show dilution 
K27me3_Alabert_summary$adjusted_level <- K27me3_Alabert_summary$adjusted_level*sim_end_val
K27me3_Alabert_summary$adjusted_sd <- K27me3_Alabert_summary$adjusted_sd*sim_end_val

K27me3_Alabert_summary[K27me3_Alabert_summary$time==24,]$adjusted_level <- K27me3_Alabert_summary[K27me3_Alabert_summary$time==24,]$adjusted_level/2
K27me3_Alabert_summary[K27me3_Alabert_summary$time==24,]$adjusted_sd <- K27me3_Alabert_summary[K27me3_Alabert_summary$time==24,]$adjusted_sd/2
K27me3_Alabert_summary[K27me3_Alabert_summary$time==48,]$adjusted_level <- K27me3_Alabert_summary[K27me3_Alabert_summary$time==48,]$adjusted_level/4
K27me3_Alabert_summary[K27me3_Alabert_summary$time==48,]$adjusted_sd <- K27me3_Alabert_summary[K27me3_Alabert_summary$time==48,]$adjusted_sd/4
 
colnames(K27me3_Alabert_summary) <- c("mod","time","species","mean","sd")

## Plot simulation with experimental data
p_silac_avg <- ggplot(data=agg,aes(x=time,y=mean,col=species)) + 
  geom_line() + 
  scale_color_manual(values = c("orange","blue2")) + 
  scale_y_continuous("K27me3 (proportion of total)") +
  theme_thesis_multiplanel + 
  geom_point(data=K27me3_Alabert_summary,size=1) + 
  geom_errorbar(data=K27me3_Alabert_summary,aes(ymin=mean-sd,ymax=mean+sd),width=2,size=0.5) +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.5,-0.3,0.5), "lines"))

id <- 1
binwidth <- 1800
loci <- 20
brks <- seq(0,22*10*3600,by=binwidth)
astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,"_",id,".txt",sep="")
Firing <- t(read.table(tDep_firing_file))
Firing <- Firing[Firing<22*10*3600]
binned <- findInterval(Firing,brks)

for (id in seq(2,loci)) {
  tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,"_",id,".txt",sep="")
  Firing <- t(read.table(tDep_firing_file))
  Firing <- Firing[Firing<22*10*3600]
  binned <- c(binned,findInterval(Firing,brks))
}

firingCounts <- as.data.frame(brks[1:length(brks)-1])
firingCounts <- cbind(firingCounts,table(binned))
colnames(firingCounts) <- c("t","bin","count")
firingCounts$t <- firingCounts$t + binwidth/2
firingCounts$count <- firingCounts$count/loci
firingCounts <- subset(firingCounts,t > 22*6*3600 & t< (22*6+60)*3600)
firingCounts$t <- firingCounts$t/3600 - 6*22

p_firing <- ggplot(data=firingCounts, aes(x=t, y=count)) +
  geom_bar(stat="identity",fill="#149F49",width=binwidth/(3600)) + theme_thesis_multiplanel + 
  scale_x_continuous(name="Time (hours)") + scale_y_continuous(name="Firing",breaks=c(0,2,4)) +
  theme(plot.margin = unit(c(-0.3,0.5,0.5,0.5), "lines"))

gp1<- ggplotGrob(p_silac_avg)
gp2<- ggplotGrob(p_firing)
gp2$widths <- gp1$widths

pdf(file="test.pdf",width=7*0.393700787,height=6*0.393700787)
grid.arrange(gp1,gp2,ncol=1, heights=c(1,0.5))
dev.off()


