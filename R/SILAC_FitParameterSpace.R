# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")
library(grid)
library(gtable)
library(DEoptim)

# Set the working directory
setwd("~/local/Modelling/me2me3_slow/")

## Note: Run simulations twice:
# 1. with -m for Silac results
# 2. with -i bal for parameter space results

s <- 60
ctrl <- 60
cc <- 15
a <- 1.0
tau <- 4.0
st <- 20

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
taustr <- paste('tau',gsub("\\.", "_",sprintf("%0.2f",tau)),sep="")

f <- 0.2
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_0_2 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_2 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,".txt",sep="")
f <- 0.4
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_0_4 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_4 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,".txt",sep="")
f <- 0.6
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_0_6 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_6 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,".txt",sep="")
f <- 0.8
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_0_8 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_8 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,".txt",sep="")
f <- 1.0
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_1_0 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_1_0 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,fstr,taustr,"st",st,".txt",sep="")


silacAvg_0_2 <- read.table(silacAvg_file_0_2,header = TRUE)
parameterSpace_0_2 <- read.table(parameterSpace_file_0_2,header = TRUE)
silacAvg_0_4 <- read.table(silacAvg_file_0_4,header = TRUE)
parameterSpace_0_4 <- read.table(parameterSpace_file_0_4,header = TRUE)
silacAvg_0_6 <- read.table(silacAvg_file_0_6,header = TRUE)
parameterSpace_0_6 <- read.table(parameterSpace_file_0_6,header = TRUE)
silacAvg_0_8 <- read.table(silacAvg_file_0_8,header = TRUE)
parameterSpace_0_8 <- read.table(parameterSpace_file_0_8,header = TRUE)
silacAvg_1_0 <- read.table(silacAvg_file_1_0,header = TRUE)
parameterSpace_1_0 <- read.table(parameterSpace_file_1_0,header = TRUE)

silacAvg <- rbind(silacAvg_0_2,silacAvg_0_4,silacAvg_0_6,silacAvg_0_8,silacAvg_1_0)
parameterSpace <- rbind(parameterSpace_0_2,parameterSpace_0_4,parameterSpace_0_6,parameterSpace_0_8,parameterSpace_1_0)

## Extract the cell-cycle-end K27me3 level for this parameter set (over all cell cycles)
endVal <- subset(parameterSpace,select=c(FIRING,P_DEMETHYLATE,P_METHYLATE,FIRING_THRESHOLD,me3_end))

## Adjust relative SILAC values by this asymptotic cell-cycle-end value
adjustedSILAC <- merge(silacAvg,endVal,by=c("FIRING","P_DEMETHYLATE","P_METHYLATE","FIRING_THRESHOLD"))
adjustedSILAC$adjusted_level <- adjustedSILAC$level/adjustedSILAC$me3_end

## Read experimental data file and adjust names
K27me3_Alabert <- readRDS(file = "R/K27me3_SILAC_expt.rds")
K27me3_Alabert$Time <- as.numeric(as.character(K27me3_Alabert$Time))
colnames(K27me3_Alabert) <- c("mod","time","label","level")
K27me3_Alabert$mod <- revalue(K27me3_Alabert$mod,c("K27me3"="me3"))
K27me3_Alabert$label <- revalue(K27me3_Alabert$label,c("H"="HEAVY","L"="LIGHT"))

## Extract the cell-cyle-end K27me3 from the experimental data and use this to 
## normalise to the experimental data
me3_light_0 <- mean(K27me3_Alabert[K27me3_Alabert$label=="LIGHT" & K27me3_Alabert$time==0,]$level)
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

## Add dummy data to experimental data frame for plotting over all parameter space
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
pars <- expand.grid(P_DEMETHYLATE=unique(adjustedSILAC$P_DEMETHYLATE),
                    P_METHYLATE=unique(adjustedSILAC$P_METHYLATE))
K27me3_Alabert_Dummy <- expand.grid.df(pars,K27me3_Alabert_summary)

ggplot(K27me3_Alabert_Dummy,aes(x=time,y=adjusted_level,col=label,group=label)) + 
  geom_point() + 
  geom_errorbar(width=3,aes(ymin=adjusted_level-adjusted_sd,
                    ymax=adjusted_level+adjusted_sd)) +
  facet_grid(P_METHYLATE ~ P_DEMETHYLATE) + 
  geom_line(data=subset(adjustedSILAC,FIRING_THRESHOLD==1.0)) + 
  scale_y_continuous(limits=c(0,2)) + 
  scale_color_manual(values=c("blue2","orange")) +
  theme_thesis_multiplanel

ggsave(file="ParameterSpaceFit.pdf",width=20,height=15,units="cm")
# 

# 
# ## me3_end
# p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
# p <- p + geom_raster(aes(fill=me3_end)) + 
#   scale_y_log10("Methylation rate (per second)",
#                 labels = trans_format("log10", math_format(10^.x))) +
#   scale_x_log10("Demethylation probability (per histone per firing)",
#                 labels = trans_format("log10", math_format(10^.x))) +
#   coord_fixed(ratio=1) + facet_grid(FIRING_K ~ .) +
#   scale_fill_gradientn(name="me3_end",colours=rev(rainbow(3)),limits=c(0,1)) +
#   theme_bw(7) + theme_thesis_multiplanel +
#   theme(plot.title = element_text(lineheight=.8, face="bold"),
#         panel.margin=unit(0,"lines"),
#         plot.margin = unit(c(0,0,0,0), "cm"),
#         axis.text.x = element_text(angle=0),
#         strip.text.x = element_text(size=7, angle=0),
#         strip.text.y = element_text(size=7, angle=0))
# p

## bistability (15 cell cycles)
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_raster(aes(fill=probU)) + 
  scale_y_log10("Methylation rate (per second)",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability (per histone per firing)",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + facet_grid(FIRING_THRESHOLD ~ .) +
  scale_fill_gradientn(name="Bistability",colours=rev(rainbow(3)),limits=c(0,1)) +
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        strip.text.x = element_text(size=7, angle=0),
        strip.text.y = element_text(size=7, angle=0))
p
ggsave("Bistability.pdf",width=7,height=18,units="cm")

## Quantify model fit to data
summary(K27me3_Alabert)
summary(adjustedSILAC)

comparison <- merge(subset(K27me3_Alabert,select=c(mod,time,label,adjusted_level)),
      subset(adjustedSILAC,select=c(mod,time,label,adjusted_level,P_METHYLATE,P_DEMETHYLATE,FIRING_THRESHOLD)),
      by=c("mod","time","label"))

names(comparison)
comparison$sq_err <- (comparison$adjusted_level.x - comparison$adjusted_level.y)^2
sum_sq_err <- aggregate(comparison$sq_err,by=list(comparison$P_METHYLATE,
                                                  comparison$P_DEMETHYLATE,
                                                  comparison$FIRING_THRESHOLD),
          FUN = sum)
colnames(sum_sq_err) <- c("P_METHYLATE","P_DEMETHYLATE","FIRING_THRESHOLD","SSE")

p <- ggplot(sum_sq_err,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_raster(aes(fill=SSE)) + 
  scale_y_log10("Methylation rate (per second)",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability (per histone per firing)",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + facet_grid(FIRING_THRESHOLD ~ .) +
  scale_fill_gradientn(name="SSE",colours=rev(rainbow(3)),trans="log") +
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        strip.text.x = element_text(size=7, angle=0),
        strip.text.y = element_text(size=7, angle=0))
p
ggsave("SSE.pdf",width=7,height=18,units="cm")

sum_sq_err[sum_sq_err$SSE==min(sum_sq_err$SSE),]

# geom_rect(data = tp,aes(fill = day),xmin = -Inf,xmax = Inf,
#           ymin = -Inf,ymax = Inf,alpha = 0.3)
