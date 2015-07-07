# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")
library(grid)
library(gtable)

# Set the working directory
setwd("~/Network/group-share/berrys/me2me3_slow/")

## Note: Run simulations twice:
# 1. with -m -t $threshold for Silac results
# 2. with -i bal -t $threshold for parameter space results

s <- 60
ctrl <- 60
cc <- 20
a <- 0.0
b <- 1.0
tau <- 4.0
st <- 22

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
bstr <- paste('b',gsub("\\.", "_",sprintf("%0.2f",b)),sep="")
taustr <- paste('tau',gsub("\\.", "_",sprintf("%0.2f",tau)),sep="")

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

scientific_labeller <- function(variable,value) {
  return(scientific_10(signif(value,2)))
}

scientific_labeller_factor <- function(variable,value) {
  return(scientific_10(signif(as.numeric(as.character(value)),2)))
}

f <- 0.2
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_0_2 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_2 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_2_bal <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,"_bal.txt",sep="")
f <- 0.3
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_0_3 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_3 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_3_bal <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,"_bal.txt",sep="")
f <- 0.4
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_0_4 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_4 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_4_bal <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,"_bal.txt",sep="")
f <- 0.6
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_0_6 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_6 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_6_bal <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,"_bal.txt",sep="")
f <- 0.8
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_0_8 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_8 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_0_8_bal <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,"_bal.txt",sep="")
f <- 1.0
fstr <- paste('fir',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
silacAvg_file_1_0 <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_1_0 <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,".txt",sep="")
parameterSpace_file_1_0_bal <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,taustr,"st",st,"_bal.txt",sep="")


silacAvg_0_2 <- read.table(silacAvg_file_0_2,header = TRUE)
parameterSpace_0_2 <- read.table(parameterSpace_file_0_2,header = TRUE)
parameterSpace_0_2_bal <- read.table(parameterSpace_file_0_2_bal,header = TRUE)
silacAvg_0_3 <- read.table(silacAvg_file_0_3,header = TRUE)
parameterSpace_0_3 <- read.table(parameterSpace_file_0_3,header = TRUE)
parameterSpace_0_3_bal <- read.table(parameterSpace_file_0_3_bal,header = TRUE)
silacAvg_0_4 <- read.table(silacAvg_file_0_4,header = TRUE)
parameterSpace_0_4 <- read.table(parameterSpace_file_0_4,header = TRUE)
parameterSpace_0_4_bal <- read.table(parameterSpace_file_0_4_bal,header = TRUE)
silacAvg_0_6 <- read.table(silacAvg_file_0_6,header = TRUE)
parameterSpace_0_6 <- read.table(parameterSpace_file_0_6,header = TRUE)
parameterSpace_0_6_bal <- read.table(parameterSpace_file_0_6_bal,header = TRUE)
silacAvg_0_8 <- read.table(silacAvg_file_0_8,header = TRUE)
parameterSpace_0_8 <- read.table(parameterSpace_file_0_8,header = TRUE)
parameterSpace_0_8_bal <- read.table(parameterSpace_file_0_8_bal,header = TRUE)
silacAvg_1_0 <- read.table(silacAvg_file_1_0,header = TRUE)
parameterSpace_1_0 <- read.table(parameterSpace_file_1_0,header = TRUE)
parameterSpace_1_0_bal <- read.table(parameterSpace_file_1_0_bal,header = TRUE)

silacAvg <- rbind(silacAvg_0_2,silacAvg_0_3,silacAvg_0_4,silacAvg_0_6,silacAvg_0_8,silacAvg_1_0)
parameterSpace <- rbind(parameterSpace_0_2,parameterSpace_0_3,parameterSpace_0_4,parameterSpace_0_6,parameterSpace_0_8,parameterSpace_1_0)
parameterSpace_bal <- rbind(parameterSpace_0_2_bal,parameterSpace_0_3_bal,parameterSpace_0_4_bal,parameterSpace_0_6_bal,parameterSpace_0_8_bal,parameterSpace_1_0_bal)

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

## bistability (all cell cycles)
p <- ggplot(subset(parameterSpace_bal,FIRING_THRESHOLD !=0.3),aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_raster(aes(fill=bistability)) + 
  scale_y_log10("Methylation rate",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability (per histone per firing)",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=15/12) + facet_grid(. ~ FIRING_THRESHOLD) +
  scale_fill_gradientn(name="Bistability",colours=rev(rainbow(3)),limits=c(0,1)) +
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        strip.text.x = element_text(size=7, angle=0),
        strip.text.y = element_text(size=7, angle=0))
p
ggsave("Bistability.pdf",width=12,height=6,units="cm")

## me3_end (all cell cycles - initially repressed)
p <- ggplot(subset(parameterSpace,FIRING_THRESHOLD !=0.3),aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_raster(aes(fill=me3_end)) + 
  scale_y_log10("Methylation rate",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability (per histone per firing)",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=15/12) + facet_grid(. ~ FIRING_THRESHOLD) +
  scale_fill_gradientn(name="End K27me3",colours=rev(rainbow(3)),limits=c(0,1)) +
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        strip.text.x = element_text(size=7, angle=0),
        strip.text.y = element_text(size=7, angle=0))
p
ggsave("me3_end.pdf",width=12,height=6,units="cm")

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

p <- ggplot(subset(sum_sq_err,FIRING_THRESHOLD!=0.3),aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_raster(aes(fill=SSE)) + 
  scale_y_log10("Methylation rate",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability (per histone per firing)",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=15/12) + facet_grid(. ~ FIRING_THRESHOLD) +
  scale_fill_gradientn(name="SSE",colours=c("red","white","blue"),
                       trans="log",breaks=c(0.1,1,10)) +
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        strip.text.x = element_text(size=7, angle=0),
        strip.text.y = element_text(size=7, angle=0))
p
ggsave("SSE.pdf",width=12,height=6,units="cm")

fitSILAC <- sum_sq_err[sum_sq_err$SSE<0.9,]
fitSILAC_parameterSpace <- merge(fitSILAC,parameterSpace_bal,by=c("P_METHYLATE","P_DEMETHYLATE","FIRING_THRESHOLD"))
fitSILAC_parameterSpace[fitSILAC_parameterSpace$bistability==max(fitSILAC_parameterSpace$bistability),]

p_dem <- sort(unique(adjustedSILAC$P_DEMETHYLATE))[seq(6,16,by=1)]
p_me <- sort(unique(adjustedSILAC$P_METHYLATE))[seq(2,8,by=1)]

SILACtoPlotFits <- subset(adjustedSILAC,P_DEMETHYLATE %in% p_dem & P_METHYLATE %in% p_me)
SSEtoPlotFits  <- subset(sum_sq_err,P_DEMETHYLATE %in% p_dem & P_METHYLATE %in% p_me)
parameterSpace_baltoPlotFits <- subset(parameterSpace_bal,P_DEMETHYLATE %in% p_dem & P_METHYLATE %in% p_me)

## Add dummy data to experimental data frame for plotting over all parameter space
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
pars <- expand.grid(P_DEMETHYLATE=unique(SILACtoPlotFits$P_DEMETHYLATE),
                    P_METHYLATE=unique(SILACtoPlotFits$P_METHYLATE))
K27me3_Alabert_Dummy <- expand.grid.df(pars,K27me3_Alabert_summary)

ggplot(data=K27me3_Alabert_Dummy) + 
  geom_rect(data = subset(parameterSpace_baltoPlotFits,FIRING_THRESHOLD==1.0),
            aes(fill = bistability),
            xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.3) +
  scale_fill_gradientn(name="Bistability",colours=c("white","#95BA79"),breaks=c(0,0.5,1),limits=c(0,1)) +  
  geom_point(aes(x=time,y=adjusted_level,col=label,group=label),size=1.2) + 
  scale_y_continuous(name="K27me3 level (normalised)",breaks=c(0,1,2),limits=c(0,2.05)) +
  scale_x_continuous(name="Time (hours)",breaks=c(0,20,40)) +
  geom_errorbar(width=3,aes(x=time,y=adjusted_level,col=label,group=label,
                            ymin=adjusted_level-adjusted_sd,
                            ymax=adjusted_level+adjusted_sd)) +
  facet_grid(P_METHYLATE ~ P_DEMETHYLATE, labeller = scientific_labeller_factor, as.table = FALSE) + 
  geom_line(data=subset(SILACtoPlotFits,FIRING_THRESHOLD==1.0),
            aes(x=time,y=adjusted_level,col=label,group=label),size=0.4) + 
  scale_color_manual(values=c("blue2","orange"),name="Histones",
                     breaks=c("LIGHT", "HEAVY"),
                     labels=c("Old", "New")) + 
  theme_thesis_multiplanel

ggsave(file="ParameterSpaceFit_fir1_0.pdf",width=13.5,height=9,units="cm")


ggplot(data=K27me3_Alabert_Dummy) + 
  geom_rect(data = subset(parameterSpace_baltoPlotFits,FIRING_THRESHOLD==0.3),
            aes(fill = bistability),
            xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.3) +
  scale_fill_gradientn(name="Bistability",colours=c("white","#95BA79"),breaks=c(0,0.5,1),limits=c(0,1)) +  
  geom_point(aes(x=time,y=adjusted_level,col=label,group=label),size=1.2) + 
  scale_y_continuous(name="K27me3 level (normalised)",breaks=c(0,1,2),limits=c(0,2.05)) +
  scale_x_continuous(name="Time (hours)",breaks=c(0,20,40)) +
  geom_errorbar(width=3,aes(x=time,y=adjusted_level,col=label,group=label,
                            ymin=adjusted_level-adjusted_sd,
                            ymax=adjusted_level+adjusted_sd)) +
  facet_grid(P_METHYLATE ~ P_DEMETHYLATE, labeller = scientific_labeller_factor, as.table = FALSE) + 
  geom_line(data=subset(SILACtoPlotFits,FIRING_THRESHOLD==0.3),
            aes(x=time,y=adjusted_level,col=label,group=label),size=0.4) + 
  scale_color_manual(values=c("blue2","orange"),name="Histones",
                     breaks=c("LIGHT", "HEAVY"),
                     labels=c("Old", "New")) + 
  theme_thesis_multiplanel

ggsave(file="ParameterSpaceFit_fir0_3.pdf",width=13.5,height=9,units="cm")

# geom_rect(data = tp,aes(fill = day),xmin = -Inf,xmax = Inf,
#           ymin = -Inf,ymax = Inf,alpha = 0.3)
