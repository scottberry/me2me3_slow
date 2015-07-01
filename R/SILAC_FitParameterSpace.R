# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")
library(grid)
library(gtable)
library(DEoptim)

# Set the working directory
setwd("~/local/Modelling/me2me3_slow/")

s <- 60
ctrl <- 60
cc <- 15
a <- 1.0
st <- 18

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")

silacAvg_file <- paste("SilacRelAverage_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
parameterSpace_file <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,".txt",sep="")
#parameterSpace_balanced_file <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,"st",st,"_balanced.txt",sep="")

silacAvg <- read.table(silacAvg_file,header = TRUE)
parameterSpace <- read.table(parameterSpace_file,header = TRUE)
#parameterSpace_balanced <- read.table(parameterSpace_balanced_file,header = TRUE)

## Extract the cell-cycle-end K27me3 level for this parameter set
endVal <- subset(parameterSpace,select=c(FIRING,P_DEMETHYLATE,P_METHYLATE,FIRING_THRESHOLD,me3_end))

## Adjust relative SILAC values by this asymptotic value
adjustedSILAC <- merge(silacAvg,endVal,by=c("FIRING","P_DEMETHYLATE","P_METHYLATE","FIRING_THRESHOLD"))

adjustedSILAC$adjusted_level <- adjustedSILAC$level/adjustedSILAC$me3_end
adjustedSILAC

## Read experimental data file
K27me3_Alabert <- readRDS(file = "R/K27me3_SILAC_expt.rds")
K27me3_Alabert$Time <- as.numeric(as.character(K27me3_Alabert$Time))
colnames(K27me3_Alabert) <- c("mod","time","label","level")
K27me3_Alabert$mod <- revalue(K27me3_Alabert$mod,c("K27me3"="me3"))
K27me3_Alabert$label <- revalue(K27me3_Alabert$label,c("H"="HEAVY","L"="LIGHT"))

## We do not know the absolute level to adjust relative to the cell-cycle end level
## which is equal to the K27me3 (light) level at time = 0

## Extract the cell-cyle-end K27me3 from the experimental data and use this to normalise
me3_light_0 <- mean(K27me3_Alabert[K27me3_Alabert$label=="LIGHT" & K27me3_Alabert$time==0,]$level)
K27me3_Alabert$adjusted_level <- K27me3_Alabert$level/me3_light_0

## Add dummy data to experimental data frame for plotting ovcer parameter space
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
pars <- expand.grid(P_DEMETHYLATE=unique(adjustedSILAC$P_DEMETHYLATE),
                    FIRING_THRESHOLD=unique(adjustedSILAC$FIRING_THRESHOLD))
K27me3_Alabert_Dummy <- expand.grid.df(pars,K27me3_Alabert)

ggplot(K27me3_Alabert_Dummy,aes(x=time,y=adjusted_level,col=label,group=label)) + 
  geom_point() + 
  facet_grid(FIRING_THRESHOLD ~ P_DEMETHYLATE) + 
  geom_line(data=adjustedSILAC) + 
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
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=FIRING_THRESHOLD))
p <- p + geom_raster(aes(fill=bistability)) + 
#  scale_y_log10("Methylation rate (per second)",
#                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability (per histone per firing)",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + #facet_grid(FIRING_THRESHOLD ~ .) +
  scale_fill_gradientn(name="Bistability",colours=rev(rainbow(3)),limits=c(0,1)) +
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        strip.text.x = element_text(size=7, angle=0),
        strip.text.y = element_text(size=7, angle=0))
p
ggsave("Threshold.pdf",width=7,height=18,units="cm")
p

## Quantify model fit to data
summary(K27me3_Alabert)
summary(adjustedSILAC)

comparison <- merge(subset(K27me3_Alabert,select=c(mod,time,label,adjusted_level)),
      subset(adjustedSILAC,select=c(mod,time,label,adjusted_level,P_DEMETHYLATE,FIRING_THRESHOLD)),
      by=c("mod","time","label"))

names(comparison)
comparison$sq_err <- (comparison$adjusted_level.x - comparison$adjusted_level.y)^2
sum_sq_err <- aggregate(comparison$sq_err,by=list(comparison$P_DEMETHYLATE,
                                    comparison$FIRING_THRESHOLD),
          FUN = sum)
colnames(sum_sq_err) <- c("P_DEMETHYLATE","FIRING_THRESHOLD","SSE")

sum_sq_err[sum_sq_err$SSE==min(sum_sq_err$SSE),]

# geom_rect(data = tp,aes(fill = day),xmin = -Inf,xmax = Inf,
#           ymin = -Inf,ymax = Inf,alpha = 0.3)
