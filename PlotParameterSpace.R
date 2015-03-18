# Remove all variables
rm(list=ls())
library(ggplot2)
library(reshape2)
library(grid)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
four_color <- c(cbPalette[7],cbPalette[2],cbPalette[4],cbPalette[6])

# Set the working directory
setwd("~/local/Modelling/TwoStateCoupled/")

file= "ParamOptimRes_s60r100000tr0_600st6.txt"
parameterSpace <- read.table(file,header = TRUE)

parameterSpace$firstPassageM[parameterSpace$firstPassageM == -1] <- max(parameterSpace$firstPassageM)
parameterSpace$firstPassageU[parameterSpace$firstPassageU == -1] <- max(parameterSpace$firstPassageU)

## gap
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=ENZYMATIC))
p <- p + geom_tile(aes(fill=gap)) + scale_y_log10() + 
  scale_fill_gradientn(colours=rev(rainbow(3))) + 
  facet_grid(R_OFF ~ FIRING,labeller=label_both) + 
  ggtitle("Average gap (P_OFF = 0.6)") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines")) + theme_bw(8)
p

## avgM
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=ENZYMATIC))
p <- p + geom_tile(aes(fill=Mavg)) + scale_y_log10() + 
  scale_fill_gradientn(colours=rev(rainbow(3)),limits=c(0,1)) + 
  facet_grid(R_OFF ~ FIRING,labeller=label_both) + 
  ggtitle("Average M level over all loci (P_OFF = 0.6)") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines")) + theme_bw(8)
p

## first passage M
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=ENZYMATIC))
p <- p + geom_tile(aes(fill=firstPassageM)) + scale_y_log10() + 
  scale_fill_gradientn(colours=rev(rainbow(3)),trans="log") + 
  facet_grid(R_OFF ~ FIRING,labeller=label_both) + 
  ggtitle("First passage time (initial M) (P_OFF = 0.6)") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines")) + theme_bw(8)
p

## first passage U
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=ENZYMATIC))
p <- p + geom_tile(aes(fill=firstPassageU)) + scale_y_log10() + 
  scale_fill_gradientn(colours=rev(rainbow(3)),trans="log") + 
  facet_grid(R_OFF ~ FIRING,labeller=label_both) + 
  ggtitle("First passage time (initial U) (P_OFF = 0.6)") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines")) + theme_bw(8)
p
