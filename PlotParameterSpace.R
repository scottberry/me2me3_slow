# Remove all variables
rm(list=ls())
library(ggplot2)
library(reshape2)
library(grid)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
four_color <- c(cbPalette[7],cbPalette[2],cbPalette[4],cbPalette[6])

# Set the working directory
setwd("~/Network/group-share/berrys/me2me3_slow/")

file= "ParamOptimRes_s60r50000tr0_000st10.txt"
parameterSpace <- read.table(file,header = TRUE)

parameterSpace$firstPassageM[parameterSpace$firstPassageM == -1] <- max(parameterSpace$firstPassageM)
parameterSpace$firstPassageU[parameterSpace$firstPassageU == -1] <- max(parameterSpace$firstPassageU)

## gap
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_tile(aes(fill=gap)) + scale_y_log10() + scale_x_log10() +
  scale_fill_gradientn(colours=rev(rainbow(3))) + 
  facet_grid(FIRING ~ . ,labeller=label_both) + 
  ggtitle("Average gap") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines")) + theme_bw(8)
p
ggsave(plot=p,file="AvgGap.pdf",width=8,height=8)


## avgM
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_tile(aes(fill=Mavg)) + scale_y_log10() + scale_x_log10() +
  scale_fill_gradientn(colours=rev(rainbow(3)),limits=c(0,1)) + 
  facet_grid(FIRING ~ .,labeller=label_both) + 
  ggtitle("Average M level over all loci") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines")) + theme_bw(8)
p
ggsave(plot=p,file="AvgM.pdf",width=8,height=8)

## sim time U state
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_tile(aes(fill=avgInitU/3600)) + scale_y_log10() + scale_x_log10() +
  scale_fill_gradientn(colours=rev(rainbow(3)),trans="log") + 
  facet_grid(FIRING ~ ., labeller=label_both) + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines")) + theme_bw(8)
p

## sim time M state
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_tile(aes(fill=avgInitM/3600)) + scale_y_log10() + scale_x_log10() +
  scale_fill_gradientn(colours=rev(rainbow(3)),trans="log") + 
  facet_grid(FIRING ~ ., labeller=label_both) + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines")) + theme_bw(8)
p

min_sim_time <- min(parameterSpace$avgInitM,parameterSpace$avgInitU)/3600

## firstpassage M
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_tile(aes(fill=firstPassageM/3600)) + scale_y_log10() + scale_x_log10() +
  scale_fill_gradientn(colours=rev(rainbow(3)),trans="log",limits=c(1,min_sim_time)) + 
  facet_grid(FIRING ~ ., labeller=label_both) + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines")) + theme_bw(8)
p

## firstpassage U
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_tile(aes(fill=firstPassageU/3600)) + scale_y_log10() + scale_x_log10() +
  scale_fill_gradientn(colours=rev(rainbow(3)),trans="log",limits=c(1,min_sim_time)) + 
  facet_grid(FIRING ~ ., labeller=label_both) + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines")) + theme_bw(8)
p

## Gap on subset
sub <- subset(parameterSpace,firstPassageU > 3600*min_sim_time & firstPassageM > 3600*min_sim_time)

## Bistability
m_b <- round(rev(unique(parameterSpace$P_METHYLATE))[c(1,3,5,7)],digits=5)
dem_b <- round(rev(unique(parameterSpace$P_DEMETHYLATE))[c(1,3,5,7)],digits=4)
p1 <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p1 <- p1 + geom_tile(aes(fill=bistability)) + scale_y_log10(breaks=m_b) + scale_x_log10(breaks=dem_b) +
  scale_fill_gradientn(name="Hours",colours=rev(rainbow(3)),limits=c(0,1)) + 
  facet_grid(FIRING ~ .,labeller=label_both) + theme_bw(7) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=90),
        strip.text.x = element_text(size=7, angle=90),
        strip.text.y = element_text(size=7, angle=0))
p1
ggsave(plot=p1,file="Bistability.pdf",width=8,height=8)

