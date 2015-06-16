# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")

# Set the working directory
setwd("~/local/Modelling/me2me3_slow/")

file= "ParamOptimRes_s60ctrl60cc50a1_00st16_Proc.txt"
parameterSpace <- read.table(file,header = TRUE)

parameterSpace$firstPassageM[parameterSpace$firstPassageM == -1] <- max(parameterSpace$firstPassageM)
parameterSpace$firstPassageU[parameterSpace$firstPassageU == -1] <- max(parameterSpace$firstPassageU)

min_sim_time <- min(parameterSpace$avgInitM,parameterSpace$avgInitU)/3600

labeller_FIRING <- function(variable,value) {
  val <- round(value,digits=5)
  fold <- round(value/0.0004,digits=1)
  return(paste("f = ",val," (",fold,"-fold)",sep=""))
}
m_b <- round(rev(unique(parameterSpace$P_METHYLATE))[seq(1,length(unique(parameterSpace$P_METHYLATE)),by=3)],digits=6)
dem_b <- round(rev(unique(parameterSpace$P_DEMETHYLATE))[seq(1,length(unique(parameterSpace$P_DEMETHYLATE)),by=3)],digits=4)

## gap
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_raster(aes(fill=gap)) + 
  scale_y_log10(breaks=m_b) + 
  scale_x_log10(breaks=dem_b) +
  coord_fixed(ratio=1) +
  scale_fill_gradientn(name="Gap",colours=rev(rainbow(3)),limits=c(0,1)) + 
  facet_grid(FIRING ~ .,labeller=labeller_FIRING) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=90),
        strip.text.x = element_text(size=7, angle=90),
        strip.text.y = element_text(size=7, angle=0))
p
ggsave(plot=p,file="AvgGap.pdf",width=5,height=8)

## Mavg
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_raster(aes(fill=Mavg)) + 
  scale_y_log10(breaks=m_b) + 
  scale_x_log10(breaks=dem_b) +
  coord_fixed(ratio=1) +
  scale_fill_gradientn(name="Average(me2/me3)",colours=rev(rainbow(3)),limits=c(0,1)) + 
  facet_grid(FIRING ~ .,labeller=labeller_FIRING) + 
  theme_bw(7) + theme_thesis_multiplanel + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=90),
        strip.text.x = element_text(size=7, angle=90),
        strip.text.y = element_text(size=7, angle=0))
p
ggsave(plot=p,file="AvgM.pdf",width=5,height=8)

## Bistability
p1 <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p1 <- p1 + geom_raster(aes(fill=bistability)) + 
  scale_y_log10(breaks=m_b) + 
  scale_x_log10(breaks=dem_b) +
  coord_fixed(ratio=1) +
  scale_fill_gradientn(name="Bistability",colours=rev(rainbow(3)),limits=c(0,1)) + 
  facet_grid(FIRING ~ .,labeller=labeller_FIRING) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=90),
        strip.text.x = element_text(size=7, angle=90),
        strip.text.y = element_text(size=7, angle=0))
p1
ggsave(plot=p1,file="Bistability.pdf",width=8,height=8)


