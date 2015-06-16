# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")

# Set the working directory
setwd("~/Network/group-share/berrys/me2me3_slow/")

# ctrl5_file = "ParamOptimRes_s60ctrl5cc50st16.txt"
# ctrl10_file = "ParamOptimRes_s60ctrl10cc50st16.txt"
# ctrl20_file = "ParamOptimRes_s60ctrl20cc50st16.txt"
# ctrl30_file = "ParamOptimRes_s60ctrl30cc50st16.txt"
ctrl60_file = "ParamOptimRes_s60ctrl60cc50a1_00st20.txt"

# parameterSpace_ctrl5 <- read.table(ctrl5_file,header = TRUE)
# parameterSpace_ctrl10 <- read.table(ctrl10_file,header = TRUE)
# parameterSpace_ctrl20 <- read.table(ctrl20_file,header = TRUE)
# parameterSpace_ctrl30 <- read.table(ctrl30_file,header = TRUE)
parameterSpace_ctrl60 <- read.table(ctrl60_file,header = TRUE)

# parameterSpace <- rbind(parameterSpace_ctrl5,
#                         parameterSpace_ctrl10,
#                         parameterSpace_ctrl20,
#                         parameterSpace_ctrl30,
#                         parameterSpace_ctrl60)

parameterSpace <- parameterSpace_ctrl60
parameterSpace$controlSites <- factor(parameterSpace$controlSites)

min_sim_time <- min(parameterSpace$avgInitM,parameterSpace$avgInitU)/3600

labeller_FIRING <- function(variable,value) {
  val <- round(value,digits=5)
  fold <- round(value/0.0004,digits=1)
  return(paste("f = ",val," (",fold,"-fold)",sep=""))
}

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

m_b <- round(rev(unique(parameterSpace$P_METHYLATE))[seq(1,length(unique(parameterSpace$P_METHYLATE)),by=5)],digits=6)
dem_b <- round(rev(unique(parameterSpace$P_DEMETHYLATE))[seq(1,length(unique(parameterSpace$P_DEMETHYLATE)),by=5)],digits=4)

## gap
p <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p <- p + geom_raster(aes(fill=gap)) + 
  scale_y_log10("Methylation rate (per second)",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability (per histone per firing)",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) +
  scale_fill_gradientn(name="Gap",colours=rev(rainbow(3)),limits=c(0,1)) + 
  facet_grid(FIRING ~ controlSites) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        strip.text.x = element_text(size=7, angle=0),
        strip.text.y = element_text(size=7, angle=0))
p
ggsave(plot=p,file="AvgGap.pdf",width=12,height=12,units="cm")

## Bistability
p1 <- ggplot(parameterSpace,aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p1 <- p1 + geom_raster(aes(fill=bistability)) + 
  scale_y_log10("Methylation rate (per second)",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability (per histone per firing)",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) +
  scale_fill_gradientn(name="Bistability",colours=rev(rainbow(3)),limits=c(0,1)) + 
  facet_grid(FIRING ~ controlSites) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        strip.text.x = element_text(size=7, angle=0),
        strip.text.y = element_text(size=7, angle=0))
p1
ggsave(plot=p1,file="Bistability.pdf",width=12,height=12,units="cm")


## Bistability
p1 <- ggplot(subset(parameterSpace,controlSites==60 & FIRING < 0.02),aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p1 <- p1 + geom_tile(aes(fill=bistability)) + 
  scale_y_log10("Methylation rate",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=0.2/0.15) +
  scale_fill_gradientn(name="Bistability",colours=rev(rainbow(3)),limits=c(0,1),breaks=seq(0,1,by=0.5)) + 
  facet_grid(. ~ FIRING) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.size = unit(0.3, "cm"))
p1
ggsave(plot=p1,file="Bistability_FIRING_noG2.pdf",width=14,height=4,units="cm")


## Lifetime
p1 <- ggplot(subset(parameterSpace,controlSites==60 & FIRING < 0.02),aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p1 <- p1 + geom_tile(aes(fill=(firstPassageM/(3600*16)/50)*(firstPassageU/(3600*16)/50))) + 
  scale_y_log10("Methylation rate",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=0.2/0.15) +
  scale_fill_gradientn(name="Combined t_FP",colours=rev(rainbow(3)),limits=c(0,1)) +
  facet_grid(. ~ FIRING) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.size = unit(0.3, "cm"))
p1
ggsave(plot=p1,file="FirstPassageM.pdf",width=14,height=4,units="cm")


saveRDS(object = subset(parameterSpace,controlSites==60 & FIRING < 0.02),file = "DNArep_noG2.rds")

