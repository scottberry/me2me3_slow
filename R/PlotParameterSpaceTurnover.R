# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")

# Set the working directory
#setwd("~/local/Modelling/me2me3_slow/")
setwd("~/Network/group-share/berrys/Results/Bursty2_08032016/")

parameterSpace_file = "ParamOptimRes_s60ctrl60cc50a1_00b1_00thresh0_33turn0_00100000Rep_st30.txt"
burstySpace_file = "BurstyOptimRes_s60ctrl60cc50a1_00b1_00thresh0_33turn0_00100000Rep_st30.txt"
parameterSpace <- read.table(parameterSpace_file,header = TRUE)
burstySpace <- read.table(burstySpace_file,header = TRUE)
min_sim_time <- min(parameterSpace$avgInitM,parameterSpace$avgInitU)/3600

parameterSpace$H3_3gap <- abs((parameterSpace$avgH3_3_M - parameterSpace$avgH3_3_U)/(parameterSpace$avgH3_3_M + parameterSpace$avgH3_3_U))

parameterSpace <- cbind(parameterSpace,subset(burstySpace,select=-c(P_DEMETHYLATE,FIRING)))
parameterSpace$burstSize <- parameterSpace$totalFiring/parameterSpace$bursts
parameterSpace$firing_cellcycle <- parameterSpace$totalFiring/21
parameterSpace$FP <- parameterSpace$firstPassageM*parameterSpace$firstPassageU/(parameterSpace$avgInitM*parameterSpace$avgInitU)

scientific_10_plotmath <- function(x) {
  gsub("e", " %*% 10^", scientific_format()(x))
}

parameterSpace$DEMETH_lab <- scientific_10_plotmath(signif(parameterSpace$P_DEMETHYLATE,digits=2))
parameterSpace <- within(parameterSpace, DEMETH_lab <- reorder(DEMETH_lab, P_DEMETHYLATE))

selectDemethylationRate <- unique(parameterSpace$P_DEMETHYLATE)[2:6]
selectK_ON_MAX <- unique(parameterSpace$K_ON_MAX)
selectK_OFF <- unique(parameterSpace$K_OFF)
parameterSpace <- subset(parameterSpace,
                         P_DEMETHYLATE %in% selectDemethylationRate & K_ON_MAX %in% selectK_ON_MAX & K_OFF %in% selectK_OFF)

## burstSize
p1 <- ggplot(subset(parameterSpace),aes(x=K_ON_MAX,y=K_OFF))
p1 <- p1 + geom_tile(aes(fill=burstSize)) + 
  scale_y_log10("K_OFF",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("K_ON_MAX",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + 
  scale_fill_viridis(name="",option="magma", trans="log10") +
  facet_grid(. ~ DEMETH_lab, labeller = label_parsed) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        #strip.background = element_blank(),
        #strip.text.y = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave(plot=p1,file="BurstSize.pdf",width=16,height=3,units="cm")

## burstDuration
p1 <- ggplot(subset(parameterSpace),aes(x=K_ON_MAX,y=K_OFF))
p1 <- p1 + geom_tile(aes(fill=burstDuration)) + 
  scale_y_log10("K_OFF",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("K_ON_MAX",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + 
  scale_fill_viridis(name="",option="magma", trans="log10") +
  facet_grid(. ~ DEMETH_lab, labeller = label_parsed) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        #strip.background = element_blank(),
        #strip.text.y = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave(plot=p1,file="BurstDuration.pdf",width=16,height=3,units="cm")

## quiescentDuration
p1 <- ggplot(subset(parameterSpace),aes(x=K_ON_MAX,y=K_OFF))
p1 <- p1 + geom_tile(aes(fill=quiescentDuration)) + 
  scale_y_log10("K_OFF",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("K_ON_MAX",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + 
  scale_fill_viridis(name="",option="magma", trans="log10") +
  facet_grid(. ~ DEMETH_lab, labeller = label_parsed) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        #strip.background = element_blank(),
        #strip.text.y = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave(plot=p1,file="QuiescentDuration.pdf",width=16,height=3,units="cm")

## firing/cell-cycle
p1 <- ggplot(subset(parameterSpace),aes(x=K_ON_MAX,y=K_OFF))
p1 <- p1 + geom_tile(aes(fill=firing_cellcycle)) + 
  scale_y_log10("K_OFF",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("K_ON_MAX",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + 
  scale_fill_viridis(name="",option="magma", trans="log10") +
  facet_grid(. ~ DEMETH_lab, labeller = label_parsed) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        #strip.background = element_blank(),
        #strip.text.y = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave(plot=p1,file="Firing_cellcycle.pdf",width=16,height=3,units="cm")

## bistability
p1 <- ggplot(subset(parameterSpace),aes(x=K_ON_MAX,y=K_OFF))
p1 <- p1 + geom_tile(aes(fill=bistability)) + 
  scale_y_log10("K_OFF",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("K_ON_MAX",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + 
  scale_fill_viridis(name="",option="viridis",limits=c(0,1),breaks=c(0,0.5,1.0)) +
  facet_grid(. ~ DEMETH_lab, labeller = label_parsed) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        #strip.background = element_blank(),
        #strip.text.y = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave(plot=p1,file="Bistability.pdf",width=16,height=3,units="cm")


## bistability
p1 <- ggplot(subset(parameterSpace),aes(x=K_ON_MAX,y=K_OFF))
p1 <- p1 + geom_tile(aes(fill=FP)) + 
  scale_y_log10("K_OFF",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("K_ON_MAX",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + 
  scale_fill_viridis(name="",option="viridis",limits=c(0,1),breaks=c(0,0.5,1.0)) +
  facet_grid(. ~ DEMETH_lab, labeller = label_parsed) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        #strip.background = element_blank(),
        #strip.text.y = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave(plot=p1,file="FP.pdf",width=16,height=3,units="cm")



## probM
p1 <- ggplot(subset(parameterSpace),aes(x=K_ON_MAX,y=K_OFF))
p1 <- p1 + geom_tile(aes(fill=probM)) + 
  scale_y_log10("K_OFF",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("K_ON_MAX",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + 
  scale_fill_viridis(name="",option="viridis",limits=c(0,1),breaks=c(0,0.5,1.0)) +
  facet_grid(. ~ DEMETH_lab, labeller = label_parsed) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        #strip.background = element_blank(),
        #strip.text.y = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave(plot=p1,file="ProbM.pdf",width=16,height=3,units="cm")


## probU
p1 <- ggplot(subset(parameterSpace),aes(x=K_ON_MAX,y=K_OFF))
p1 <- p1 + geom_tile(aes(fill=probU)) + 
  scale_y_log10("K_OFF",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("K_ON_MAX",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) + 
  scale_fill_viridis(name="",option="viridis",limits=c(0,1),breaks=c(0,0.5,1.0)) +
  facet_grid(. ~ DEMETH_lab, labeller = label_parsed) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        #strip.background = element_blank(),
        #strip.text.y = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave(plot=p1,file="ProbU.pdf",width=16,height=3,units="cm")
