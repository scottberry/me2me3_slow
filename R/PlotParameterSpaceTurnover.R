# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")

# Set the working directory
setwd("~/local/Modelling/me2me3_slow/")

parameterSpace_file = "ParamOptimRes_s60ctrl60cc50a1_00b1_00thresh1_00turn0_00100000Rep_st30.txt"
parameterSpace <- read.table(parameterSpace_file,header = TRUE)
min_sim_time <- min(parameterSpace$avgInitM,parameterSpace$avgInitU)/3600

labeller_FIRING <- function(variable,value) {
  val <- round(value,digits=5)
  fold <- round(value/0.0004,digits=1)
  return(paste("f = ",val," (",fold,"-fold)",sep=""))
}

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

parameterSpace$H3_3gap <- abs((parameterSpace$avgH3_3_M - parameterSpace$avgH3_3_U)/(parameterSpace$avgH3_3_M + parameterSpace$avgH3_3_U))

## Bistability
p1 <- ggplot(subset(parameterSpace,controlSites==60 & FIRING < 0.02),aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p1 <- p1 + geom_tile(aes(fill=histoneGap)) + 
  scale_y_log10("Methylation rate",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=0.15/0.12) +
  scale_fill_gradientn(name="",colours=rev(rainbow(3)),limits=c(0,1),breaks=seq(0,1,by=0.5)) + 
  facet_grid(. ~ FIRING) + 
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
p1
ggsave(plot=p1,file="Bistability.pdf",width=16,height=3,units="cm")


## Histone Lifetime
p1 <- ggplot(subset(parameterSpace,controlSites==60 & FIRING < 0.02),aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p1 <- p1 + geom_tile(aes(fill=1/(3600*totTurnover))) + 
  scale_y_log10("Methylation rate",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=0.15/0.12) +
  scale_fill_gradientn(name="",colours=c("red","white","blue"),trans="log",breaks=c(1,10,100,1000)) + 
  facet_grid(. ~ FIRING) + 
  theme_bw(7) + theme_thesis_multiplanel +
  theme(plot.title = element_text(lineheight=.8, face="bold"),
        panel.margin=unit(0,"lines"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle=0),
        #strip.text.x = element_text(size=7, angle=0),
        #strip.text.y = element_text(size=7, angle=0),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        legend.key.size = unit(0.3, "cm"))
p1
ggsave(plot=p1,file="HistoneLifetime.pdf",width=16,height=3,units="cm")

