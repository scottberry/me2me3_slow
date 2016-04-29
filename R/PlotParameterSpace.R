# Remove all variables
rm(list=ls())

source("~/local/Thesis/R/ThesisTheme.R")

# Set the working directory
#setwd("~/local/Modelling/me2me3_slow/")
setwd("~/local/Manuscripts/TranscriptionOpposing2015/Figures/Processivity/st30/ProcMeth_st30_28042016/")

s <- 60
ctrl <- 60
cc <- 50
a <- 1.0
b <- 1.0
f <- 0.333
rep <- "Rep"
st <- 30
turn <- 0.001

astr <- paste('a',gsub("\\.", "_",sprintf("%0.2f",a)),sep="")
bstr <- paste('b',gsub("\\.", "_",sprintf("%0.2f",b)),sep="")
fstr <- paste('thresh',gsub("\\.", "_",sprintf("%0.2f",f)),sep="")
turnstr <- paste('turn',gsub("\\.", "_",sprintf("%0.8f",turn)),sep="")

parameterSpace_file <- paste("ParamOptimRes_s",s,"ctrl",ctrl,"cc",cc,astr,bstr,fstr,turnstr,rep,"_st",st,".txt",sep="")
parameterSpace <- read.table(parameterSpace_file,header=TRUE)
min_sim_time <- min(parameterSpace$avgInitM,parameterSpace$avgInitU)/3600

labeller_FIRING <- function(variable,value) {
  val <- round(value,digits=5)
  fold <- round(value/0.0004,digits=1)
  return(paste("f = ",val," (",fold,"-fold)",sep=""))
}

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

## Bistability
p1 <- ggplot(subset(parameterSpace,FIRING < 0.02),aes(x=P_DEMETHYLATE,y=P_METHYLATE))
p1 <- p1 + geom_tile(aes(fill=bistability)) + 
  scale_y_log10("Methylation rate",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10("Demethylation probability",
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed(ratio=1) +
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
ggsave(plot=p1,file="Bistability.pdf",width=14,height=4,units="cm")
# 
# ## Lifetime
# p1 <- ggplot(subset(parameterSpace,controlSites==60 & FIRING < 0.02),aes(x=P_DEMETHYLATE,y=P_METHYLATE))
# p1 <- p1 + geom_tile(aes(fill=(firstPassageM/(3600*22)/51)*(firstPassageU/(3600*22)/51))) + 
#   scale_y_log10("Methylation rate",
#                 labels = trans_format("log10", math_format(10^.x))) +
#   scale_x_log10("Demethylation probability",
#                 labels = trans_format("log10", math_format(10^.x))) +
#   coord_fixed(ratio=0.2/0.15) +
#   scale_fill_gradientn(name="First Passage",colours=rev(rainbow(3)),limits=c(0,1),breaks=seq(0,1,by=0.5)) +
#   facet_grid(. ~ FIRING) + 
#   theme_bw(7) + theme_thesis_multiplanel +
#   theme(plot.title = element_text(lineheight=.8, face="bold"),
#         panel.margin=unit(0,"lines"),
#         plot.margin = unit(c(0,0,0,0), "cm"),
#         axis.text.x = element_text(angle=0),
#         #strip.text.x = element_text(size=7, angle=0),
#         #strip.text.y = element_text(size=7, angle=0),
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         legend.key.size = unit(0.3, "cm"))
# p1
# ggsave(plot=p1,file="FirstPassage.pdf",width=14,height=4,units="cm")
# 
