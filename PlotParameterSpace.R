# Remove all variables
rm(list=ls())
library(ggplot2)
library(reshape2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
four_color <- c(cbPalette[7],cbPalette[2],cbPalette[4],cbPalette[6])

# Set the working directory
setwd("~/local/Modelling/TwoStateCoupled/")

file <- "ParameterTable.txt"
parameterSpace <- read.table(file,header = TRUE)

p <- ggplot(subset(parameterSpace,Mavg<0.6 & Mavg>0.4),aes(x=TRANSCRIPTION,y=ENZYMATIC))
p <- p + geom_tile(aes(fill=gap)) + 
  scale_fill_gradient() + 
  facet_wrap(R_OFF ~ FIRING)
p
