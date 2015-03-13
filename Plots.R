# Remove all variables
rm(list=ls())
library(ggplot2)
library(reshape2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
four_color <- c(cbPalette[7],cbPalette[2],cbPalette[4],cbPalette[6])

# Set the working directory
setwd("~/local/Modelling/TwoStateCoupled/")

s <- 60
r <- format(100000,scientific=FALSE)

time_file <- paste("t_s",s,"r",r,".txt",sep="")
tDep_methyl_file <- paste("Meth_t_s",s,"r",r,".txt",sep="")
tDep_protein_file <- paste("RepBound_t_s",s,"r",r,".txt",sep="")
tDep_firing_file <- paste("Firing_t_s",s,"r",r,".txt",sep="")

time <- t(read.table(time_file))
PRC2 <- t(read.table(tDep_protein_file))
K27me3 <- t(read.table(tDep_methyl_file))
Firing <- t(read.table(tDep_firing_file))/86400

par(mfrow=c(3,1),mar=c(2,4,1,1)+0.1)
plot(time/86400,PRC2,type="l",ylim=c(0,1),col="blue3")
plot(time/86400,K27me3,type="l",ylim=c(0,1),col="red3")
plot(Firing,rep(1,length(Firing)),type="p",ylim=c(0,1),col="black")

# s

# plot(rowMeans(U),ylim=c(0,1),type="l",xlim=c(0,86400), col="orange3")
# lines(rowMeans(B), col="black")

# 
# hl <- read.table("histoneStateLifetime.txt")
# hl <- as.matrix(hl)
# rownames(hl) <- as.character(seq(1:ncol(hl))*0.01-0.01)
# colnames(hl) <- as.character(seq(1:ncol(hl))*0.01-0.01)
# hl_df <- melt(hl,varnames=c("modRecruit","enzymatic"),value.name="histoneLifetime")
# 
# theme_set(theme_bw(8))
# p <- ggplot(hl_df, aes(x=modRecruit,y=enzymatic))
# p <- p + geom_tile(aes(fill=histoneLifetime/3600)) + 
#   scale_fill_continuous(name="Histone state lifetime (hours)") +
#   scale_x_continuous( expand = c(0,0)) +
#   scale_y_continuous( expand = c(0,0))
# p
# 
# pl <- read.table("proteinStateLifetime.txt")
# pl <- as.matrix(pl)
# rownames(pl) <- as.character(seq(1:nrow(pl))*0.01-0.01)
# colnames(pl) <- as.character(seq(1:ncol(pl))*0.01-0.01)
# pl_df <- melt(pl,varnames=c("modRecruit","enzymatic"),value.name="proteinLifetime")
# 
# theme_set(theme_bw(8))
# p <- ggplot(pl_df, aes(x=modRecruit,y=enzymatic))
# p <- p + geom_tile(aes(fill=proteinLifetime/3600)) + 
#   scale_fill_continuous(name="Protein state lifetime (hours)") +
#   scale_x_continuous( expand = c(0,0)) +
#   scale_y_continuous( expand = c(0,0))
# p
# 
# gap <- read.table("gap.txt")
# g <- as.matrix(gap)
# rownames(g) <- as.character(seq(1:nrow(g))*0.01-0.01)
# colnames(g) <- as.character(seq(1:ncol(g))*0.01-0.01)
# g_df <- melt(g,varnames=c("modRecruit","enzymatic"),value.name="gap")
# 
# theme_set(theme_bw(8))
# p <- ggplot(g_df, aes(x=modRecruit,y=enzymatic))
# p <- p + geom_tile(aes(fill=gap)) + 
#   scale_fill_gradient(name="Gap",limits=c(0.98,1)) +
#   scale_x_continuous( expand = c(0,0)) +
#   scale_y_continuous( expand = c(0,0))
# p
# 
