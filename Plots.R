# Remove all variables
rm(list=ls())
library(ggplot2)
library(reshape2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
four_color <- c(cbPalette[7],cbPalette[2],cbPalette[4],cbPalette[6])

# Set the working directory
setwd("~/Network/group-share/berrys/me2me3_slow/")

s <- 60
ctrl <- 60
cc <- 100 
st <- 1

time_file <- paste("t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")
tDep_me0_file <- paste("me0_t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")
tDep_me1_file <- paste("me1_t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")
tDep_me2_file <- paste("me2_t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")
tDep_me3_file <- paste("me3_t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")
tDep_firing_file <- paste("Firing_t_s",s,"ctrl",ctrl,"cc",cc,"st",st,".txt",sep="")

time <- read.table(time_file) ; colnames(time) <- "time"
me0 <- read.table(tDep_me0_file) ; colnames(me0) <- "level"
me1 <- read.table(tDep_me1_file) ; colnames(me1) <- "level"
me2 <- read.table(tDep_me2_file) ; colnames(me2) <- "level"
me3 <- read.table(tDep_me3_file) ; colnames(me3) <- "level"
Firing <- read.table(tDep_firing_file)

par(mfrow=c(5,1),mar=c(2,4,0,0)+0.5,oma=c(3,3,0,0))
plot(time$time/3600,me0$level,type="l",ylim=c(0,1),col="blue3")
plot(time$time/3600,me1$level,type="l",ylim=c(0,1),col="blue3")
plot(time$time/3600,me2$level,type="l",ylim=c(0,1),col="red3")
plot(time$time/3600,me3$level,type="l",ylim=c(0,1),col="red3")
hist(Firing$V1/3600,seq(0,(max(Firing)/3600 + 1),by=1),ylab="Firing/hour",main="")
mtext("time (hours)",side=1,line=0,outer=TRUE)
mtext("modification level",side=2,line=0,outer=TRUE)



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
