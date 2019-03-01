############################################
# Author: Pichai Raman
# Let's do kmeans and get difference of means
# maybe that's why kidney does so much better
# Date: 01/31/2018
#############################################

library(ggplot2)

load("../data/ALLDATA.RData")

set.seed(30)
getClustMeanDiff <- function(x) {
  x <- as.numeric(x)
  x <- x+abs(rnorm(length(x), mean=.5, sd=.000001)) #added for when its all 0
  out <- kmeans(x, 2)
  abs(log2(out$centers[1]/out$centers[2]))
}

hnMD <- apply(exprs_hn, FUN=getClustMeanDiff, MARGIN=1)
ovMD <- apply(exprs_ov, FUN=getClustMeanDiff, MARGIN=1)
prMD <- apply(exprs_pr, FUN=getClustMeanDiff, MARGIN=1)
kiMD <- apply(exprs_ki, FUN=getClustMeanDiff, MARGIN=1)

data <- data.frame(c(hnMD, ovMD, prMD, kiMD), c(rep("HN", 20531), rep("OV", 20531), rep("PR", 20531), rep("KI", 20531)))
colnames(data) <- c("Mean Difference", "cancer")

tmp <- data
tmp <- tmp[tmp[,1]>1,]
output <- data.frame(t(table(tmp[,2])/sum(table(tmp[,2]))), i=1)

for(i in 2:12) {
  tmp <- data
  tmp <- tmp[tmp[,1]>i,]
  tmpOut <- data.frame(t(table(tmp[,2])/sum(table(tmp[,2]))), i)
  output <- rbind(output, tmpOut)
}

jpeg("../Figures/MeanDiffPerClusteredGenePerCancer.jpg", width=4320, height=4320,  res=486)
p <- ggplot(output, aes(i, Freq, color=Var2)) + 
  geom_point() + 
  geom_smooth(se=F) + theme_bw() + 
  xlab("Log Fold Change Threshold") + 
  ylab("Percentage of Number of samples greater than threshold") + 
  ggtitle("Number of highly variable genes")
p <- p + labs(color="Cancer Type")
p
dev.off()

