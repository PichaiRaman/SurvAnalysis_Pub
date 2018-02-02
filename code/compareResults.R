###################################################
# Author: Pichai Raman, Komal Rathi
# Function: images and plots from merged data
# Date: 01/31/2018
###################################################

#######################################################################
#1) Let's look at the correlation for discovery and validation for all cancer data sets
#2) Let's do ROC curves for each of the methods compared to gold standard lists
#######################################################################

# call Libraries
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gplots)
library(pheatmap)
library(GSEABase)
library(patchwork)
library(ggpubr)

source("RocOn.R")

############################################
# Part 1 - Reproducibility
############################################
# read in the DiscValid Data set
dvData <- read.delim("../data/MergedData/resultsDiscValid.txt")
allOdds <- c(1:dim(dvData)[2])
allOdds <- allOdds[allOdds%%2==T]
compareCor <- data.frame()

myCors <- c()
for(i in 1:length(allOdds)) {
  tmpNum <- allOdds[i]
  out <- cor.test(dvData[,tmpNum], dvData[,(tmpNum+1)], method="spearman")$estimate
  myCors <-c(myCors, out)
}
correlationResults <- data.frame(colnames(dvData)[allOdds], myCors)
colnames(correlationResults) <- c("Study", "Correlation")
correlationResults[,1] <- gsub("1", "", correlationResults[,1])
correlationResults[,"Algorithm"] <- rep(c("C-index","Cox Regression","D-index","Distribution Specific Cut","K-means","Kaplan Scan","Quantile Cut 25-75","Quantile Cut Median"), each = 4)
correlationResults[,"Cancer"] <- rep(c("Head & Neck", "Kidney", "Ovarian", "Prostate"), 8)

# bargraph of results
jpeg("../Figures/DV_Barchart.jpg", width = 7000, height = 2000,  res = 550)
DV_Barchart <- ggplot(correlationResults, aes(factor(Cancer), Correlation, fill=Cancer)) + 
  geom_bar(stat="identity", width=0.9) + 
  facet_grid(.~Algorithm) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -75, hjust = 0)) + 
  theme(legend.position = "none") + xlab("Cancer") +
  coord_flip() + 
  theme(axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 16, colour = "black"),
        strip.text = element_text(size = 14, colour = "black"))
DV_Barchart
dev.off()

# boxplot of results
jpeg("../Figures/DV_Boxplot.jpg", width = 4320, height = 1440,  res = 324)
ggplot(correlationResults, aes(factor(Algorithm), Correlation, fill=Algorithm)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle = 30, hjust = 1, vjust = 1)) + 
  theme(legend.position = "none") + xlab("Cancer") + 
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        strip.text = element_text(size = 14, colour = "black"))
dev.off()

# scatter plots
# first have to convert to a stacked format
dvDataFS <- dvData[, c(1,2)]
dvDataFS[,"comparison"] <- gsub("1", "", colnames(dvData[1]))
colnames(dvDataFS)[1:2] <- c("Disc", "Valid")


myVector <- c(3:64)[c(3:64)%%2==1]

for(i in 1:length(myVector)) {
  j <- myVector[i]
  tmpdvDataFS <- dvData[, c(j,j+1)]
  tmpdvDataFS[,"comparison"] <- gsub("1", "", colnames(dvData[j]))
  colnames(tmpdvDataFS)[1:2] <- c("Disc", "Valid")
  dvDataFS <- rbind(dvDataFS, tmpdvDataFS)
}

dvDataFS[,"Algorithm"] <- rep(c("C-index","Cox Regression","D-index","Distribution Specific Cut","K-means","Kaplan Scan","Quantile Cut 25-75","Quantile Cut Median"), each = 82124)
dvDataFS[,"Cancer"] <- rep(c(rep("Head & Neck", 20531), rep("Kidney", 20531), rep("Ovarian", 20531), rep("Prostate", 20531)), 8)

jpeg("../Figures/DV_Scatterplot.jpg", width = 7000, height = 5000,  res = 550)
DV_SCPlot <- ggplot(dvDataFS, aes((-1)*log10(Disc), (-1)*log10(Valid))) + 
  geom_point(size=.5) + 
  facet_grid("Cancer ~ Algorithm") + 
  theme_bw() + 
  scale_x_continuous(limits=c(0,10), breaks=c(0,5,10)) + 
  scale_y_continuous(limits=c(0,10)) + 
  xlab("Partition 1 (-log10 P-value)") + 
  ylab("Partition 2 (-log10 P-value)") + 
  geom_smooth(method="lm") + 
  theme(axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 16, colour = "black"),
        strip.text = element_text(size = 14, colour = "black"))
DV_SCPlot 
dev.off()

# patch the two plots
Fig2 <- DV_Barchart + DV_SCPlot + plot_layout(ncol = 1, heights = c(1, 3)) + plot_annotation(tag_levels = "A", tag_suffix = ")") 
ggsave("../Figures/Figure2.png", plot = Fig2, height = 15, width = 20, units = 'in', dpi = 750)

# last plot is plot of overlaps at various p-values
runHypGeom <- function(set, genes, n = 20000) {
  # number of white balls
  x <- length(intersect(genes, set))
  
  # white balls
  m <- length(genes)
  
  # black balls
  n2 <- n-m
  
  # balls drawn from the urn
  k <- length(set)

  out <- phyper(x-1, m, n2, k, lower.tail=F)
  setSize <- k
  overLap <- x
  numGenes <- m
  
  myRet <- c(setSize, numGenes, overLap, out)
  return(myRet)
}

myVector <- c(3:64)[c(3:64)%%2==1]

cutoffVect <- c(.0001, .001, .005, .01, .05)

tmpOut <- c()
tmpDF <- dvData[, c(1,2)]
for(k in 1:length(cutoffVect)) {
  myCut <- cutoffVect[k]
  tmpDFRep1 <- rownames(tmpDF[tmpDF[,1]<myCut,])
  tmpDFRep2 <- rownames(tmpDF[tmpDF[,2]<myCut,])
  output <- runHypGeom(tmpDFRep1,tmpDFRep2, 20531)
  output <- c(gsub("1", "", colnames(dvData[1])), myCut, output)
  tmpOut <- c(tmpOut, output)
}
overlapDF <- matrix(tmpOut, nrow=5, byrow=T)
overlapDF <- data.frame(overlapDF)
colnames(overlapDF) <- c("Comparison", "Cutoff", "Num Discovery", "Num Validation", "Overlap", "PValue")

for(i in 1:length(myVector)) {
  j <- myVector[i]
  tmpDF <- dvData[, c(j,j+1)]
  
  tmpOut <- c()
  for(k in 1:length(cutoffVect)) {
    myCut <- cutoffVect[k]
    tmpDFRep1 <- rownames(tmpDF[tmpDF[,1]<myCut,])
    tmpDFRep2 <- rownames(tmpDF[tmpDF[,2]<myCut,])
    output <- runHypGeom(tmpDFRep1,tmpDFRep2, 20531)
    output <- c(gsub("1", "", colnames(dvData[j])), myCut, output)
    tmpOut <- c(tmpOut, output)
  }
  tmpOut <- matrix(tmpOut, nrow=5, byrow=T)
  tmpOut <- data.frame(tmpOut)
  colnames(tmpOut) <- c("Comparison", "Cutoff", "Num Discovery", "Num Validation", "Overlap", "PValue")
  overlapDF <- rbind(overlapDF, tmpOut)
}


overlapDF[,"Algorithm"] <- rep(c("C-index","Cox Regression","D-index","Distribution Specific Cut","K-means","Kaplan Scan","Quantile Cut 25-75","Quantile Cut Median"), each = 20)
overlapDF[,"Cancer"] <- rep(c(rep("Head & Neck", 5), rep("Kidney", 5), rep("Ovarian", 5), rep("Prostate", 5)), 8)
overlapDF[,"PValue"] <- as.numeric(as.character(overlapDF[,"PValue"]))

overlapDF[,"Overlap"] <- as.numeric(as.character(overlapDF[,"Overlap"]))
overlapDF[,"Num Discovery"] <- as.numeric(as.character(overlapDF[,"Num Discovery"]))
overlapDF[,"Num Validation"] <- as.numeric(as.character(overlapDF[,"Num Validation"]))

overlapDF[,"InterPercent"] <- overlapDF[,"Overlap"]/(overlapDF[,"Num Discovery"]+overlapDF[,"Num Validation"]-overlapDF[,"Overlap"])

overlapDF[,"Cutoff"] <- factor(overlapDF[,"Cutoff"], levels=c(0.05,0.01,0.005,0.001,1e-04))

jpeg("../Figures/DV_BarchartCutoff.jpg", width = 4500, height = 2880,  res = 324)
ggplot(overlapDF, aes(Cutoff, InterPercent)) + 
  geom_bar(stat="identity") + 
  facet_grid(Cancer~Algorithm) + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -75, hjust = 0)) + 
  theme(legend.position = "none") + 
  xlab("Cancer") + ylab("Overlap Percentage")
dev.off()

################################################
# Part 2 - Cancer Driver Detection / Correlation
################################################

##################################################
# Code : ROC Curve in ggplot2
# No one has done this yet…seriously…what the f*#(
# 10/13/2015
##################################################

# read in results matrix
resA <- read.delim("../data/MergedData/resultsPrimary.txt")

# load in all gene signatures
load("../data/CancerGeneLists/GeneSigDB_GS.RData")

# main Lists
# Ovarian Cancer Gene lists Ovarian_Spentzos04_115genes
ovList <- geneIds(genesigdbSymbol[[649]])

# Kidney Cancer gene list Kidney_Zhao06_259genes
kiList <- geneIds(genesigdbSymbol[[1135]])

# Head & Neck gene list HeadandNeck_Chung06_42genes
hnList <- geneIds(genesigdbSymbol[[1449]])

# Prostate cancer gene list HeadandNeck_Chung06_42genes
prList <- geneIds(genesigdbSymbol[[256]])


createROCFrame <- function(data, myCol, cList) {
  data <- data[myCol]
  data[,1] <- (-1)*log10(data[,1])
  data[,"gene"] <- rownames(data)
  data[,"symbol"] <- sapply(data[,"gene"], FUN=getGeneName)
  data[,"labs"]<- as.numeric(data[,"symbol"]%in%cList)
  data <- data[,c(1,4)]
  colnames(data)[1] <- "preds"
  return(data)
}


getGeneName <- function(x) {
  x <- strsplit(x, split="\\|")[[1]][1]
  return(x)
}


hn <- rbind(data.frame(createROCFrame(resA, 1, hnList), method = "C-index"),
            data.frame(createROCFrame(resA, 5, hnList), method = "Cox Regression"),
            data.frame(createROCFrame(resA, 9, hnList), method = "D-index"),
            data.frame(createROCFrame(resA, 13, hnList), method = "Distribution Specific Cut"),
            data.frame(createROCFrame(resA, 17, hnList), method = "Kaplan-Scan"),
            data.frame(createROCFrame(resA, 21, hnList), method = "K-Means"),
            data.frame(createROCFrame(resA, 25, hnList), method = "Quantile 25th-75th"),
            data.frame(createROCFrame(resA, 29, hnList), method = "Median"))
hnOutput <- roconMult(hn, myTitle = "Head & Neck")
hnOutputForPub <- roconMult2(hn, myTitle = "Head & Neck")

ki <- rbind(data.frame(createROCFrame(resA, 2, kiList), method = "C-index"),
            data.frame(createROCFrame(resA, 6, kiList), method = "Cox Regression"),
            data.frame(createROCFrame(resA, 10, kiList), method = "D-index"),
            data.frame(createROCFrame(resA, 14, kiList), method = "Distribution Specific Cut"),
            data.frame(createROCFrame(resA, 18, kiList), method = "Kaplan-Scan"),
            data.frame(createROCFrame(resA, 22, kiList), method = "K-Means"),
            data.frame(createROCFrame(resA, 26, kiList), method = "Quantile 25th-75th"),
            data.frame(createROCFrame(resA, 30, kiList), method = "Median"))
kiOutput <- roconMult(ki, myTitle = "Kidney")
kiOutputForPub <- roconMult2(ki, myTitle = "Kidney")


ov <- rbind(data.frame(createROCFrame(resA, 3, ovList), method = "C-index"),
            data.frame(createROCFrame(resA, 7, ovList), method = "Cox Regression"),
            data.frame(createROCFrame(resA, 11, ovList), method = "D-index"),
            data.frame(createROCFrame(resA, 15, ovList), method = "Distribution Specific Cut"),
            data.frame(createROCFrame(resA, 19, ovList), method = "Kaplan-Scan"),
            data.frame(createROCFrame(resA, 23, ovList), method = "K-Means"),
            data.frame(createROCFrame(resA, 27, ovList), method = "Quantile 25th-75th"),
            data.frame(createROCFrame(resA, 31, ovList), method = "Median"))
ovOutput <- roconMult(ov, myTitle = "Ovarian")
ovOutputForPub <- roconMult2(ov, myTitle = "Ovarian")


pr <- rbind(data.frame(createROCFrame(resA, 4, prList), method = "C-index"),
            data.frame(createROCFrame(resA, 8, prList), method = "Cox Regression"),
            data.frame(createROCFrame(resA, 12, prList), method = "D-index"),
            data.frame(createROCFrame(resA, 16, prList), method = "Distribution Specific Cut"),
            data.frame(createROCFrame(resA, 20, prList), method = "Kaplan-Scan"),
            data.frame(createROCFrame(resA, 24, prList), method = "K-Means"),
            data.frame(createROCFrame(resA, 28, prList), method = "Quantile 25th-75th"),
            data.frame(createROCFrame(resA, 32, prList), method = "Median"))
prOutput <- roconMult(pr, myTitle = "Prostate")
prOutputForPub <- roconMult2(pr, myTitle = "Prostate")


# individual Figures
jpeg("../Figures/ROC_HeadNeck.jpg", width=2880, height=2880,  res=324)
hnOutput[[2]]
dev.off()
jpeg("../Figures/ROC_Kidney.jpg", width=2880, height=2880,  res=324)
kiOutput[[2]]
dev.off()
jpeg("../Figures/ROC_Ovarian.jpg", width=2880, height=2880,  res=324)
ovOutput[[2]]
dev.off()
jpeg("../Figures/ROC_Prostate.jpg", width=2880, height=2880,  res=324)
prOutput[[2]]
dev.off()

# all ROC Plots
jpeg("../Figures/ROC_All.jpg", width=4320, height=4320,  res=486)
grid.arrange(hnOutput[[2]] + theme(legend.position="none") + ggtitle("Head & Neck"), kiOutput[[2]] + theme(legend.position="none") + ggtitle("Kidney"), ovOutput[[2]] + theme(legend.position="none") + ggtitle("Ovarian"), prOutput[[2]] + theme(legend.position="none") + ggtitle("Prostate"), nrow=2)
dev.off()

# with legend
jpeg("../Figures/Figure5.png", width = 4700, height = 4220,  res = 486)
ggarrange(hnOutputForPub[[2]], kiOutputForPub[[2]], ovOutputForPub[[2]], prOutputForPub[[2]], 
          nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
dev.off()

# Table for All Plots
allAUCTable <- cbind(hnOutput[[3]][,c(2,1)], kiOutput[[3]][1], ovOutput[[3]][1], prOutput[[3]][1])
allAUCTable[,"Average"] <- rowMeans(allAUCTable[2:5])
colnames(allAUCTable)<- c("Method", "Head & Neck", "Kidney", "Ovarian", "Prostate", "Mean")
allAUCTable <- allAUCTable[order(-allAUCTable[,"Mean"]),]
write.table(allAUCTable, "../Tables/AUC_Table.txt", sep="\t", row.names=F)

# Anova
library(reshape2)
allAUCTableTS <- melt(allAUCTable)
pval <- summary(aov(lm(value ~ Method, data=allAUCTableTS)))[[1]][[5]][1]
write.table(data.frame(TukeyHSD(aov(lm(value ~ Method, data=allAUCTableTS)))[[1]]), "../Tables/AUC_Tukey.txt", sep="\t", row.names=T)

#########################################
# Code : Finally correlation between methods
# Heatmaps
#########################################

# there are some NAs and a zero value in c-index (total 345 rows)
# pretty heatmap
row.annotation <- data.frame(ann = colnames(resA))
row.annotation <- cbind(row.annotation, colsplit(row.annotation$ann, '_', c('method','cancer')))
row.annotation$cancer[row.annotation$cancer == "ki"] <- "ki = Kidney"
row.annotation$cancer[row.annotation$cancer == "hn"] <- "hn = Head & Neck"
row.annotation$cancer[row.annotation$cancer == "ov"] <- "ov = Ovarian"
row.annotation$cancer[row.annotation$cancer == "pr"] <- "pr = Prostate"
row.annotation$method[row.annotation$method == "cindex"] <- "cindex = C-index"
row.annotation$method[row.annotation$method == "coxreg"] <- "coxreg = Cox Regression"
row.annotation$method[row.annotation$method == "dindex"] <- "dindex = D-index"
row.annotation$method[row.annotation$method == "dmod"] <- "dmod = Distribution Specific Cutoff"
row.annotation$method[row.annotation$method == "kmeans"] <- "kmeans = K-Means"
row.annotation$method[row.annotation$method == "kmScan"] <- "kmScan = Kaplan Scan"
row.annotation$method[row.annotation$method == "qCut2575"] <- "qCut2575 = Quantile-cut 25-75"
row.annotation$method[row.annotation$method == "qCut50"] <- "qCut50 = Quantile-cut 50"
rownames(row.annotation) <- row.annotation$ann
row.annotation$ann <- NULL

resA[is.na(resA)] <- 1
resA[resA == 0] <- 1
resALog <- (-1)*log10(resA)
jpeg("../Figures/Figure3.jpg", width=5500, height=4120,  res=486)
pheatmap(cor(resALog), annotation_row = row.annotation)
dev.off()






















