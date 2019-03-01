#########################################
# Author: Pichai Raman
# Function: Code : in silico data, ROC
# Date: 10/13/2015
#########################################

# Packages
source("rocon.R")
library(gridExtra)
library(ggplot2)

# Read in data
ResnoNoise <- read.delim("../data/SyntheticDataMerge/Res_0Noise.txt", stringsAsFactors = F)
Res10Noise <- read.delim("../data/SyntheticDataMerge/Res_0.1Noise.txt", stringsAsFactors = F)
Res25Noise <- read.delim("../data/SyntheticDataMerge/Res_0.25Noise.txt", stringsAsFactors = F)
Res50Noise <- read.delim("../data/SyntheticDataMerge/Res_0.5Noise.txt", stringsAsFactors = F)
Res75Noise <- read.delim("../data/SyntheticDataMerge/Res_0.75Noise.txt", stringsAsFactors = F)
Res100Noise <- read.delim("../data/SyntheticDataMerge/Res_1Noise.txt", stringsAsFactors = F)
Res150Noise <- read.delim("../data/SyntheticDataMerge/Res_1.5Noise.txt", stringsAsFactors = F)

# Read in positive controls
posControls <- read.delim("../data/SyntheticData/PositiveControls.txt")
posControlInd <- posControls[,1]
posControlList <- rownames(ResnoNoise)[posControlInd]

# Let's plot ROC curves for all per noise level
createROCFrame <- function(data, myCol, cList) {
  data <- data[myCol]
  data[,1] <- (-1)*log10(data[,1])
  data[,"gene"] <- rownames(data)
  data[,"symbol"] <- sapply(data[,"gene"], FUN=getGeneName)
  data[,"labs"]<- as.numeric(data[,"gene"]%in%cList)
  data <- data[,c(1,4)]
  colnames(data)[1] <- "preds"
  return(data)
}


getGeneName <- function(x) {
  x <- strsplit(x, split="\\.")[[1]][1]
  return(x)
}

##################################################
# Create ROC Objects for different noise levels
##################################################

noNoiseDF <- rbind(data.frame(createROCFrame(ResnoNoise, 2, posControlList), method="K-Means"),
                   data.frame(createROCFrame(ResnoNoise, 3, posControlList), method="Cox Regression"),
                   data.frame(createROCFrame(ResnoNoise, 4, posControlList), method="Median"),
                   data.frame(createROCFrame(ResnoNoise, 5, posControlList), method="Quantile 25th-75th"),
                   data.frame(createROCFrame(ResnoNoise, 7, posControlList), method="Kaplan-Scan"),
                   data.frame(createROCFrame(ResnoNoise, 8, posControlList), method="C-index"),
                   data.frame(createROCFrame(ResnoNoise, 9, posControlList), method="D-index"),
                   data.frame(createROCFrame(ResnoNoise, 10, posControlList), method="Distribution Specific Cut"))
noNoiseROC <- roconMult2(noNoiseDF, myTitle="ROC No Noise")

Noise10PercDF <- rbind(data.frame(createROCFrame(Res10Noise, 2, posControlList), method="K-Means"),
                       data.frame(createROCFrame(Res10Noise, 3, posControlList), method="Cox Regression"),
                       data.frame(createROCFrame(Res10Noise, 4, posControlList), method="Median"),
                       data.frame(createROCFrame(Res10Noise, 5, posControlList), method="Quantile 25th-75th"),
                       data.frame(createROCFrame(Res10Noise, 7, posControlList), method="Kaplan-Scan"),
                       data.frame(createROCFrame(Res10Noise, 8, posControlList), method="C-index"),
                       data.frame(createROCFrame(Res10Noise, 9, posControlList), method="D-index"),
                       data.frame(createROCFrame(Res10Noise, 10, posControlList), method="Distribution Specific Cut"))
Noise10PercROC <- roconMult2(Noise10PercDF, myTitle="ROC 10% Noise")

Noise25PercDF <- rbind(data.frame(createROCFrame(Res25Noise, 2, posControlList), method="K-Means"),
                       data.frame(createROCFrame(Res25Noise, 3, posControlList), method="Cox Regression"),
                       data.frame(createROCFrame(Res25Noise, 4, posControlList), method="Median"),
                       data.frame(createROCFrame(Res25Noise, 5, posControlList), method="Quantile 25th-75th"),
                       data.frame(createROCFrame(Res25Noise, 7, posControlList), method="Kaplan-Scan"),
                       data.frame(createROCFrame(Res25Noise, 8, posControlList), method="C-index"),
                       data.frame(createROCFrame(Res25Noise, 9, posControlList), method="D-index"),
                       data.frame(createROCFrame(Res25Noise, 10, posControlList), method="Distribution Specific Cut"))
Noise25PercROC <- roconMult2(Noise25PercDF, myTitle="ROC 25% Noise")

Noise50PercDF <- rbind(data.frame(createROCFrame(Res50Noise, 2, posControlList), method="K-Means"),
                       data.frame(createROCFrame(Res50Noise, 3, posControlList), method="Cox Regression"),
                       data.frame(createROCFrame(Res50Noise, 4, posControlList), method="Median"),
                       data.frame(createROCFrame(Res50Noise, 5, posControlList), method="Quantile 25th-75th"),
                       data.frame(createROCFrame(Res50Noise, 7, posControlList), method="Kaplan-Scan"),
                       data.frame(createROCFrame(Res50Noise, 8, posControlList), method="C-index"),
                       data.frame(createROCFrame(Res50Noise, 9, posControlList), method="D-index"),
                       data.frame(createROCFrame(Res50Noise, 10, posControlList), method="Distribution Specific Cut"))
Noise50PercROC <- roconMult2(Noise50PercDF, myTitle="ROC 50% Noise")

Noise75PercDF <- rbind(data.frame(createROCFrame(Res75Noise, 2, posControlList), method="K-Means"),
                       data.frame(createROCFrame(Res75Noise, 3, posControlList), method="Cox Regression"),
                       data.frame(createROCFrame(Res75Noise, 4, posControlList), method="Median"),
                       data.frame(createROCFrame(Res75Noise, 5, posControlList), method="Quantile 25th-75th"),
                       data.frame(createROCFrame(Res75Noise, 7, posControlList), method="Kaplan-Scan"),
                       data.frame(createROCFrame(Res75Noise, 8, posControlList), method="C-index"),
                       data.frame(createROCFrame(Res75Noise, 9, posControlList), method="D-index"),
                       data.frame(createROCFrame(Res75Noise, 10, posControlList), method="Distribution Specific Cut"))
Noise75PercROC <- roconMult2(Noise75PercDF, myTitle="ROC 70% Noise")

Noise100PercDF <- rbind(data.frame(createROCFrame(Res100Noise, 2, posControlList), method="K-Means"),
                        data.frame(createROCFrame(Res100Noise, 3, posControlList), method="Cox Regression"),
                        data.frame(createROCFrame(Res100Noise, 4, posControlList), method="Median"),
                        data.frame(createROCFrame(Res100Noise, 5, posControlList), method="Quantile 25th-75th"),
                        data.frame(createROCFrame(Res100Noise, 7, posControlList), method="Kaplan-Scan"),
                        data.frame(createROCFrame(Res100Noise, 8, posControlList), method="C-index"),
                        data.frame(createROCFrame(Res100Noise, 9, posControlList), method="D-index"),
                        data.frame(createROCFrame(Res100Noise, 10, posControlList), method="Distribution Specific Cut"))
Noise100PercROC <- roconMult2(Noise100PercDF, myTitle="ROC 100% Noise")

Noise150PercDF <- rbind(data.frame(createROCFrame(Res150Noise, 2, posControlList), method="K-Means"),
                        data.frame(createROCFrame(Res150Noise, 3, posControlList), method="Cox Regression"),
                        data.frame(createROCFrame(Res150Noise, 4, posControlList), method="Median"),
                        data.frame(createROCFrame(Res150Noise, 5, posControlList), method="Quantile 25th-75th"),
                        data.frame(createROCFrame(Res150Noise, 7, posControlList), method="Kaplan-Scan"),
                        data.frame(createROCFrame(Res150Noise, 8, posControlList), method="C-index"),
                        data.frame(createROCFrame(Res150Noise, 9, posControlList), method="D-index"),
                        data.frame(createROCFrame(Res150Noise, 10, posControlList), method="Distribution Specific Cut"))
Noise150PercROC <- roconMult2(Noise150PercDF, myTitle="ROC 100% Noise")


##################################################
# END Create ROC Objects
##################################################

# All ROC Plots
jpeg("../Figures/ROC_All_inSilico.jpg", width=4320, height=4320,  res=486)
grid.arrange(noNoiseROC[[2]]+ theme(legend.position="none")+ggtitle("No Noise"), Noise50PercROC[[2]]+ theme(legend.position="none")+ggtitle("50% Noise"), Noise100PercROC[[2]]+ theme(legend.position="none")+ggtitle("100% Noise"), Noise150PercROC[[2]]+ theme(legend.position="none")+ggtitle("150% Noise"), nrow=2)
dev.off()

#Table for All Plots
allAUCTable <- rbind(data.frame(noNoiseROC[[3]][,c(2,1)], NoiseLevel=0), 
                     data.frame(Noise10PercROC[[3]][,c(2,1)], NoiseLevel=.10),
                     data.frame(Noise25PercROC[[3]][,c(2,1)], NoiseLevel=.25),
                     data.frame(Noise50PercROC[[3]][,c(2,1)], NoiseLevel=.50),
                     data.frame(Noise75PercROC[[3]][,c(2,1)], NoiseLevel=.75),
                     data.frame(Noise100PercROC[[3]][,c(2,1)], NoiseLevel=1.0),
                     data.frame(Noise150PercROC[[3]][,c(2,1)], NoiseLevel=1.5))

jpeg("../Figures/AUC_VS_NOISE_inSilico.jpg", width=4320, height=4320,  res=486)
aucTable <- ggplot(allAUCTable, aes(NoiseLevel, AUCValue, color=method))+geom_point()+geom_smooth(se=F)+theme_bw() + ggtitle(" ")
aucTable
dev.off()

# for publication
nonoise <- noNoiseROC[[2]] + theme(legend.position="none") + ggtitle("No Noise") + theme(plot.margin = unit(c(0.7,0.2,0.3,0.2), "cm"))
perc50 <- Noise50PercROC[[2]] + theme(legend.position="none") + ggtitle("50% Noise") + theme(plot.margin = unit(c(0.7,0.2,0.3,0.2), "cm"))
perc100 <- Noise100PercROC[[2]] + theme(legend.position="none") + ggtitle("100% Noise") + theme(plot.margin = unit(c(0.4,0.2,0.2,0.2), "cm"))
perc150 <- Noise150PercROC[[2]] + theme(legend.position="none") + ggtitle("150% Noise") + theme(plot.margin = unit(c(0.4,0.2,0.2,0.2), "cm"))

jpeg("../Figures/Figure6.jpg", width=5500, height=3000,  res=486)
grid.arrange(nonoise, perc50, perc100, perc150, aucTable,
layout_matrix = rbind(c(1, 2, 5), c(3, 4, 5)), widths = c(1, 1, 2))
dev.off()

