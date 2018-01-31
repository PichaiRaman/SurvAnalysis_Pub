####################################################
# Author: Pichai Raman, Komal Rathi
# Function: Compare various methods of dealing with
# continuous variables during a survival analysis. 
# Date: 01/30/2018
####################################################

setwd('~/Projects/SurvAnalysis_Pub/code/')

# load libraries
library(survcomp)
library(reshape2)
library(GGally)
library(ggplot2)
library(lattice)
library(pheatmap)
library(mixtools)
library(snow)
library(survival)

# source R functions
source("../code/kmeansSA.R")
source("../code/KaplanScan.R")
source("../code/quantileCut.R")
source("../code/coxReg.R")
source("../code/helper.R")
source("../code/index.R")

# read data
print("Starting Read of Data")
load('../data/ALLDATA.RData')

# choose genes with most variability
geneCV <- apply(exprs_ov, FUN = max, MARGIN = 1)
geneCV <- sort(geneCV, T)

# start analysis
print("Starting Survival Analysis")
numGenes <- length(geneCV) # 20531

# cancers used
cancers <- c('hn','ov','ki','pr')

clus <- makeCluster(10)
print("Started cluster")

clusterExport(clus, "ov")
clusterExport(clus, "pr")
clusterExport(clus, "ki")
clusterExport(clus, "hn")
clusterExport(clus, "kmeansSA")
clusterExport(clus, "numGenes")
clusterExport(clus, "quantCutSA")
clusterExport(clus, "createSurvivalFrame")
clusterExport(clus, "kmScan")
clusterExport(clus, "kapmPlot")
clusterExport(clus, "Surv")
clusterExport(clus, "survdiff")
clusterExport(clus, "survfit")
clusterExport(clus, "coxph")
clusterExport(clus, "normalmixEM")
clusterExport(clus, "concordance.index")
clusterExport(clus, "D.index")

for(i in 1:length(cancers)){
  exprs <- paste0('exprs_', cancers[i])
  exprs <- get(exprs)
  annot <- paste0('annot_', cancers[i])
  annot <- get(annot)
  myData <- get(cancers[i])
  
  # run algorithms
  # kmeans
  fname <- paste0('../data/survAnalysisResults/kmeans_', cancers[i],'.txt')
  if(!file.exists(fname)){
    kmeans <- parSapply(clus, names(geneCV)[1:numGenes], FUN = kmeansSA, myData, tVar = "TimeVar", eVar = "eventVar")
    kmeans <- data.frame(t(data.frame(kmeans)))
    colnames(kmeans) <- c("Gene", "P.Value")
    write.table(kmeans, fname, sep = "\t", row.names = F)
    print(paste("Kmeans done for", cancers[i], sep = " "))
  } else {
    print("File exists")
  }
  
  # coxReg
  fname <- paste0('../data/survAnalysisResults/coxreg_', cancers[i], '.txt')
  if(!file.exists(fname)){
    cox.reg <- parSapply(clus, names(geneCV)[1:numGenes], FUN = coxReg, myData)
    cox.reg <- data.frame(t(data.frame(cox.reg)))
    colnames(cox.reg) <- c("Gene", "P.Value")
    write.table(cox.reg, fname, sep = "\t", row.names = F)
    print(paste("Cox regression done for", cancers[i], sep = " "))
  } else {
    print("File exists")
  }
  
  # Quantile Technique, cutting at median
  fname <- paste0('../data/survAnalysisResults/qCut50_', cancers[i], '.txt')
  if(!file.exists(fname)){
    qCut50 <- parSapply(clus, names(geneCV)[1:numGenes], FUN = quantCutSA, myData, F, quantLow = .50,  quantHigh = .50, tVar = "TimeVar", eVar = "eventVar")
    qCut50 <- data.frame(t(data.frame(qCut50)))
    colnames(qCut50) <- c("Gene", "P.Value")
    write.table(qCut50, fname, sep = "\t", row.names = F)
    print(paste("Quantile 50 done for", cancers[i], sep = " "))
  } else {
    print("File exists")
  }
  
  # Quantile Technique, cutting at 25 and 75
  fname <- paste0('../data/survAnalysisResults/qCut2575_', cancers[i], '.txt')
  if(!file.exists(fname)){
    qCut2575 <- parSapply(clus, names(geneCV)[1:numGenes], FUN = quantCutSA, myData, F, quantLow = .25,  quantHigh = .75, tVar = "TimeVar", eVar = "eventVar")
    qCut2575 <- data.frame(t(data.frame(qCut2575)))
    colnames(qCut2575) <- c("Gene", "P.Value")
    write.table(qCut2575, fname, sep = "\t", row.names = F)
    print(paste("Quantile 75 done for", cancers[i], sep = " "))
  } else {
    print("File exists")
  }
  
  # KM Scan Technique
  fname <- paste0('../data/survAnalysisResults/kmScan_', cancers[i], '.txt')
  if(!file.exists(fname)){
    kmScan <- parSapply(clus, names(geneCV)[1:numGenes], FUN = kapmPlot, myData, F, tVar = "TimeVar", eVar = "eventVar")
    kmScan <- data.frame(t(data.frame(kmScan)))
    colnames(kmScan) <- c("Gene", "P.Value", "Adj.P.Value")
    write.table(kmScan, fname, sep = "\t", row.names = F)
    print(paste("KM Scan done for",cancers[i],sep = " "))
  } else {
    print("File exists")
  }

  # c-index 
  fname <- paste0('../data/survAnalysisResults/cindex_', cancers[i], '.txt')
  if(!file.exists(fname)){
    cindex <- parSapply(clus, names(geneCV)[1:numGenes], FUN = c.index, myData)
    cindex <- data.frame(t(data.frame(cindex)))
    colnames(cindex) <- c("Gene", "P.Value")
    write.table(cindex, fname, sep = "\t", row.names = F)
    print(paste("c-index done for", cancers[i], sep = " "))
  } else {
    print("File exists")
  }
  
  # d-index 
  fname <- paste0('../data/survAnalysisResults/dindex_', cancers[i], '.txt')
  if(!file.exists(fname)){
    dindex <- parSapply(clus, names(geneCV)[1:numGenes], FUN = d.index, myData)
    dindex <- data.frame(t(data.frame(dindex)))
    colnames(dindex) <- c("Gene", "P.Value")
    write.table(dindex, fname, sep = "\t", row.names = F)
    print(paste("d-index done for", cancers[i], sep = " "))
  } else {
    print("File exists")
  }

  print(paste('All done for',cancers[i],sep = ' '))
}

stopCluster(clus)
print("Stopped cluster")

