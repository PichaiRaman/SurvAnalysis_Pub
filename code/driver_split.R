##################################################
# Author: Pichai Raman, Komal Rathi
# Function: Compare various methods of dealing with
# continuous variables during a survival analysis. 
# Date: 01/30/2018
##################################################

# call libraries
library(snow)
library(mixtools)
library(ggplot2)

##########################################
# source all code here
# methods of survival analysis 
source("../code/KaplanScan.R")
source("../code/quantileCut.R")
source("../code/coxReg.R")
source("../code/kmeansSA.R")
source("../code/index.R")

# some functions to help out 
source("../code/helper.R")

##########################################

##########################################
# read in data DONE
print("Starting Read of Data")
load("../data/ALLDATA.RData")

# fix colnames
# ovarian
colnames(exprs_ov) <- gsub("\\.", "-", colnames(exprs_ov))
rownames(exprs_ov) <- gsub(" ", "", rownames(exprs_ov))
ov <- list(exprs_ov, annot_ov)

# prostate
colnames(exprs_pr) <- gsub("\\.", "-", colnames(exprs_pr))
rownames(exprs_pr) <- gsub(" ", "", rownames(exprs_pr))
pr <- list(exprs_pr, annot_pr)

# head and neck
colnames(exprs_hn) <- gsub("\\.", "-", colnames(exprs_hn))
rownames(exprs_hn) <- gsub(" ", "", rownames(exprs_hn))
hn <- list(exprs_hn, annot_hn)

# kidney
colnames(exprs_ki) <- gsub("\\.", "-", colnames(exprs_ki))
rownames(exprs_ki) <- gsub(" ", "", rownames(exprs_ki))
ki <- list(exprs_ki, annot_ki)

print("Finished Reading Data")

##########################################
# split the datasets in two
cancers <- c('hn','ov','ki','pr')

split.file <- function(annot, exprs, cancer){
  len1 <- dim(annot[annot[,"eventVar"]==1,])[1]
  len0 <- dim(annot[annot[,"eventVar"]==0,])[1]
  an_1_p1 <- rownames(annot[annot[,"eventVar"]==1,])[1:round((len1/2))]
  an_1_p2 <- setdiff(rownames(annot[annot[,"eventVar"]==1,]), an_1_p1)
  an_0_p1 <- rownames(annot[annot[,"eventVar"]==0,])[1:round((len0/2))]
  an_0_p2 <- setdiff(rownames(annot[annot[,"eventVar"]==0,]), an_0_p1)
  an_1 <- annot[c(an_1_p1,an_0_p1),]
  an_2 <- annot[c(an_1_p2,an_0_p2),]
  exp_1 <- exprs[,c(an_1_p1,an_0_p1)]
  exp_2 <- exprs[,c(an_1_p2,an_0_p2)]
  # fname1 <- paste0("../data/splitData/", cancer, "ca_annot_1.txt") # no need to write it out
  # fname2 <- paste0("../data/splitData/", cancer, "ca_annot_2.txt")
  # fname3 <- paste0("../data/splitData/", cancer, "ca_exp_1.txt")
  # fname4 <- paste0("../data/splitData/", cancer, "ca_exp_2.txt")
  # write.table(an_1, fname1, sep="\t", row.names=T)
  # write.table(an_2, fname2, sep="\t", row.names=T)
  # write.table(exp_1, fname3, sep="\t", row.names=T)
  # write.table(exp_2, fname4, sep="\t", row.names=T)
  
  assign(x = paste0(cancer,"1"), list(exp_1, an_1), envir = globalenv())
  assign(x = paste0(cancer,"2"), list(exp_2, an_2), envir = globalenv())
}

split.file(annot = annot_ov, exprs = exprs_ov, cancer = "ov")
split.file(annot = annot_hn, exprs = exprs_hn, cancer = "hn")
split.file(annot = annot_ki, exprs = exprs_ki, cancer = "ki")
split.file(annot = annot_pr, exprs = exprs_pr, cancer = "pr")

##########################################
# choose genes with most variability 
geneCV <- apply(exprs_ov, FUN=max, MARGIN=1)
geneCV <- sort(geneCV, T)
numGenes <- length(geneCV)

##########################################
# apply functions
print("Starting Survival Analysis")

clus <- makeCluster(10)
print("Started cluster")
clusterExport(clus, "ov1")
clusterExport(clus, "pr1")
clusterExport(clus, "ki1")
clusterExport(clus, "hn1")
clusterExport(clus, "ov2")
clusterExport(clus, "pr2")
clusterExport(clus, "ki2")
clusterExport(clus, "hn2")
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
clusterExport(clus, "concordance.index")
clusterExport(clus, "D.index")


# 
cancers.split <- sort(apply(expand.grid(cancers, c(1:2)), 1, paste, collapse = "", sep = ""))

for(i in 1:length(cancers.split)){
  myData <- get(cancers.split[i])
  
  # run algorithms
  # kmeans
  fname <- paste0('../data/SurvAnalysisResultsRepro/kmeans_', cancers.split[i],'.txt')
  if(!file.exists(fname)){
    kmeans <- sapply(names(geneCV)[1:numGenes], FUN = kmeansSA, myData, tVar = "TimeVar", eVar = "eventVar")
    kmeans <- data.frame(t(data.frame(kmeans)))
    colnames(kmeans) <- c("Gene", "P.Value")
    write.table(kmeans, fname, sep = "\t", row.names = F)
    print(paste("Kmeans done for", cancers.split[i], sep = " "))
  } else {
    print("File exists")
  }
  
  # coxReg
  fname <- paste0('../data/SurvAnalysisResultsRepro/coxreg_', cancers.split[i], '.txt')
  if(!file.exists(fname)){
    cox.reg <- sapply(names(geneCV)[1:numGenes], FUN = coxReg, myData)
    cox.reg <- data.frame(t(data.frame(cox.reg)))
    colnames(cox.reg) <- c("Gene", "P.Value")
    write.table(cox.reg, fname, sep = "\t", row.names = F)
    print(paste("Cox regression done for", cancers.split[i], sep = " "))
  } else {
    print("File exists")
  }
  
  # Quantile Technique, cutting at median
  fname <- paste0('../data/SurvAnalysisResultsRepro/qCut50_', cancers.split[i], '.txt')
  if(!file.exists(fname)){
    qCut50 <- sapply(names(geneCV)[1:numGenes], FUN = quantCutSA, myData, F, quantLow = .50,  quantHigh = .50, tVar = "TimeVar", eVar = "eventVar")
    qCut50 <- data.frame(t(data.frame(qCut50)))
    colnames(qCut50) <- c("Gene", "P.Value")
    write.table(qCut50, fname, sep = "\t", row.names = F)
    print(paste("Quantile 50 done for", cancers.split[i], sep = " "))
  } else {
    print("File exists")
  }
  
  # Quantile Technique, cutting at 25 and 75
  fname <- paste0('../data/SurvAnalysisResultsRepro/qCut2575_', cancers.split[i], '.txt')
  if(!file.exists(fname)){
    qCut2575 <- sapply(names(geneCV)[1:numGenes], FUN = quantCutSA, myData, F, quantLow = .25,  quantHigh = .75, tVar = "TimeVar", eVar = "eventVar")
    qCut2575 <- data.frame(t(data.frame(qCut2575)))
    colnames(qCut2575) <- c("Gene", "P.Value")
    write.table(qCut2575, fname, sep = "\t", row.names = F)
    print(paste("Quantile 75 done for", cancers.split[i], sep = " "))
  } else {
    print("File exists")
  }
  
  # KM Scan Technique
  fname <- paste0('../data/SurvAnalysisResultsRepro/kmScan_', cancers.split[i], '.txt')
  if(!file.exists(fname)){
    kmScan <- sapply(names(geneCV)[1:numGenes], FUN = kapmPlot, myData, F, tVar = "TimeVar", eVar = "eventVar")
    kmScan <- data.frame(t(data.frame(kmScan)))
    colnames(kmScan) <- c("Gene", "P.Value", "Adj.P.Value")
    write.table(kmScan, fname, sep = "\t", row.names = F)
    print(paste("KM Scan done for",cancers.split[i],sep = " "))
  } else {
    print("File exists")
  }
  
  # c-index 
  fname <- paste0('../data/SurvAnalysisResultsRepro/cindex_', cancers.split[i], '.txt')
  if(!file.exists(fname)){
    cindex <- sapply(names(geneCV)[1:numGenes], FUN = c.index, myData)
    cindex <- data.frame(t(data.frame(cindex)))
    colnames(cindex) <- c("Gene", "P.Value")
    write.table(cindex, fname, sep = "\t", row.names = F)
    print(paste("c-index done for", cancers.split[i], sep = " "))
  } else {
    print("File exists")
  }
  
  # d-index 
  fname <- paste0('../data/SurvAnalysisResultsRepro/dindex_', cancers.split[i], '.txt')
  if(!file.exists(fname)){
    dindex <- sapply(names(geneCV)[1:numGenes], FUN = d.index, myData)
    dindex <- data.frame(t(data.frame(dindex)))
    colnames(dindex) <- c("Gene", "P.Value")
    write.table(dindex, fname, sep = "\t", row.names = F)
    print(paste("d-index done for", cancers.split[i], sep = " "))
  } else {
    print("File exists")
  }
  
  print(paste('All done for',cancers.split[i],sep = ' '))
}

stopCluster(clus)
print("Stopped cluster")


