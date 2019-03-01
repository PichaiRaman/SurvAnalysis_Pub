##############################################
# Author: Pichai Raman, Komal Rathi
# Function: Code to simulate RNA-Seq data and see
# how well methods can deal with noise 
# in the system
# Date: 01/30/2018
##############################################

# call libraries
library(SimSeq)
source("~/Projects/Survival_Analysis_RNASeq/code/rocon.R")

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


#########################################
# read in data DONE
load("../data/ALLDATA.RData")

#########################################
# generate simulated data
set.seed(50)
simDat <- SimData(counts = round(as.matrix(exprs_hn)), treatment=round(rnorm(520, mean=.5, sd=.01)), n.genes=5000, n.diff=0, k.ind=80, sort.method="unpaired")

# counts
simExprs <- data.frame(simDat$counts)
simExprs <- simExprs[,1:100]

# metaData
simMeta <- data.frame(colnames(simExprs), simDat$treatment[1:100])
set.seed(50)

survTimes1 <- sort(annot_ki[annot_ki["eventVar"]==0, "TimeVar"])[(nrow(annot_ki[annot_ki["eventVar"]==0,])-79): nrow(annot_ki[annot_ki["eventVar"]==0,])]
survTimes2 <- sort(annot_ki[annot_ki["eventVar"]==1, "TimeVar"])[c(1:20)]
survTimes <- c(survTimes1, survTimes2)

simMeta[,"TimeVar"] <- survTimes
colnames(simMeta)[1:2] <- c("Sample", "eventVar")

# now Let's generate positive controls
genPosControlMat <- function(x) {
  tmpMean <- mean(x)
  numUp <- abs(round(rnorm(1, mean=16, sd=6)))
  survTimesNorm <- abs(c(rnorm((100-numUp), mean=.1, sd=.01),rnorm(numUp, mean=1,sd=.01)))
  multiplier <- runif(1, min=1, max=8)
  xN <- round(multiplier*x*survTimesNorm)
}

# get genes with only values > 100
posControlGenes <- rownames(simExprs)%in%sample(rownames(simExprs)[apply(simExprs, FUN=mean, MARGIN=1)>100], 250)
simExprs[posControlGenes,] <- simExprs[posControlGenes,]+t(apply(simExprs[posControlGenes,], FUN=genPosControlMat, MARGIN=1))
simObjNoNoise <- list(simExprs, simMeta)

addNoise <- function(myPerc) {
  myMeans <-as.numeric(apply(simExprs, FUN=mean, MARGIN=1))*myPerc
  mySD <-as.numeric(apply(simExprs, FUN=sd, MARGIN=1))*myPerc
  
  myDat <- rnorm(ncol(simExprs), mean=myMeans[1], sd=mySD[1])
  for(i in 2:length(myMeans))
  {
    tmpMat <- rnorm(ncol(simExprs), mean=myMeans[i], sd=mySD[i])
    myDat <- rbind(myDat, tmpMat)
  }
  return(as.matrix(myDat))
}


CreateMatrix <- function(simObj) {
  kmeans_sim <- sapply(rownames(simExprs), FUN= kmeansSA, simObj, tVar="TimeVar", eVar="eventVar")
  kmeans_sim <- data.frame(t(data.frame(kmeans_sim)))
  colnames(kmeans_sim) <- c("Gene", "P.Value")
  print("done kmeans")
  
  coxReg_sim <- sapply(rownames(simExprs), FUN= coxReg, simObj)
  coxReg_sim <- data.frame(t(data.frame(coxReg_sim)))
  colnames(coxReg_sim) <- c("Gene", "P.Value")
  print("done cox")
  
  qCut50_sim <- sapply(rownames(simExprs), FUN=quantCutSA, simObj, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar")
  qCut50_sim <- data.frame(t(data.frame(qCut50_sim)))
  colnames(qCut50_sim) <- c("Gene", "P.Value")
  print("done median cut")
  
  qCut2575_sim <- sapply(rownames(simExprs), FUN=quantCutSA, simObj, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar")
  qCut2575_sim <- data.frame(t(data.frame(qCut2575_sim)))
  colnames(qCut2575_sim) <- c("Gene", "P.Value")
  print("done quantile cut")
  
  kmScan_sim <- sapply(rownames(simExprs), FUN=kapmPlot, simObj, F, tVar="TimeVar", eVar="eventVar")
  kmScan_sim <- data.frame(t(data.frame(kmScan_sim)))
  colnames(kmScan_sim) <- c("Gene", "P.Value", "Adj.P.Value")
  print("done kmScan")
  
  cindex_sim <- sapply(rownames(simExprs), FUN=c.index, simObj)
  cindex_sim <- data.frame(t(data.frame(cindex_sim)))
  colnames(cindex_sim) <- c("Gene", "P.Value")
  print("done c-index")
  
  dindex_sim <- sapply(rownames(simExprs), FUN=d.index, simObj)
  dindex_sim <- data.frame(t(data.frame(dindex_sim)))
  colnames(dindex_sim) <- c("Gene", "P.Value")
  print("done d-index")
  
  coxReg_sim <- coxReg_sim[rownames(kmeans_sim),]
  qCut50_sim <- qCut50_sim[rownames(kmeans_sim),]
  qCut2575_sim <- qCut2575_sim[rownames(kmeans_sim),]
  kmScan_sim <- kmScan_sim[rownames(kmeans_sim),]
  cindex_sim <- cindex_sim[rownames(kmeans_sim),]
  dindex_sim <- dindex_sim[rownames(kmeans_sim),]
  
  myDF <- data.frame(kmeans_sim, coxReg_sim[2], qCut50_sim[2], qCut2575_sim[2], kmScan_sim[2:3], cindex_sim[2], dindex_sim[2])
  colnames(myDF) <- c("gene", "kmeans", "coxreg", "qcut50", "qcut2575", "kmScanP", "kmScanQ", "c-index", "d-index") 
  
  return(myDF)
}

# ROC frame
createROCFrame <- function(data, myCol, cList) {
  data <- data[myCol]
  data[,1] <- (-1)*log10(data[,1])
  data[,"gene"] <- rownames(data)
  data[,"symbol"] <- rownames(data)
  data[,"labs"]<- as.numeric(data[,"gene"] %in% cList)
  data <- data[,c(1,4)]
  colnames(data)[1] <- "preds"
  return(data)
}

# write no noise
write.table(data.frame(posControlGenes), "../data/SyntheticData/PositiveControls.txt", sep="\t", row.names=F)
posControlList <- gsub("\\|", ".", rownames(simExprs)[posControlGenes])

iter <- c(0,.1,.25,.5,.75, 1, 1.5, 2, 3, 4, 5)
for(i in 1:length(iter)) {
  resTitle <- paste("../data/SyntheticData/Res_", iter[i], "Noise.txt", sep="")
  if(!file.exists(resTitle)){
    simExprsTmp <- round(abs(addNoise(iter[i])+simExprs)) 
    simObjTmpNoise <- list(simExprsTmp, simMeta)
    ResTmpNoise <- CreateMatrix(simObjTmpNoise)
    write.table(ResTmpNoise, resTitle, sep="\t", row.names=T)
    print("done no noise") 
  }
}

