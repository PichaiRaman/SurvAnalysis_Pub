##############################################
#Code to simulate RNA-Seq data and see
#how well methods can deal with noise 
#in the system
#
##############################################

#Call libraries
library("SimSeq");

##########################################
#Source all code here
#Methods of survival analysis 
source("../code/KaplanScan.R");
source("../code/quantileCut.R");
source("../code/coxReg.R");
source("../code/kmeansSA.R");

#Some functions to help out 
source("../code/helper.R");

##########################################


#########################################
#Read in data DONE
load("../data/ALLDATA.RData");

#########################################
#Generate simulated data
simDat <- SimData(counts = round(exprs_hn), treatment=annot_hn[,"eventVar"], norm.factors=2*(annot_hn[,"eventVar"]+.5), n.genes=5000, n.diff=250, k.ind=100, sort.method="unpaired");

#Counts
simExprs <- data.frame(simDat$counts);

#metaData
simMeta <- data.frame(colnames(simExprs), simDat$treatment);
set.seed(100);
#survTimes <- c(sample(annot_hn[annot_hn["eventVar"]==0, "TimeVar"], 100), sample(annot_hn[annot_hn["eventVar"]==1, "TimeVar"], 100));

survTimes1 <- sort(annot_hn[annot_hn["eventVar"]==0, "TimeVar"])[(nrow(annot_hn[annot_hn["eventVar"]==0,])-99): nrow(annot_hn[annot_hn["eventVar"]==0,])];
survTimes2 <- sort(annot_hn[annot_hn["eventVar"]==1, "TimeVar"])[c(1:100)];
survTimes <- c(survTimes1, survTimes2);

simMeta[,"TimeVar"] <- survTimes;
colnames(simMeta)[1:2] <- c("Sample", "eventVar");

simObj <- list(simExprs, simMeta);

numGenes <- 10;

CreateMatrix <- function(simObj)
{
kmeans_sim <- sapply(rownames(simExprs), FUN= kmeansSA, simObj, tVar="TimeVar", eVar="eventVar");
kmeans_sim <- data.frame(t(data.frame(kmeans_sim)));
colnames(kmeans_sim) <- c("Gene", "P.Value");

coxReg_sim <- sapply(rownames(simExprs), FUN= coxReg, simObj);
coxReg_sim <- data.frame(t(data.frame(coxReg_sim)));
colnames(coxReg_sim) <- c("Gene", "P.Value");

qCut50_sim <- sapply(rownames(simExprs), FUN=quantCutSA, simObj, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_sim <- data.frame(t(data.frame(qCut50_sim)));
colnames(qCut50_sim) <- c("Gene", "P.Value");

qCut2575_sim <- sapply(rownames(simExprs), FUN=quantCutSA, simObj, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_sim <- data.frame(t(data.frame(qCut2575_sim)));
colnames(qCut2575_sim) <- c("Gene", "P.Value");

kmScan_sim <- sapply(rownames(simExprs), FUN=kapmPlot, simObj, F, tVar="TimeVar", eVar="eventVar");
kmScan_sim <- data.frame(t(data.frame(kmScan_sim)));
colnames(kmScan_sim) <- c("Gene", "P.Value", "Adj.P.Value");

coxReg_sim <- coxReg_sim[rownames(kmeans_sim),];
qCut50_sim <- qCut50_sim[rownames(kmeans_sim),];
qCut2575_sim <- qCut2575_sim[rownames(kmeans_sim),];
kmScan_sim <- kmScan_sim[rownames(kmeans_sim),];

myDF <- data.frame(kmeans_sim, coxReg_sim[2], qCut50_sim[2], qCut2575_sim[2], kmScan_sim[2:3]);
colnames(myDF) <- c("gene", "kmeans", "coxreg", "qcut50", "qcut2575", "kmScanP", "kmScanQ"); 
return(myDF)
}











