##########################################
#Objective : Compare various methods of dealing with
#continuous variables during a survival analysis. 
#
##########################################


#Call libraries
library("GGally");
library("ggplot2");
library("lattice");
library("pheatmap");
library("mixtools");
library("snow");


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

##########################################
#Read in data DONE
print("Starting Read of Data");

#ovarian
annot_ov <- read.delim("../data/raw/ovca/annot.txt");
exprs_ov <- read.delim("../data/raw/ovca/exprs.txt")
ov <- list(exprs_ov, annot_ov);

#prostate
annot_pr <- read.delim("../data/raw/prca/annot.txt");
exprs_pr <- read.delim("../data/raw/prca/exprs.txt")
pr <- list(exprs_pr, annot_pr);

#head and neck
annot_hn <- read.delim("../data/raw/hnca/annot.txt");
exprs_hn <- read.delim("../data/raw/hnca/exprs.txt")
hn <- list(exprs_hn, annot_hn);

#Kidney
annot_ki <- read.delim("../data/raw/kica/annot.txt");
exprs_ki <- read.delim("../data/raw/kica/exprs.txt")
ki <- list(exprs_ki, annot_ki);
##########################################

##########################################
#Choose genes with most variability DONE
geneCV <- apply(exprs_ov, FUN=max, MARGIN=1);
geneCV <- sort(geneCV, T);

##########################################

##########################################
#Run Surival analysis using various techniques DONE 
#
#numGenes here is the number of genes we want to run survival analysis on, we should probably
#do all genes but kmScan takes a bit of time so for testing purposes let's set n to a small number
#
print("Starting Survival Analysis");
numGenes <- 20531;
clus <- makeCluster(20);
print("Started cluster");

clusterExport(clus, "ov");
clusterExport(clus, "pr");
clusterExport(clus, "ki");
clusterExport(clus, "hn");
clusterExport(clus, "kmeansSA");
clusterExport(clus, "numGenes");
clusterExport(clus, "quantCutSA");
clusterExport(clus, "createSurvivalFrame");
clusterExport(clus, "kmScan");
clusterExport(clus, "kapmPlot");
clusterExport(clus, "Surv");
clusterExport(clus, "survdiff");
clusterExport(clus, "survfit");
clusterExport(clus, "coxph");
clusterExport(clus, "normalmixEM");



#Cox Regression
coxReg_pr <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, pr);
coxReg_pr <- data.frame(t(data.frame(coxReg_pr)));
colnames(coxReg_pr) <- c("Gene", "P.Value");
print("Done Prostate");

coxReg_ov <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, ov);
coxReg_ov <- data.frame(t(data.frame(coxReg_ov)));
colnames(coxReg_ov) <- c("Gene", "P.Value");
print("Done Ovarian");

coxReg_ki <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, ki);
coxReg_ki <- data.frame(t(data.frame(coxReg_ki)));
colnames(coxReg_ki) <- c("Gene", "P.Value");
print("Done Kidney");

coxReg_hn <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, hn);
coxReg_hn <- data.frame(t(data.frame(coxReg_hn)));
colnames(coxReg_hn) <- c("Gene", "P.Value");
print("Done Head and Neck");
print("Finished cox regression");


#Quantile Technique, cutting at median
qCut50_ov <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ov, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_ov <- data.frame(t(data.frame(qCut50_ov)));
colnames(qCut50_ov) <- c("Gene", "P.Value");
print("Done Ovarian");

qCut50_pr <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, pr, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_pr <- data.frame(t(data.frame(qCut50_pr)));
colnames(qCut50_pr) <- c("Gene", "P.Value");
print("Done Prostate");

qCut50_ki <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ki, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_ki <- data.frame(t(data.frame(qCut50_ki)));
colnames(qCut50_ki) <- c("Gene", "P.Value");
print("Done Kidney");

qCut50_hn <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, hn, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_hn <- data.frame(t(data.frame(qCut50_hn)));
colnames(qCut50_hn) <- c("Gene", "P.Value");
print("Done Head and Neck");
print("Finished median quantile cut");


#Quantile Technique, cutting at 25 and 75
qCut2575_ov <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ov, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_ov <- data.frame(t(data.frame(qCut2575_ov)));
colnames(qCut2575_ov) <- c("Gene", "P.Value");
print("Done Ovarian");

qCut2575_pr <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, pr, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_pr <- data.frame(t(data.frame(qCut2575_pr)));
colnames(qCut2575_pr) <- c("Gene", "P.Value");
print("Done Prostate");

qCut2575_ki <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ki, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_ki <- data.frame(t(data.frame(qCut2575_ki)));
colnames(qCut2575_ki) <- c("Gene", "P.Value");
print("Done Kidney");

qCut2575_hn <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, hn, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_hn <- data.frame(t(data.frame(qCut2575_hn)));
colnames(qCut2575_hn) <- c("Gene", "P.Value");
print("Done Head and Neck");
print("Finished 75th 25th quantile cut");

#Kmneans
kmeans_ov <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, ov, tVar="TimeVar", eVar="eventVar");
kmeans_ov <- data.frame(t(data.frame(kmeans_ov)));
colnames(kmeans_ov) <- c("Gene", "P.Value");
print("Done Ovarian");

kmeans_pr <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, pr, tVar="TimeVar", eVar="eventVar");
kmeans_pr <- data.frame(t(data.frame(kmeans_pr)));
colnames(kmeans_pr) <- c("Gene", "P.Value");
print("Done Prostate");

kmeans_ki <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, ki, tVar="TimeVar", eVar="eventVar");
kmeans_ki <- data.frame(t(data.frame(kmeans_ki)));
colnames(kmeans_ki) <- c("Gene", "P.Value");
print("Done Kidney");

kmeans_hn <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, hn, tVar="TimeVar", eVar="eventVar");
kmeans_hn <- data.frame(t(data.frame(kmeans_hn)));
colnames(kmeans_hn) <- c("Gene", "P.Value");

print("Finished GMM");

#KM Scan Technique
kmScan_ov <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, ov, F, tVar="TimeVar", eVar="eventVar");
kmScan_ov <- data.frame(t(data.frame(kmScan_ov)));
colnames(kmScan_ov) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Ovarian");

kmScan_pr <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, pr, F, tVar="TimeVar", eVar="eventVar");
kmScan_pr <- data.frame(t(data.frame(kmScan_pr)));
colnames(kmScan_pr) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Prostate");

kmScan_ki <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, ki, F, tVar="TimeVar", eVar="eventVar");
kmScan_ki <- data.frame(t(data.frame(kmScan_ki)));
colnames(kmScan_ov) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Kidney");

kmScan_hn <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, hn, F, tVar="TimeVar", eVar="eventVar");
kmScan_hn <- data.frame(t(data.frame(kmScan_hn)));
colnames(kmScan_hn) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Head and Neck");
print("Finished KM Scan");

stopCluster(clus);
print("Stopped cluster");

##########################################


##########################################
#Merge all results together into one data frame and a matrix for convenience DONE
print("Starting Data Merge");

results_ov <- cbind(kmScan_ov[3], qCut50_ov[2], qCut2575_ov[2], coxReg_ov[2], kmeans_ov[2]);
colnames(results_ov) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value", "Cox.P.Value", "kmeans.P.Value")
resultsMat_ov <- sapply(results_ov, FUN=factToNum);
rownames(resultsMat_ov) <- rownames(resultsMat_ov);
write.table(resultsMat_ov, "allresults_ov.txt", sep="\t", row.names=T);

results_pr <- cbind(kmScan_pr[3], qCut50_pr[2], qCut2575_pr[2], coxReg_pr[2], kmeans_pr[2]);
colnames(results_pr) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value","Cox.P.Value", "kmeans.P.Value")
resultsMat_pr <- sapply(results_pr, FUN=factToNum);
rownames(resultsMat_pr) <- rownames(resultsMat_pr);
write.table(resultsMat_pr, "allresults_pr.txt", sep="\t", row.names=T);

results_ki <- cbind(kmScan_ki[3], qCut50_ki[2], qCut2575_ki[2], coxReg_ki[2], kmeans_ki[2]);
colnames(results_ki) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value", "Cox.P.Value", "kmeans.P.Value")
resultsMat_ki <- sapply(results_ki, FUN=factToNum);
rownames(resultsMat_ki) <- rownames(resultsMat_ki);
write.table(resultsMat_ki, "allresults_ki.txt", sep="\t", row.names=T);

results_hn <- cbind(kmScan_hn[3], qCut50_hn[2], qCut2575_hn[2], coxReg_hn[2], kmeans_hn[2]);
colnames(results_hn) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value", "Cox.P.Value", "kmeans.P.Value")
resultsMat_hn <- sapply(results_hn, FUN=factToNum);
rownames(resultsMat_hn) <- rownames(resultsMat_hn);
write.table(resultsMat_hn, "allresults_hn.txt", sep="\t", row.names=T);

##########################################

