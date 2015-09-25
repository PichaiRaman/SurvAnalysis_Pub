
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
#Create two datasets from 1

#Ovarian
len1 <- dim(annot_ov[annot_ov[,"eventVar"]==1,])[1];
len0 <- dim(annot_ov[annot_ov[,"eventVar"]==0,])[1];
ov_1 <- sample(rownames(annot_ov[annot_ov[,"eventVar"]==1,]), len1/2);
ov_2 <- sample(rownames(annot_ov[annot_ov[,"eventVar"]==0,]), len0/2);
ov_an_1 <- annot_ov[ov_an_1,]
ov_an_2 <- annot_ov[ov_an_2,];
ov_exp_1 <- exprs_ov[,ov_an_1];
ov_exp_2 <- exprs_ov[,ov_an_2];
ov_1 <- list(ov_exp_1, ov_an_1);
ov_2 <- list(ov_exp_2, ov_an_2);

#Prostate
len1 <- dim(annot_pr[annot_pr[,"eventVar"]==1,])[1];
len0 <- dim(annot_pr[annot_pr[,"eventVar"]==0,])[1];
pr_1 <- sample(rownames(annot_pr[annot_pr[,"eventVar"]==1,]), len1/2);
pr_2 <- sample(rownames(annot_pr[annot_pr[,"eventVar"]==0,]), len0/2);
pr_an_1 <- annot_pr[pr_an_1,]
pr_an_2 <- annot_pr[pr_an_2,];
pr_exp_1 <- exprs_pr[,pr_an_1];
pr_exp_2 <- exprs_pr[,pr_an_2];
pr_1 <- list(pr_exp_1, pr_an_1);
pr_2 <- list(pr_exp_2, pr_an_2);

#Kidney
len1 <- dim(annot_ki[annot_ki[,"eventVar"]==1,])[1];
len0 <- dim(annot_ki[annot_ki[,"eventVar"]==0,])[1];
ki_1 <- sample(rownames(annot_ki[annot_ki[,"eventVar"]==1,]), len1/2);
ki_2 <- sample(rownames(annot_ki[annot_ki[,"eventVar"]==0,]), len0/2);
ki_an_1 <- annot_ki[ki_an_1,]
ki_an_2 <- annot_ki[ki_an_2,];
ki_exp_1 <- exprs_ki[,ki_an_1];
ki_exp_2 <- exprs_ki[,ki_an_2];
ki_1 <- list(ki_exp_1, ki_an_1);
ki_2 <- list(ki_exp_2, ki_an_2);

#Head and Neck
len1 <- dim(annot_hn[annot_hn[,"eventVar"]==1,])[1];
len0 <- dim(annot_hn[annot_hn[,"eventVar"]==0,])[1];
hn_1 <- sample(rownames(annot_hn[annot_hn[,"eventVar"]==1,]), len1/2);
hn_2 <- sample(rownames(annot_hn[annot_hn[,"eventVar"]==0,]), len0/2);
hn_an_1 <- annot_hn[hn_an_1,]
hn_an_2 <- annot_hn[hn_an_2,];
hn_exp_1 <- exprs_hn[,hn_an_1];
hn_exp_2 <- exprs_hn[,hn_an_2];
hn_1 <- list(hn_exp_1, hn_an_1);
hn_2 <- list(hn_exp_2, hn_an_2);






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
clus <- makeCluster(10);
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

#Kmneans
kmeans_ov1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, ov1, tVar="TimeVar", eVar="eventVar");
kmeans_ov1 <- data.frame(t(data.frame(kmeans_ov1)));
colnames(kmeans_ov1) <- c("Gene", "P.Value");
print("Done Ovarian");

kmeans_pr1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, pr1, tVar="TimeVar", eVar="eventVar");
kmeans_pr1 <- data.frame(t(data.frame(kmeans_pr1)));
colnames(kmeans_pr1) <- c("Gene", "P.Value");
print("Done Prostate");

kmeans_ki1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, ki1, tVar="TimeVar", eVar="eventVar");
kmeans_ki1 <- data.frame(t(data.frame(kmeans_ki1)));
colnames(kmeans_ki1) <- c("Gene", "P.Value");
print("Done Kidney");

kmeans_hn1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, hn1, tVar="TimeVar", eVar="eventVar");
kmeans_hn1 <- data.frame(t(data.frame(kmeans_hn1)));
colnames(kmeans_hn1) <- c("Gene", "P.Value");

print("Finished KMeans 1");

#Kmneans
kmeans_ov2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, ov2, tVar="TimeVar", eVar="eventVar");
kmeans_ov2 <- data.frame(t(data.frame(kmeans_ov2)));
colnames(kmeans_ov2) <- c("Gene", "P.Value");
print("Done Ovarian");

kmeans_pr2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, pr2, tVar="TimeVar", eVar="eventVar");
kmeans_pr2 <- data.frame(t(data.frame(kmeans_pr2)));
colnames(kmeans_pr2) <- c("Gene", "P.Value");
print("Done Prostate");

kmeans_ki2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, ki2, tVar="TimeVar", eVar="eventVar");
kmeans_ki2 <- data.frame(t(data.frame(kmeans_ki2)));
colnames(kmeans_ki2) <- c("Gene", "P.Value");
print("Done Kidney");

kmeans_hn2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, hn2, tVar="TimeVar", eVar="eventVar");
kmeans_hn2 <- data.frame(t(data.frame(kmeans_hn2)));
colnames(kmeans_hn2) <- c("Gene", "P.Value");

print("Finished KMeans 2");



#Cox Regression
coxReg_pr1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, pr1);
coxReg_pr1 <- data.frame(t(data.frame(coxReg_pr1)));
colnames(coxReg_pr1) <- c("Gene", "P.Value");
print("Done Prostate");

coxReg_ov1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, ov1);
coxReg_ov1 <- data.frame(t(data.frame(coxReg_ov1)));
colnames(coxReg_ov1) <- c("Gene", "P.Value");
print("Done Ovarian");

coxReg_ki1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, ki1);
coxReg_ki1 <- data.frame(t(data.frame(coxReg_ki1)));
colnames(coxReg_ki1) <- c("Gene", "P.Value");
print("Done Kidney");

coxReg_hn1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, hn1);
coxReg_hn1 <- data.frame(t(data.frame(coxReg_hn1)));
colnames(coxReg_hn1) <- c("Gene", "P.Value");
print("Done Head and Neck");
print("Finished cox regression 1");

coxReg_pr2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, pr2);
coxReg_pr2 <- data.frame(t(data.frame(coxReg_pr2)));
colnames(coxReg_pr2) <- c("Gene", "P.Value");
print("Done Prostate");

coxReg_ov2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, ov2);
coxReg_ov2 <- data.frame(t(data.frame(coxReg_ov2)));
colnames(coxReg_ov2) <- c("Gene", "P.Value");
print("Done Ovarian");

coxReg_ki2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, ki2);
coxReg_ki2 <- data.frame(t(data.frame(coxReg_ki2)));
colnames(coxReg_ki2) <- c("Gene", "P.Value");
print("Done Kidney");

coxReg_hn2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, hn)2;
coxReg_hn2 <- data.frame(t(data.frame(coxReg_hn2)));
colnames(coxReg_hn2) <- c("Gene", "P.Value");
print("Done Head and Neck");
print("Finished cox regression 2");


#Quantile Technique, cutting at median
qCut50_ov1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ov1, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_ov1 <- data.frame(t(data.frame(qCut50_ov1)));
colnames(qCut50_ov1) <- c("Gene", "P.Value");
print("Done Ovarian");

qCut50_pr1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, pr1, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_pr1 <- data.frame(t(data.frame(qCut50_pr1)));
colnames(qCut50_pr1) <- c("Gene", "P.Value");
print("Done Prostate");

qCut50_ki1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ki1, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_ki1 <- data.frame(t(data.frame(qCut50_ki1)));
colnames(qCut50_ki1) <- c("Gene", "P.Value");
print("Done Kidney");

qCut50_hn1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, hn1, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_hn1 <- data.frame(t(data.frame(qCut50_hn1)));
colnames(qCut50_hn1) <- c("Gene", "P.Value");
print("Done Head and Neck");
print("Finished median quantile cut 1");


qCut50_ov2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ov2, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_ov2 <- data.frame(t(data.frame(qCut50_ov2)));
colnames(qCut50_ov2) <- c("Gene", "P.Value");
print("Done Ovarian");

qCut50_pr2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, pr2, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_pr2 <- data.frame(t(data.frame(qCut50_pr2)));
colnames(qCut50_pr2) <- c("Gene", "P.Value");
print("Done Prostate");

qCut50_ki2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ki2, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_ki2 <- data.frame(t(data.frame(qCut50_ki2)));
colnames(qCut50_ki2) <- c("Gene", "P.Value");
print("Done Kidney");

qCut50_hn2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, hn2, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_hn2 <- data.frame(t(data.frame(qCut50_hn2)));
colnames(qCut50_hn2) <- c("Gene", "P.Value");
print("Done Head and Neck");
print("Finished median quantile cut 2");



#Quantile Technique, cutting at 25 and 75
qCut2575_ov1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ov1, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_ov1 <- data.frame(t(data.frame(qCut2575_ov1)));
colnames(qCut2575_ov1) <- c("Gene", "P.Value");
print("Done Ovarian");

qCut2575_pr1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, pr1, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_pr1 <- data.frame(t(data.frame(qCut2575_pr1)));
colnames(qCut2575_pr1) <- c("Gene", "P.Value");
print("Done Prostate");

qCut2575_ki1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ki1, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_ki1 <- data.frame(t(data.frame(qCut2575_ki1)));
colnames(qCut2575_ki1) <- c("Gene", "P.Value");
print("Done Kidney");

qCut2575_hn1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, hn1, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_hn1 <- data.frame(t(data.frame(qCut2575_hn1)));
colnames(qCut2575_hn1) <- c("Gene", "P.Value");
print("Done Head and Neck");
print("Finished 75th 25th quantile cut 1");

qCut2575_ov2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ov2, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_ov2 <- data.frame(t(data.frame(qCut2575_ov2)));
colnames(qCut2575_ov2) <- c("Gene", "P.Value");
print("Done Ovarian");

qCut2575_pr2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, pr2, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_pr2 <- data.frame(t(data.frame(qCut2575_pr2)));
colnames(qCut2575_pr2) <- c("Gene", "P.Value");
print("Done Prostate");

qCut2575_ki2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ki2, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_ki2 <- data.frame(t(data.frame(qCut2575_ki2)));
colnames(qCut2575_ki2) <- c("Gene", "P.Value");
print("Done Kidney");

qCut2575_hn2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, hn2, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_hn2 <- data.frame(t(data.frame(qCut2575_hn2)));
colnames(qCut2575_hn) <- c("Gene", "P.Value");
print("Done Head and Neck");
print("Finished 75th 25th quantile cut 2");



#KM Scan Technique
kmScan_ov1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, ov1, F, tVar="TimeVar", eVar="eventVar");
kmScan_ov1 <- data.frame(t(data.frame(kmScan_ov1)));
colnames(kmScan_ov1) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Ovarian");

kmScan_pr1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, pr1, F, tVar="TimeVar", eVar="eventVar");
kmScan_pr1 <- data.frame(t(data.frame(kmScan_pr1)));
colnames(kmScan_pr1) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Prostate");

kmScan_ki1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, ki1, F, tVar="TimeVar", eVar="eventVar");
kmScan_ki1 <- data.frame(t(data.frame(kmScan_ki1)));
colnames(kmScan_ki1) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Kidney");

kmScan_hn1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, hn1, F, tVar="TimeVar", eVar="eventVar");
kmScan_hn1 <- data.frame(t(data.frame(kmScan_hn1)));
colnames(kmScan_hn1) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Head and Neck");
print("Finished KM Scan 1");


#KM Scan Technique
kmScan_ov2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, ov2, F, tVar="TimeVar", eVar="eventVar");
kmScan_ov2 <- data.frame(t(data.frame(kmScan_ov2)));
colnames(kmScan_ov2) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Ovarian");

kmScan_pr2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, pr2, F, tVar="TimeVar", eVar="eventVar");
kmScan_pr2 <- data.frame(t(data.frame(kmScan_pr2)));
colnames(kmScan_pr2) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Prostate");

kmScan_ki2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, ki2, F, tVar="TimeVar", eVar="eventVar");
kmScan_ki2 <- data.frame(t(data.frame(kmScan_ki2)));
colnames(kmScan_ki2) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Kidney");

kmScan_hn2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, hn2, F, tVar="TimeVar", eVar="eventVar");
kmScan_hn2 <- data.frame(t(data.frame(kmScan_hn2)));
colnames(kmScan_hn2) <- c("Gene", "P.Value", "Adj.P.Value");
print("Done Head and Neck");
print("Finished KM Scan 2");


stopCluster(clus);
print("Stopped cluster");

##########################################


##########################################
#Merge all results together into one data frame and a matrix for convenience DONE
print("Starting Data Merge");

results_ov1 <- cbind(kmScan_ov1[3], qCut50_ov1[2], qCut2575_ov1[2], coxReg_ov1[2], kmeans_ov1[2]);
colnames(results_ov1) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value", "Cox.P.Value", "kmeans.P.Value")
resultsMat_ov1 <- sapply(results_ov, FUN=f1actToNum);
rownames(resultsMat_ov1) <- rownames(resultsMat_ov1);
write.table(resultsMat_ov1, "allresults_ov1.txt", sep="\t", row.names=T);

results_pr1 <- cbind(kmScan_pr1[3], qCut50_pr1[2], qCut2575_pr1[2], coxReg_pr1[2], kmeans_pr1[2]);
colnames(results_pr1) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value","Cox.P.Value", "kmeans.P.Value")
resultsMat_pr1 <- sapply(results_pr1, FUN=factToNum);
rownames(resultsMat_pr1) <- rownames(resultsMat_pr1);
write.table(resultsMat_pr1, "allresults_pr1.txt", sep="\t", row.names=T);

results_ki1 <- cbind(kmScan_ki1[3], qCut50_ki1[2], qCut2575_ki1[2], coxReg_ki1[2], kmeans_ki1[2]);
colnames(results_ki1) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value", "Cox.P.Value", "kmeans.P.Value")
resultsMat_ki1 <- sapply(results_ki1, FUN=factToNum);
rownames(resultsMat_ki1) <- rownames(resultsMat_ki1);
write.table(resultsMat_ki1, "allresults_ki1.txt", sep="\t", row.names=T);

results_hn1 <- cbind(kmScan_hn1[3], qCut50_hn1[2], qCut2575_hn1[2], coxReg_hn1[2], kmeans_hn1[2]);
colnames(results_hn1) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value", "Cox.P.Value", "kmeans.P.Value")
resultsMat_hn1 <- sapply(results_hn1, FUN=factToNum);
rownames(resultsMat_hn1) <- rownames(resultsMat_hn1);
write.table(resultsMat_hn1, "allresults_hn1.txt", sep="\t", row.names=T);


results_ov2 <- cbind(kmScan_ov2[3], qCut50_ov2[2], qCut2575_ov2[2], coxReg_ov2[2], kmeans_ov2[2]);
colnames(results_ov2) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value", "Cox.P.Value", "kmeans.P.Value")
resultsMat_ov2 <- sapply(results_ov2, FUN=factToNum);
rownames(resultsMat_ov2) <- rownames(resultsMat_ov2);
write.table(resultsMat_ov2, "allresults_ov2.txt", sep="\t", row.names=T);

results_pr2 <- cbind(kmScan_pr2[3], qCut50_pr2[2], qCut2575_pr2[2], coxReg_pr2[2], kmeans_pr2[2]);
colnames(results_pr2) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value","Cox.P.Value", "kmeans.P.Value")
resultsMat_pr2 <- sapply(results_pr2, FUN=factToNum);
rownames(resultsMat_pr2) <- rownames(resultsMat_pr2);
write.table(resultsMat_pr2, "allresults_pr2.txt", sep="\t", row.names=T);

results_ki2 <- cbind(kmScan_ki2[3], qCut50_ki2[2], qCut2575_ki2[2], coxReg_ki2[2], kmeans_ki2[2]);
colnames(results_ki2) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value", "Cox.P.Value", "kmeans.P.Value")
resultsMat_ki2 <- sapply(results_ki2, FUN=factToNum);
rownames(resultsMat_ki2) <- rownames(resultsMat_ki2);
write.table(resultsMat_ki2, "allresults_ki2.txt", sep="\t", row.names=T);

results_hn2 <- cbind(kmScan_hn2[3], qCut50_hn2[2], qCut2575_hn2[2], coxReg_hn2[2], kmeans_hn2[2]);
colnames(results_hn2) <- c("Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value", "Cox.P.Value", "kmeans.P.Value")
resultsMat_hn2 <- sapply(results_hn2, FUN=factToNum);
rownames(resultsMat_hn2) <- rownames(resultsMat_hn2);
write.table(resultsMat_hn2, "allresults_hn2.txt", sep="\t", row.names=T);


















##########################################

Status API Training Shop Blog About Pricing
© 2015 GitHub, Inc. Terms Privacy Security Contact Help
