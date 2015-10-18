
##########################################
#Objective : Compare various methods of dealing with
#continuous variables during a survival analysis. 
#
##########################################


#Call libraries
library(snow);
library(mixtools);
library(ggplot2);


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
colnames(exprs_ov) <- gsub("\\.", "-", colnames(exprs_ov))
rownames(exprs_ov) <- gsub(" ", "", rownames(exprs_ov));
ov <- list(exprs_ov, annot_ov);

#prostate
annot_pr <- read.delim("../data/raw/prca/annot.txt");
exprs_pr <- read.delim("../data/raw/prca/exprs.txt")
colnames(exprs_pr) <- gsub("\\.", "-", colnames(exprs_pr));
rownames(exprs_pr) <- gsub(" ", "", rownames(exprs_pr));
pr <- list(exprs_pr, annot_pr);

#head and neck
annot_hn <- read.delim("../data/raw/hnca/annot.txt");
exprs_hn <- read.delim("../data/raw/hnca/exprs.txt")
colnames(exprs_hn) <- gsub("\\.", "-", colnames(exprs_hn))
rownames(exprs_hn) <- gsub(" ", "", rownames(exprs_hn));
hn <- list(exprs_hn, annot_hn);

#Kidney
annot_ki <- read.delim("../data/raw/kica/annot.txt");
exprs_ki <- read.delim("../data/raw/kica/exprs.txt")
colnames(exprs_ki) <- gsub("\\.", "-", colnames(exprs_ki))
rownames(exprs_ki) <- gsub(" ", "", rownames(exprs_ki));
ki <- list(exprs_ki, annot_ki);

print("Finished Reading Data");
##########################################

##########################################
#Create two datasets from 1

#Ovarian
len1 <- dim(annot_ov[annot_ov[,"eventVar"]==1,])[1];
len0 <- dim(annot_ov[annot_ov[,"eventVar"]==0,])[1];
ov_an_1_p1 <- rownames(annot_ov[annot_ov[,"eventVar"]==1,])[1:round((len1/2))]
ov_an_1_p2 <- setdiff(rownames(annot_ov[annot_ov[,"eventVar"]==1,]), ov_an_1_p1)
ov_an_0_p1 <- rownames(annot_ov[annot_ov[,"eventVar"]==0,])[1:round((len0/2))]
ov_an_0_p2 <- setdiff(rownames(annot_ov[annot_ov[,"eventVar"]==0,]), ov_an_0_p1)
ov_an_1 <- annot_ov[c(ov_an_1_p1,ov_an_0_p1),]
ov_an_2 <- annot_ov[c(ov_an_1_p2,ov_an_0_p2),];
ov_exp_1 <- exprs_ov[,c(ov_an_1_p1,ov_an_0_p1)];
ov_exp_2 <- exprs_ov[,c(ov_an_1_p2,ov_an_0_p2)];
ov1 <- list(ov_exp_1, ov_an_1);
ov2 <- list(ov_exp_2, ov_an_2);

#Prostate
len1 <- dim(annot_pr[annot_pr[,"eventVar"]==1,])[1];
len0 <- dim(annot_pr[annot_pr[,"eventVar"]==0,])[1];
pr_an_1_p1 <- rownames(annot_pr[annot_pr[,"eventVar"]==1,])[1:round((len1/2))]
pr_an_1_p2 <- setdiff(rownames(annot_pr[annot_pr[,"eventVar"]==1,]), pr_an_1_p1)
pr_an_0_p1 <- rownames(annot_pr[annot_pr[,"eventVar"]==0,])[1:round((len0/2))]
pr_an_0_p2 <- setdiff(rownames(annot_pr[annot_pr[,"eventVar"]==0,]), pr_an_0_p1)
pr_an_1 <- annot_pr[c(pr_an_1_p1,pr_an_0_p1),]
pr_an_2 <- annot_pr[c(pr_an_1_p2,pr_an_0_p2),];
pr_exp_1 <- exprs_pr[,c(pr_an_1_p1,pr_an_0_p1)];
pr_exp_2 <- exprs_pr[,c(pr_an_1_p2,pr_an_0_p2)];
pr1 <- list(pr_exp_1, pr_an_1);
pr2 <- list(pr_exp_2, pr_an_2);

#Kidney
len1 <- dim(annot_ki[annot_ki[,"eventVar"]==1,])[1];
len0 <- dim(annot_ki[annot_ki[,"eventVar"]==0,])[1];
ki_an_1_p1 <- rownames(annot_ki[annot_ki[,"eventVar"]==1,])[1:round((len1/2))]
ki_an_1_p2 <- setdiff(rownames(annot_ki[annot_ki[,"eventVar"]==1,]), ki_an_1_p1)
ki_an_0_p1 <- rownames(annot_ki[annot_ki[,"eventVar"]==0,])[1:round((len0/2))]
ki_an_0_p2 <- setdiff(rownames(annot_ki[annot_ki[,"eventVar"]==0,]), ki_an_0_p1)
ki_an_1 <- annot_ki[c(ki_an_1_p1,ki_an_0_p1),]
ki_an_2 <- annot_ki[c(ki_an_1_p2,ki_an_0_p2),];
ki_exp_1 <- exprs_ki[,c(ki_an_1_p1,ki_an_0_p1)];
ki_exp_2 <- exprs_ki[,c(ki_an_1_p2,ki_an_0_p2)];
ki1 <- list(ki_exp_1, ki_an_1);
ki2 <- list(ki_exp_2, ki_an_2);

#Head and Neck
len1 <- dim(annot_hn[annot_hn[,"eventVar"]==1,])[1];
len0 <- dim(annot_hn[annot_hn[,"eventVar"]==0,])[1];
hn_an_1_p1 <- rownames(annot_hn[annot_hn[,"eventVar"]==1,])[1:round((len1/2))]
hn_an_1_p2 <- setdiff(rownames(annot_hn[annot_hn[,"eventVar"]==1,]), hn_an_1_p1)
hn_an_0_p1 <- rownames(annot_hn[annot_hn[,"eventVar"]==0,])[1:round((len0/2))]
hn_an_0_p2 <- setdiff(rownames(annot_hn[annot_hn[,"eventVar"]==0,]), hn_an_0_p1)
hn_an_1 <- annot_hn[c(hn_an_1_p1,hn_an_0_p1),]
hn_an_2 <- annot_hn[c(hn_an_1_p2,hn_an_0_p2),];
hn_exp_1 <- exprs_hn[,c(hn_an_1_p1,hn_an_0_p1)];
hn_exp_2 <- exprs_hn[,c(hn_an_1_p2,hn_an_0_p2)];
hn1 <- list(hn_exp_1, hn_an_1);
hn2 <- list(hn_exp_2, hn_an_2);






##########################################






##########################################
#Choose genes with most variability DONE
geneCV <- apply(exprs_ov, FUN=max, MARGIN=1);
geneCV <- sort(geneCV, T);
geneCV <- gsub(" ", "", geneCV);

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

clusterExport(clus, "ov1");
clusterExport(clus, "pr1");
clusterExport(clus, "ki1");
clusterExport(clus, "hn1");
clusterExport(clus, "ov2");
clusterExport(clus, "pr2");
clusterExport(clus, "ki2");
clusterExport(clus, "hn2");
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

#Kmneans
kmeans_ov1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, ov1, tVar="TimeVar", eVar="eventVar");
kmeans_ov1 <- data.frame(t(data.frame(kmeans_ov1)));
colnames(kmeans_ov1) <- c("Gene", "P.Value");
write.table(kmeans_ov1, "kmeans_ov1.txt", sep="\t", row.names=F);
print("Done Ovarian 1");

kmeans_pr1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, pr1, tVar="TimeVar", eVar="eventVar");
kmeans_pr1 <- data.frame(t(data.frame(kmeans_pr1)));
colnames(kmeans_pr1) <- c("Gene", "P.Value");
write.table(kmeans_pr1, "kmeans_pr1.txt", sep="\t", row.names=F);
print("Done Prostate 1");

kmeans_ki1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, ki1, tVar="TimeVar", eVar="eventVar");
kmeans_ki1 <- data.frame(t(data.frame(kmeans_ki1)));
colnames(kmeans_ki1) <- c("Gene", "P.Value");
write.table(kmeans_ki1, "kmeans_ki1.txt", sep="\t", row.names=F);
print("Done Kidney 1");

kmeans_hn1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, hn1, tVar="TimeVar", eVar="eventVar");
kmeans_hn1 <- data.frame(t(data.frame(kmeans_hn1)));
colnames(kmeans_hn1) <- c("Gene", "P.Value");
write.table(kmeans_hn1, "kmeans_hn1.txt", sep="\t", row.names=F);
print("Done Head Neck 1");
print("Finished KMeans 1");

#Kmneans
kmeans_ov2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, ov2, tVar="TimeVar", eVar="eventVar");
kmeans_ov2 <- data.frame(t(data.frame(kmeans_ov2)));
colnames(kmeans_ov2) <- c("Gene", "P.Value");
write.table(kmeans_ov2, "kmeans_ov2.txt", sep="\t", row.names=F);
print("Done Ovarian 2");

kmeans_pr2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, pr2, tVar="TimeVar", eVar="eventVar");
kmeans_pr2 <- data.frame(t(data.frame(kmeans_pr2)));
colnames(kmeans_pr2) <- c("Gene", "P.Value");
write.table(kmeans_pr2, "kmeans_pr2.txt", sep="\t", row.names=F);
print("Done Prostate 2");

kmeans_ki2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, ki2, tVar="TimeVar", eVar="eventVar");
kmeans_ki2 <- data.frame(t(data.frame(kmeans_ki2)));
colnames(kmeans_ki2) <- c("Gene", "P.Value");
write.table(kmeans_ki2, "kmeans_ki2.txt", sep="\t", row.names=F);
print("Done Kidney 2");

kmeans_hn2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= kmeansSA, hn2, tVar="TimeVar", eVar="eventVar");
kmeans_hn2 <- data.frame(t(data.frame(kmeans_hn2)));
colnames(kmeans_hn2) <- c("Gene", "P.Value");
write.table(kmeans_hn2, "kmeans_hn2.txt", sep="\t", row.names=F);
print("Done Head Neck 2");
print("Finished KMeans 2");



#Cox Regression
coxReg_pr1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, pr1);
coxReg_pr1 <- data.frame(t(data.frame(coxReg_pr1)));
colnames(coxReg_pr1) <- c("Gene", "P.Value");
write.table(coxReg_pr1, "coxreg_pr1.txt", sep="\t", row.names=F);
print("Done Prostate 1");

coxReg_ov1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, ov1);
coxReg_ov1 <- data.frame(t(data.frame(coxReg_ov1)));
colnames(coxReg_ov1) <- c("Gene", "P.Value");
write.table(coxReg_ov1, "coxreg_ov1.txt", sep="\t", row.names=F);
print("Done Ovarian 1");

coxReg_ki1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, ki1);
coxReg_ki1 <- data.frame(t(data.frame(coxReg_ki1)));
colnames(coxReg_ki1) <- c("Gene", "P.Value");
write.table(coxReg_ki1, "coxreg_ki1.txt", sep="\t", row.names=F);
print("Done Kidney 1");

coxReg_hn1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, hn1);
coxReg_hn1 <- data.frame(t(data.frame(coxReg_hn1)));
colnames(coxReg_hn1) <- c("Gene", "P.Value");
write.table(coxReg_hn1, "coxreg_hn1.txt", sep="\t", row.names=F);
print("Done Head and Neck 1");
print("Finished cox regression 1");

coxReg_pr2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, pr2);
coxReg_pr2 <- data.frame(t(data.frame(coxReg_pr2)));
colnames(coxReg_pr2) <- c("Gene", "P.Value");
write.table(coxReg_pr2, "coxreg_pr2.txt", sep="\t", row.names=F);
print("Done Prostate 2");

coxReg_ov2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, ov2);
coxReg_ov2 <- data.frame(t(data.frame(coxReg_ov2)));
colnames(coxReg_ov2) <- c("Gene", "P.Value");
write.table(coxReg_ov2, "coxreg_ov2.txt", sep="\t", row.names=F);
print("Done Ovarian 2");

coxReg_ki2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, ki2);
coxReg_ki2 <- data.frame(t(data.frame(coxReg_ki2)));
colnames(coxReg_ki2) <- c("Gene", "P.Value");
write.table(coxReg_ki2, "coxreg_ki2.txt", sep="\t", row.names=F);
print("Done Kidney 2");

coxReg_hn2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN= coxReg, hn2);
coxReg_hn2 <- data.frame(t(data.frame(coxReg_hn2)));
colnames(coxReg_hn2) <- c("Gene", "P.Value");
write.table(coxReg_hn2, "coxreg_hn2.txt", sep="\t", row.names=F);
print("Done Head and Neck 2");
print("Finished cox regression 2");


#Quantile Technique, cutting at median
qCut50_ov1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ov1, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_ov1 <- data.frame(t(data.frame(qCut50_ov1)));
colnames(qCut50_ov1) <- c("Gene", "P.Value");
write.table(qCut50_ov1, "qCut50_ov1.txt", sep="\t", row.names=F);
print("Done Ovarian 1");

qCut50_pr1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, pr1, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_pr1 <- data.frame(t(data.frame(qCut50_pr1)));
colnames(qCut50_pr1) <- c("Gene", "P.Value");
write.table(qCut50_pr1, "qCut50_pr1.txt", sep="\t", row.names=F);
print("Done Prostate 1");

qCut50_ki1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ki1, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_ki1 <- data.frame(t(data.frame(qCut50_ki1)));
colnames(qCut50_ki1) <- c("Gene", "P.Value");
write.table(qCut50_ki1, "qCut50_ki1.txt", sep="\t", row.names=F);
print("Done Kidney 1");

qCut50_hn1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, hn1, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_hn1 <- data.frame(t(data.frame(qCut50_hn1)));
colnames(qCut50_hn1) <- c("Gene", "P.Value");
write.table(qCut50_hn1, "qCut50_hn1.txt", sep="\t", row.names=F);
print("Done Head and Neck 1");
print("Finished median quantile cut 1");


qCut50_ov2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ov2, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_ov2 <- data.frame(t(data.frame(qCut50_ov2)));
colnames(qCut50_ov2) <- c("Gene", "P.Value");
write.table(qCut50_ov2, "qCut50_ov2.txt", sep="\t", row.names=F);
print("Done Ovarian 2");

qCut50_pr2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, pr2, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_pr2 <- data.frame(t(data.frame(qCut50_pr2)));
colnames(qCut50_pr2) <- c("Gene", "P.Value");
write.table(qCut50_pr2, "qCut50_pr2.txt", sep="\t", row.names=F);
print("Done Prostate 2");

qCut50_ki2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ki2, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_ki2 <- data.frame(t(data.frame(qCut50_ki2)));
colnames(qCut50_ki2) <- c("Gene", "P.Value");
write.table(qCut50_ki2, "qCut50_ki2.txt", sep="\t", row.names=F);
print("Done Kidney 2");

qCut50_hn2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, hn2, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_hn2 <- data.frame(t(data.frame(qCut50_hn2)));
colnames(qCut50_hn2) <- c("Gene", "P.Value");
write.table(qCut50_hn2, "qCut50_hn2.txt", sep="\t", row.names=F);
print("Done Head and Neck 2");
print("Finished median quantile cut 2");



#Quantile Technique, cutting at 25 and 75
qCut2575_ov1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ov1, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_ov1 <- data.frame(t(data.frame(qCut2575_ov1)));
colnames(qCut2575_ov1) <- c("Gene", "P.Value");
write.table(qCut2575_ov1, "qCut2575_ov1 .txt", sep="\t", row.names=F);
print("Done Ovarian 1");

qCut2575_pr1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, pr1, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_pr1 <- data.frame(t(data.frame(qCut2575_pr1)));
colnames(qCut2575_pr1) <- c("Gene", "P.Value");
write.table(qCut2575_pr1, "qCut2575_pr1.txt", sep="\t", row.names=F);
print("Done Prostate 1");

qCut2575_ki1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ki1, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_ki1 <- data.frame(t(data.frame(qCut2575_ki1)));
colnames(qCut2575_ki1) <- c("Gene", "P.Value");
write.table(qCut2575_ki1, "qCut2575_ki1.txt", sep="\t", row.names=F);
print("Done Kidney 1");

qCut2575_hn1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, hn1, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_hn1 <- data.frame(t(data.frame(qCut2575_hn1)));
colnames(qCut2575_hn1) <- c("Gene", "P.Value");
write.table(qCut2575_hn1, "qCut2575_hn1.txt", sep="\t", row.names=F);
print("Done Head and Neck 1");
print("Finished 75th 25th quantile cut 1");

qCut2575_ov2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ov2, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_ov2 <- data.frame(t(data.frame(qCut2575_ov2)));
colnames(qCut2575_ov2) <- c("Gene", "P.Value");
write.table(qCut2575_ov2, "qCut2575_ov2.txt", sep="\t", row.names=F);
print("Done Ovarian 2");

qCut2575_pr2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, pr2, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_pr2 <- data.frame(t(data.frame(qCut2575_pr2)));
colnames(qCut2575_pr2) <- c("Gene", "P.Value");
write.table(qCut2575_pr2, "qCut2575_pr2.txt", sep="\t", row.names=F);
print("Done Prostate 2");

qCut2575_ki2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, ki2, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_ki2 <- data.frame(t(data.frame(qCut2575_ki2)));
colnames(qCut2575_ki2) <- c("Gene", "P.Value");
write.table(qCut2575_ki2, "qCut2575_ki2.txt", sep="\t", row.names=F);
print("Done Kidney 2");

qCut2575_hn2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=quantCutSA, hn2, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_hn2 <- data.frame(t(data.frame(qCut2575_hn2)));
colnames(qCut2575_hn2) <- c("Gene", "P.Value");
write.table(qCut2575_hn2, "qCut2575_hn.txt", sep="\t", row.names=F);
print("Done Head and Neck 2");
print("Finished 75th 25th quantile cut 2");



#KM Scan Technique
kmScan_ov1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, ov1, F, tVar="TimeVar", eVar="eventVar");
kmScan_ov1 <- data.frame(t(data.frame(kmScan_ov1)));
colnames(kmScan_ov1) <- c("Gene", "P.Value", "Adj.P.Value");
write.table(kmScan_ov1, "kmScan_ov1", sep="\t", row.names=F);
print("Done Ovarian 1");

kmScan_pr1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, pr1, F, tVar="TimeVar", eVar="eventVar");
kmScan_pr1 <- data.frame(t(data.frame(kmScan_pr1)));
colnames(kmScan_pr1) <- c("Gene", "P.Value", "Adj.P.Value");
write.table(kmScan_pr1, "kmScan_pr1", sep="\t", row.names=F);
print("Done Prostate 1");

kmScan_ki1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, ki1, F, tVar="TimeVar", eVar="eventVar");
kmScan_ki1 <- data.frame(t(data.frame(kmScan_ki1)));
colnames(kmScan_ki1) <- c("Gene", "P.Value", "Adj.P.Value");
write.table(kmScan_ki1, "kmScan_ki1", sep="\t", row.names=F);
print("Done Kidney 1");

kmScan_hn1 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, hn1, F, tVar="TimeVar", eVar="eventVar");
kmScan_hn1 <- data.frame(t(data.frame(kmScan_hn1)));
colnames(kmScan_hn1) <- c("Gene", "P.Value", "Adj.P.Value");
write.table(kmScan_hn1, "kmScan_hn1", sep="\t", row.names=F);
print("Done Head and Neck");
print("Finished KM Scan 1");


#KM Scan Technique
kmScan_ov2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, ov2, F, tVar="TimeVar", eVar="eventVar");
kmScan_ov2 <- data.frame(t(data.frame(kmScan_ov2)));
colnames(kmScan_ov2) <- c("Gene", "P.Value", "Adj.P.Value");
write.table(kmScan_ov2, "kmScan_ov2", sep="\t", row.names=F);
print("Done Ovarian 2");

kmScan_pr2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, pr2, F, tVar="TimeVar", eVar="eventVar");
kmScan_pr2 <- data.frame(t(data.frame(kmScan_pr2)));
colnames(kmScan_pr2) <- c("Gene", "P.Value", "Adj.P.Value");
write.table(kmScan_pr2, "kmScan_pr2", sep="\t", row.names=F);
print("Done Prostate 2");

kmScan_ki2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, ki2, F, tVar="TimeVar", eVar="eventVar");
kmScan_ki2 <- data.frame(t(data.frame(kmScan_ki2)));
colnames(kmScan_ki2) <- c("Gene", "P.Value", "Adj.P.Value");
write.table(kmScan_ki2, "kmScan_ki2", sep="\t", row.names=F);
print("Done Kidney 2");

kmScan_hn2 <- parSapply(clus, names(geneCV)[1:numGenes], FUN=kapmPlot, hn2, F, tVar="TimeVar", eVar="eventVar");
kmScan_hn2 <- data.frame(t(data.frame(kmScan_hn2)));
colnames(kmScan_hn2) <- c("Gene", "P.Value", "Adj.P.Value");
write.table(kmScan_hn2, "kmScan_hn2", sep="\t", row.names=F);
print("Done Head and Neck 2");
print("Finished KM Scan 2");


stopCluster(clus);
print("Stopped cluster");

##########################################
