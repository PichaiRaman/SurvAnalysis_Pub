##########################################
#Objective : Compare various methods of dealing with
#continuous variables during a survival analysis. 
#
##########################################


#Call libraries
library("GGally");
library("ggplot2");
library("lattice");


##########################################
#Source all code here
#Methods of survival analysis 
source("KaplanScan.R");
source("quantileCut.R");
source("zCut.R");
source("helper.R");

##########################################

##########################################
#Read in data
annot_nbv <- read.delim("../data/Neuroblastoma_Versteeg/annot.txt")
exprs_nbv <- read.delim("../data/Neuroblastoma_Versteeg/exprs.txt")
nbv <- list(exprs_nbv, annot_nbv);
##########################################

##########################################
#Choose genes with most variability
geneCV <- apply(exprs_nbv, FUN=myCV, MARGIN=1);
geneCV <- sort(geneCV, T);

##########################################

##########################################
#Run Surival analysis using various techniques
#
#numGenes here is the number of genes we want to run survival analysis on, we should probably
#do all genes but kmScan takes a bit of time so for testing purposes let's set n to a small number
#
numGenes <- 100;

#KM Scan Technique
kmScan_nbv <- sapply(names(geneCV)[1:numGenes], FUN=kapmPlot, nbv, F, tVar="nti_surv_overall", eVar="nti_event_overall_num");
kmScan_nbv <- data.frame(t(data.frame(kmScan_nbv)));
colnames(kmScan_nbv) <- c("Gene", "P.Value", "Adj.P.Value");

#Quantile Technique, cutting at median
qCut50_nbv <- sapply(names(geneCV)[1:numGenes], FUN=quantCutSA, nbv, F, quantLow=.50,  quantHigh=.50, tVar="nti_surv_overall", eVar="nti_event_overall_num");
qCut50_nbv <- data.frame(t(data.frame(qCut50_nbv)));
colnames(qCut50_nbv) <- c("Gene", "P.Value");

#Quantile Technique, cutting at 25 and 75
qCut2575_nbv <- sapply(names(geneCV)[1:numGenes], FUN=quantCutSA, nbv, F, quantLow=.25,  quantHigh=.75, tVar="nti_surv_overall", eVar="nti_event_overall_num");
qCut2575_nbv <- data.frame(t(data.frame(qCut2575_nbv)));
colnames(qCut2575_nbv) <- c("Gene", "P.Value");

#Z Threshold Technique, cutting at 1 and -1
zCut1_nbv <- sapply(names(geneCV)[1:numGenes], FUN=zCutSA, nbv, F, zLow=-1,  zHigh=1, tVar="nti_surv_overall", eVar="nti_event_overall_num");
zCut1_nbv <- data.frame(t(data.frame(zCut1_nbv)));
colnames(zCut1_nbv) <- c("Gene", "P.Value");

##########################################


##########################################
#Merge all results together into one data frame and a matrix for convenience

results <- cbind(kmScan_nbv[2:3], qCut50_nbv[2], qCut2575_nbv[2], zCut1_nbv[2]);
colnames(results) <- c("Km.Scan.P.Value", "Km.Scan.Adj.P.Value", "Qcut50.P.Value", "Qcut2575.P.Value", "Zcut1.P.Value")
resultsMat <- sapply(results, FUN=factToNum);
rownames(resultsMat) <- rownames(results);
write.table(results, "allresults.txt", sep="\t", row.names=T);
##########################################


##########################################
#Compare results
#
#1. Compare how close methods are via correlations, venn diagrams, heatmaps
#2. ROC Curves with some positive controls to determine best methods
#3. Probably some way of looking at how robust each method is to different distributions, outliers, etc..
#
#
#


#Similarity of results
png("heatmap correlation");
pheatmap(cor(-log10(resultsMat)));
dev.off();

png("Scatterplot correlations");
splom(-log10(resultsMat));
dev.off();







##########################################













