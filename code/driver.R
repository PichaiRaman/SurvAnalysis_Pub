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


##########################################
#Source all code here
#Methods of survival analysis 
source("KaplanScan.R");
source("quantileCut.R");
source("zCut.R");

#Some functions to help out 
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

resultsMat_TS <- melt(resultsMat);

write.table(results, "allresults.txt", sep="\t", row.names=T);
##########################################


##########################################
#Compare results
#
#1. Compare how close methods are via correlations, venn diagrams, heatmaps
#2. ROC Curves with some positive controls to determine best methods
#3. Accuracy of methods via a venn diagram and fister test
##########################################


#############################################
#2.Comparison of methods
#############################################


#1. Heatmap
png("heatmap_correlation");
pheatmap(cor(-log10(resultsMat)), color = colorRampPalette(c("green", "black", "red"))(50));
dev.off();

#2. Correlations
png("Scatterplot_correlations");
splom(-log10(resultsMat));
dev.off();

#3. Venn Diagram
createInput <- function(data, pvalCutoff=.01)
{
    output <- list();
    for(i in 1:length(resultsMat))
    {
        tmp <- rownames(data[data[,i]<pvalCutoff,]);
        output[[paste(colnames(data))[i]]] <- tmp;
    }
    return(output);
}

png("venn_correlation");
venn(createInput(resultsMat));
dev.off();


#############################################
#2. Now we need ROC curves for a list of genes
#############################################

#Plot a single ROC
plotROC <- function(data, myCol, cancerSet)
{
    notInCancerSet <- setdiff(rownames(data), cancerSet);
    myX <- data[notInCancerSet,myCol];
    myY <- data[cancerSet,myCol];
    myROC <- roc.curve( myX, myY, curve=TRUE );
    plot( myROC, color = "red", auc.main=T);
}

plotROC(data, colnames(data)[5], cancerSet);


#Plot all ROC's for the matrix
plotROCAll <- function(data, cancerSet)
{
    notInCancerSet <- setdiff(rownames(data), cancerSet);
    myX <- data[notInCancerSet,1];
    myY <- data[cancerSet,1];
    myROC <- roc.curve( myX, myY, curve=TRUE );
    print(colnames(data)[i]);
    print(myROC);
    plot( myROC, color = "red", auc.main=F);
    for(i in 2:length(data))
    {
        myCol <- i;
        myX <- data[notInCancerSet,myCol];
        myY <- data[cancerSet,myCol];
        myROC <- roc.curve( myX, myY, curve=TRUE );
        print(colnames(data)[i]);
        print(myROC);
        plot( myROC, color = i, auc.main=F, add=T);
    }
    
    
}

#plotROC(data, colnames(data)[5], cancerSet);

#plotROCAll(data, cancerSet);


#############################################
#3. Accuracy Method by Venn
#############################################

#Run Hypergeometric
intHypGeo <- function(uniSize,catSize,listSize,numMatched) {
    
    pv =  0.5*dhyper(numMatched,listSize,uniSize-listSize,catSize) + sum(dhyper((numMatched+1):catSize,listSize,uniSize-listSize,catSize));
    pv;
}  


accuVenn <- function(resultA, resultB, myCol, pvalCutoff=.01)
{
    g1 <- rownames(resultA[resultA[,myCol]<pvalCutoff,]);
    g2 <- rownames(resultB[resultB[,myCol]<pvalCutoff,]);
    myList <- list("Group1" = g1, "Group2"=g2)
    venn(myList, main="hi");
    intHypGeo(dim(resultA)[1], length(g1), length(g2), length(intersect(g1,g2)));
    
}




















