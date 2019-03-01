
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
png("heatmap_correlation_ov");
pheatmap(cor(-log10(resultsMat_ov)), color = colorRampPalette(c("green", "black", "red"))(50));
dev.off();
png("heatmap_correlation_pr");
pheatmap(cor(-log10(resultsMat_pr)), color = colorRampPalette(c("green", "black", "red"))(50));
dev.off();
png("heatmap_correlation_ki");
pheatmap(cor(-log10(resultsMat_ki)), color = colorRampPalette(c("green", "black", "red"))(50));
dev.off();
png("heatmap_correlation_ov");
pheatmap(cor(-log10(resultsMat_hn)), color = colorRampPalette(c("green", "black", "red"))(50));
dev.off();


#2. Correlations
png("Scatterplot_correlations_ov");
splom(-log10(resultsMat_ov));
dev.off();
png("Scatterplot_correlations_pr");
splom(-log10(resultsMat_pr));
dev.off();
png("Scatterplot_correlations_ki");
splom(-log10(resultsMat_ki));
dev.off();
png("Scatterplot_correlations_hn");
splom(-log10(resultsMat_hn));
dev.off();

#3. Venn Diagram
createInput <- function(data, pvalCutoff=.01)
{
    output <- list();
    for(i in 1:length(resultsMat_ov))
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






Status API Training Shop Blog About Pricing
© 2015 GitHub, Inc. Terms Privacy Security Contact Help
