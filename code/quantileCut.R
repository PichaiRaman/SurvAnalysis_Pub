########################################
#Use a quantile-based cut point to dichotimize continuous
#Variable prior to survival analysis
########################################


#Call libraries
library("survival");


#This will yield a p-value and show a plot
quantCutSA <- function(genes, myData, createPlot=T, quantLow=.25,  quantHigh=.75, tVar="time", eVar="event")
{
    #Get metadata
    
    tmpMeta <- myData[[2]];
    
    #Now pick a gene/genes &
    if(length(genes)==1)
    {
        myGene <- myData[[1]][genes,];
        tmpMeta[,"Gene"] <- as.numeric(myGene);
        tmpMeta <- tmpMeta[order(tmpMeta[,"Gene"]),]
    }
    
    if(length(genes)>1)
    {
        myGene <- myData[[1]][genes,];
        myGene <- apply(myGene, FUN=standScore, MARGIN=1);
        myGene <- apply(myGene, FUN=max, MARGIN=1);
        tmpMeta[,"Gene"] <- as.numeric(myGene);
        tmpMeta <- tmpMeta[order(tmpMeta[,"Gene"]),]
        genes <- paste(genes, collapse=" / ");
    }
    
    #Run quantile
    tmpMetaScan <- tmpMeta;
    cutOffs <- quantile(tmpMeta[,"Gene"], c(quantLow, quantHigh));
    hCut <- cutOffs[[2]]
    lCut <-  cutOffs[[1]]
    tmpMetaScan[,"GeneBin"] <- -2;
    tmpMetaScan[tmpMetaScan[, "Gene"]>hCut,"GeneBin"] <- 1;
    tmpMetaScan[tmpMetaScan[, "Gene"]<=lCut,"GeneBin"] <- 0;
    tmpMetaScan <- tmpMetaScan[ tmpMetaScan[,"GeneBin"]>(-1),]
    myReturn <- c(genes, 1);
    if(length(unique(as.character(tmpMetaScan[,"GeneBin"])))==2)
    {
    #time
    timeVar <- tmpMetaScan[,tVar];
    #event
    eventVar <- tmpMetaScan[,eVar];
    
    #createsurvival
    t.Surv <- Surv(timeVar, eventVar);
    t.survfit <- survfit(t.Surv~GeneBin, data=tmpMetaScan);
    myP <- pchisq(survdiff(t.Surv~GeneBin, data=tmpMetaScan, rho=0)$chisq, df=1, lower=F);
    out <- c(genes, myP)
    
    
    #Change strata names
    myLowName <- paste("Low : n = ", t.survfit$n[[1]], sep="");
    myHighName <- paste("High : n = ", t.survfit$n[[2]], sep="");
    names(t.survfit$strata) <- c(myLowName, myHighName)
    t.survframe <- createSurvivalFrame(t.survfit)
    
    if(createPlot==T)
    {
        tmpTitle <- paste("KM Plot -", genes, "\nP-val(Adj) :", format(out[2], scientific=T, digits=3));
        myReturn <- qplot_survival(t.survframe, f.CI=F, myTitle=tmpTitle)+theme_bw()+scale_colour_manual(values=c("red", "blue") );
    }
    
    if(createPlot==F)
    {
        
        myReturn <- c(genes, out[2]);
        
    }
    }
    myReturn;
    
}
