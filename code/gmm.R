
#This will yield a p-value and show a plot
gmmSA <- function(genes, myData, createPlot=T, tVar="time", eVar="event")
{
    #Get metadata
    
    tmpMeta <- myData[[2]];
    myGene <- myData[[1]][genes,];
    tmpMeta[,"Gene"] <- as.numeric(myGene);
    out <- c(genes, 1);
    if(median(tmpMeta[,"Gene"])!=0)
    {
    tmpMix <- normalmixEM(tmpMeta[,"Gene"], maxit=1000, maxrestarts=20);
    myOut <- data.frame(tmpMix$posterior);
    tmpMeta[,"Groups"] <- myOut[,1]> myOut[,2];
    tmpMeta <- tmpMeta[order(tmpMeta[,"Gene"]),]
  
    
    

    #time
    timeVar <- tmpMeta[,tVar];
    #event
    eventVar <- tmpMeta[,eVar];
    
    #createsurvival
    t.Surv <- Surv(timeVar, eventVar);
    t.survfit <- survfit(t.Surv~Groups, data= tmpMeta);
    myP <- pchisq(survdiff(t.Surv~ Groups, data= tmpMeta, rho=0)$chisq, df=1, lower=F);
    out <- c(genes, myP)   
    }
    return(out);
}
