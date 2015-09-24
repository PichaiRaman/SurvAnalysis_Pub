
#This will yield a p-value and show a plot
kmeansSA <- function(genes, myData, createPlot=T, tVar="time", eVar="event")
{
    #Get metadata
    print(genes);
    tmpMeta <- myData[[2]];
    myGene <- myData[[1]][genes,];
    tmpMeta[,"Gene"] <- as.numeric(myGene);
    temp <- table(tmpMeta[,"Gene"])
    out <- c(genes, 1);
    if(names(temp)[temp == max(temp)]!=0)
    {
    set.seed(314);
    tmpMeta[,"Groups"] <- kmeans(tmpMeta[,"Gene"], centers=2)[[1]]-1
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
