


coxReg <- function(genes, myData)
{
    #Get metadata
    print(genes);
    tmpMeta <- myData[[2]];
    
    myGene <- myData[[1]][genes,];
    tmpMeta[,"Gene"] <- as.numeric(myGene);
    tmpMeta <- tmpMeta[order(tmpMeta[,"Gene"]),]
    
    coxExpAnalysis <- coxph(formula = Surv(TimeVar, eventVar) ~ Gene, data = tmpMeta)
    pVal <- summary(coxExpAnalysis)[7][[1]][5];
    out <- c(genes, pVal);
    out;
    
}
