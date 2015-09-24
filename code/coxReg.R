


coxReg <- function(genes, myData)
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
    tmpMeta[,"Gene"] <- log2(tmpMeta[,"Gene"]+1);
    coxExpAnalysis <- coxph(formula = Surv(TimeVar, eventVar) ~ Gene, data = tmpMeta)
    pVal <- summary(coxExpAnalysis)[7][[1]][5];
    out <- c(genes, pVal);
    }
    out;
    
}
