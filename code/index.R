library(survcomp)######################################################
# Author: Komal Rathi
# Function: function to calculate c-index and d-index
# Date: 01/29/2018
######################################################

# formatted like other functions
# c-index
c.index <- function(genes, myData){
  tmpMeta <- myData[[2]]
  myGene <- myData[[1]][genes,]
  tmpMeta[,"Gene"] <- as.numeric(myGene)
  temp <- table(tmpMeta[,"Gene"])
  out <- c(genes, 1)
  if(names(temp)[temp == max(temp)]!=0){
    tt <- concordance.index(x = tmpMeta$Gene, surv.time = tmpMeta$TimeVar, surv.event = tmpMeta$eventVar, method = "noether", na.rm = TRUE)
    cindex <- c("cindex" = tt$c.index, "cindex.se" = tt$se, "lower" = tt$lower, "upper" = tt$upper, "p-value" = tt$p.value)
    out <- c(genes, cindex[['p-value']])
  }
  out
}

# d-index
d.index <- function(genes, myData){
  tmpMeta <- myData[[2]]
  myGene <- myData[[1]][genes,]
  tmpMeta[,"Gene"] <- as.numeric(myGene)
  temp <- table(tmpMeta[,"Gene"])
  out <- c(genes, 1)
  if(names(temp)[temp == max(temp)]!=0){
    tt <- D.index(x = tmpMeta$Gene, surv.time = tmpMeta$TimeVar, surv.event = tmpMeta$eventVar, na.rm = TRUE)
    dindex <- c("dindex" = tt$d.index, "dindex.se" = tt$se, "lower" = tt$lower, "upper" = tt$upper, "p-value" = tt$p.value)
    out <- c(genes, dindex[['p-value']])
  }
  out
}

