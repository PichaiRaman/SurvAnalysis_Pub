###################################
# Author: Pichai Raman
# Function: Code to respond to reviewer & show cut-offs are sensitive to composition of datasets
# Date: 1/31/2018
###################################

# Call libaries
library(tidyverse)
library(car)

# Read datasets
load('../data/ALLDATA.RData')

generateDist <- function(myGene, myDat, myDatName) {
  meanVec <- c()
  q75Vec <- c()
  q25Vec <- c()
  
  for(i in 1:500) {
    tmp <- log2(sample(as.numeric(myDat[myGene,]), 50)+1);
    mean <- median(tmp);
    meanVec <- c(meanVec, mean);
    q75 <- quantile(tmp)[[4]]
    q75Vec <- c(q75Vec, q75);
    q25 <- quantile(tmp)[[2]]
    q25Vec <- c(q25Vec, q25);
  }
  meanVec <- meanVec-mean(meanVec)
  q75Vec <- q75Vec-mean(q75Vec)
  q25Vec <- q25Vec-mean(q25Vec)
  output <- data.frame(meanVec, q75Vec, q25Vec, myGene, myDatName);
  colnames(output)[1:3] <- c("Median", "Q3", "Q1")
  return(output);
}

generateAll <- function(myGene, returnPlot=TRUE) {
  geneVectHN <- generateDist(myGene, exprs_hn, "HN")
  geneVectOV <- generateDist(myGene, exprs_ov, "OV")
  geneVectPR <- generateDist(myGene, exprs_pr, "PR")
  geneVectKI <- generateDist(myGene, exprs_ki, "KI")
  geneVect <- rbind(geneVectHN, geneVectOV, geneVectPR, geneVectKI)
  ltp_mean <- with(geneVect, leveneTest(Median, myDatName))[[3]][1]
  ltp_q25 <- with(geneVect, leveneTest(Q1, myDatName))[[3]][1]
  ltp_q75 <- with(geneVect, leveneTest(Q3, myDatName))[[3]][1]
  pOut <- c(ltp_mean, ltp_q25, ltp_q75)	
  geneVectTS <- gather(geneVect, key="MeasureType", value="Value", -myDatName, -myGene)
  geneVectTS$MeasureType[geneVectTS$MeasureType == "Median"] <- paste0("Median (P-Value: ", format(ltp_mean, scientific = T, digits = 3), ")")
  geneVectTS$MeasureType[geneVectTS$MeasureType == "Q1"] <- paste0("Q1 (P-Value: ", format(ltp_q25, scientific = T, digits = 3), ")")
  geneVectTS$MeasureType[geneVectTS$MeasureType == "Q3"] <- paste0("Q3 (P-Value: ", format(ltp_q75, scientific = T, digits = 3), ")")
  p <- ggplot(geneVectTS, aes(Value)) + geom_histogram(bins=100) + 
    geom_density(alpha=.2) + theme_bw() + 
    theme(axis.text = element_text(size = 14, colour = "black"),
          axis.title = element_text(size = 14, colour = "black"),
          strip.text = element_text(size = 12, colour = "black"))
  p <- p + geom_vline(xintercept=0) + facet_grid(myDatName ~ MeasureType) + xlab("Mean-centered Statistic (Log2 FPKM)")
  if(returnPlot)
  {
    return(list(p, pOut))
  }
  if(returnPlot)
  {
    return(pOut)
  }
}

POLA1 <- generateAll("POLA1|5422")
ggsave("../Figures/SampleMeanMedianQuantile.png", width=12, height=10, plot = POLA1[[1]])

