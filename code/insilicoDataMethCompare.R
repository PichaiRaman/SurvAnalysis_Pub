##############################################
#Code to simulate RNA-Seq data and see
#how well methods can deal with noise 
#in the system
#
##############################################

#Call libraries
library("SimSeq");
source("~/rocon.R");


##########################################
#Source all code here
#Methods of survival analysis 
source("../code/KaplanScan.R");
source("../code/quantileCut.R");
source("../code/coxReg.R");
source("../code/kmeansSA.R");

#Some functions to help out 
source("../code/helper.R");

##########################################


#########################################
#Read in data DONE
load("../data/ALLDATA.RData");

#########################################
#Generate simulated data
set.seed(100);
simDat <- SimData(counts = round(as.matrix(exprs_hn)), treatment=round(rnorm(520, mean=.5, sd=.01)), n.genes=5000, n.diff=0, k.ind=85, sort.method="unpaired");

#Counts
simExprs <- data.frame(simDat$counts);
simExprs <- simExprs[,1:100];

#metaData
simMeta <- data.frame(colnames(simExprs), simDat$treatment[1:100]);
set.seed(100);

survTimes1 <- sort(annot_ki[annot_ki["eventVar"]==0, "TimeVar"])[(nrow(annot_ki[annot_ki["eventVar"]==0,])-84): nrow(annot_ki[annot_ki["eventVar"]==0,])];
survTimes2 <- sort(annot_ki[annot_ki["eventVar"]==1, "TimeVar"])[c(1:15)];
survTimes <- c(survTimes1, survTimes2);

simMeta[,"TimeVar"] <- survTimes;
colnames(simMeta)[1:2] <- c("Sample", "eventVar");

#Now Let's generate positive controls
genPosControlMat <- function(x)
{
tmpMean <- mean(x)
tmpSrv <- survTimes[86:100];
survTimesNorm <- c(rnorm(85, mean=.25 sd=.1),rnorm(15, mean=1 sd=.1))
multiplier <- abs(rnorm(1, mean=1, sd=1));
xN <- round(multiplier*x*survTimesNorm);
}
#Get genes with only values > 100
posControlGenes <- rownames(simExprs)%in%sample(rownames(simExprs)[apply(simExprs, FUN=mean, MARGIN=1)>100], 100)
simExprs[posControlGenes,] <- simExprs[posControlGenes,]+t(apply(simExprs[posControlGenes,], FUN=genPosControlMat, MARGIN=1));
simObjNoNoise <- list(simExprs, simMeta);

addNoise <- function(myPerc)
{
myMeans <-as.numeric(apply(simExprs, FUN=mean, MARGIN=2))*myPerc;
mySD <-as.numeric(apply(simExprs, FUN=sd, MARGIN=2))*myPerc;

myDat <- rnorm(nrow(simExprs), mean=myMeans[1], sd=mySD[1])
for(i in 2:length(myMeans))
{
tmpMat <- rnorm(nrow(simExprs), mean=myMeans[i], sd=mySD[i]);
myDat <- cbind(myDat, tmpMat);
}
return(as.matrix(myDat));
}


CreateMatrix <- function(simObj)
{
kmeans_sim <- sapply(rownames(simExprs), FUN= kmeansSA, simObj, tVar="TimeVar", eVar="eventVar");
kmeans_sim <- data.frame(t(data.frame(kmeans_sim)));
colnames(kmeans_sim) <- c("Gene", "P.Value");
print("done kmeans");
coxReg_sim <- sapply(rownames(simExprs), FUN= coxReg, simObj);
coxReg_sim <- data.frame(t(data.frame(coxReg_sim)));
colnames(coxReg_sim) <- c("Gene", "P.Value");
print("done cox");

qCut50_sim <- sapply(rownames(simExprs), FUN=quantCutSA, simObj, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_sim <- data.frame(t(data.frame(qCut50_sim)));
colnames(qCut50_sim) <- c("Gene", "P.Value");
print("done median cut");

qCut2575_sim <- sapply(rownames(simExprs), FUN=quantCutSA, simObj, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_sim <- data.frame(t(data.frame(qCut2575_sim)));
colnames(qCut2575_sim) <- c("Gene", "P.Value");
print("done quantile cut");

#kmScan_sim <- sapply(rownames(simExprs), FUN=kapmPlot, simObj, F, tVar="TimeVar", eVar="eventVar");
#kmScan_sim <- data.frame(t(data.frame(kmScan_sim)));
#colnames(kmScan_sim) <- c("Gene", "P.Value", "Adj.P.Value");

coxReg_sim <- coxReg_sim[rownames(kmeans_sim),];
qCut50_sim <- qCut50_sim[rownames(kmeans_sim),];
qCut2575_sim <- qCut2575_sim[rownames(kmeans_sim),];
#kmScan_sim <- kmScan_sim[rownames(kmeans_sim),];

#myDF <- data.frame(kmeans_sim, coxReg_sim[2], qCut50_sim[2], qCut2575_sim[2], kmScan_sim[2:3]);
#colnames(myDF) <- c("gene", "kmeans", "coxreg", "qcut50", "qcut2575", "kmScanP", "kmScanQ"); 
myDF <- data.frame(kmeans_sim, coxReg_sim[2], qCut50_sim[2], qCut2575_sim[2]);
colnames(myDF) <- c("gene", "kmeans", "coxreg", "qcut50", "qcut2575"); 

return(myDF)
}

#ROC frame
createROCFrame <- function(data, myCol, cList)
{
    data <- data[myCol]
    data[,1] <- (-1)*log10(data[,1]);
    data[,"gene"] <- rownames(data);
    data[,"symbol"] <- rownames(data);
    data[,"labs"]<- as.numeric(data[,"gene"]%in%cList)
    data <- data[,c(1,4)];
    colnames(data)[1] <- "preds";
    return(data);
}

#Write no noise
write.table(data.frame(posControlGenes), "PositiveControls.txt", sep="\t", row.names=F);
posControlList <- gsub("\\|", ".", rownames(simExprs)[posControlGenes])

#Run no noise
ResNoNoise <- CreateMatrix(simObjNoNoise);
ResNoNoise[,3] <- as.numeric(as.character(ResNoNoise[,3]));
ResNoNoise[,4] <- as.numeric(as.character(ResNoNoise[,4]));
ResNoNoise[,5] <- as.numeric(as.character(ResNoNoise[,5]));
ResNoNoise[,2] <- as.numeric(as.character(ResNoNoise[,2]));
noNoiseDF <- rbind(data.frame(createROCFrame(ResNoNoise, 3, posControlList), method="Cox Regression"),
data.frame(createROCFrame(ResNoNoise, 2, posControlList), method="K-Means"),
data.frame(createROCFrame(ResNoNoise, 5, posControlList), method="Quantile 25th-75th"),
data.frame(createROCFrame(ResNoNoise, 4, posControlList), method="Median"))
noNoiseROC <- roconMult(noNoiseDF, myTitle="ROC No Noise");
write.table(ResNoNoise, "ResNoNoise.txt", sep="\t", row.names=T);
print("done no noise");

#Add 05 % noise
simExprs05 <- round(abs(addNoise(.05)+simExprs)); 
simObj05Noise <- list(simExprs05, simMeta);
Res05Noise <- CreateMatrix(simObj05Noise);
Res05Noise[,3] <- as.numeric(as.character(Res05Noise[,3]));
Res05Noise[,4] <- as.numeric(as.character(Res05Noise[,4]));
Res05Noise[,5] <- as.numeric(as.character(Res05Noise[,5]));
Res05Noise[,2] <- as.numeric(as.character(Res05Noise[,2]));
Noise05PercDF <- rbind(data.frame(createROCFrame(Res05Noise, 3, posControlList), method="Cox Regression"),
data.frame(createROCFrame(Res05Noise, 2, posControlList), method="K-Means"),
data.frame(createROCFrame(Res05Noise, 5, posControlList), method="Quantile 25th-75th"),
data.frame(createROCFrame(Res05Noise, 4, posControlList), method="Median"))
Noise05PercROC <- roconMult(Noise05PercDF, myTitle="ROC 05% Noise");
write.table(Res05Noise, "Res05Noise.txt", sep="\t", row.names=T);
print("done 30 noise");


#Add 7.5 % noise
simExprs075 <- round(abs(addNoise(.075)+simExprs)); 
simObj075Noise <- list(simExprs075, simMeta);
Res075Noise <- CreateMatrix(simObj075Noise);
Res075Noise[,3] <- as.numeric(as.character(Res075Noise[,3]));
Res075Noise[,4] <- as.numeric(as.character(Res075Noise[,4]));
Res075Noise[,5] <- as.numeric(as.character(Res075Noise[,5]));
Res075Noise[,2] <- as.numeric(as.character(Res075Noise[,2]));
Noise075PercDF <- rbind(data.frame(createROCFrame(Res075Noise, 3, posControlList), method="Cox Regression"),
data.frame(createROCFrame(Res075Noise, 2, posControlList), method="K-Means"),
data.frame(createROCFrame(Res075Noise, 5, posControlList), method="Quantile 25th-75th"),
data.frame(createROCFrame(Res075Noise, 4, posControlList), method="Median"))
Noise075PercROC <- roconMult(Noise075PercDF, myTitle="ROC 50% Noise");
write.table(Res50Noise, "Res50Noise.txt", sep="\t", row.names=T);
print("done 50 noise");

#Add 10 % noise
simExprs10 <- round(abs(addNoise(.1)+simExprs)); 
simObj10Noise <- list(simExprs10, simMeta);
Res10Noise <- CreateMatrix(simObj10Noise);
Res10Noise[,3] <- as.numeric(as.character(Res10Noise[,3]));
Res10Noise[,4] <- as.numeric(as.character(Res10Noise[,4]));
Res10Noise[,5] <- as.numeric(as.character(Res10Noise[,5]));
Res10Noise[,2] <- as.numeric(as.character(Res10Noise[,2]));
Noise10PercDF <- rbind(data.frame(createROCFrame(Res10Noise, 3, posControlList), method="Cox Regression"),
data.frame(createROCFrame(Res10Noise, 2, posControlList), method="K-Means"),
data.frame(createROCFrame(Res10Noise, 5, posControlList), method="Quantile 25th-75th"),
data.frame(createROCFrame(Res10Noise, 4, posControlList), method="Median"))
Noise10PercROC <- roconMult(Noise10PercDF, myTitle="ROC 10% Noise");
write.table(Res10Noise, "Res10Noise.txt", sep="\t", row.names=T);
print("done 10 noise");

#Add 15 % noise
simExprs15 <- round(abs(addNoise(.15)+simExprs)); 
simObj15Noise <- list(simExprs15, simMeta);
Res15Noise <- CreateMatrix(simObj15Noise);
Res15Noise[,3] <- as.numeric(as.character(Res15Noise[,3]));
Res15Noise[,4] <- as.numeric(as.character(Res15Noise[,4]));
Res15Noise[,5] <- as.numeric(as.character(Res15Noise[,5]));
Res15Noise[,2] <- as.numeric(as.character(Res15Noise[,2]));
Noise15PercDF <- rbind(data.frame(createROCFrame(Res15Noise, 3, posControlList), method="Cox Regression"),
data.frame(createROCFrame(Res15Noise, 2, posControlList), method="K-Means"),
data.frame(createROCFrame(Res15Noise, 5, posControlList), method="Quantile 25th-75th"),
data.frame(createROCFrame(Res15Noise, 4, posControlList), method="Median"))
Noise15PercROC <- roconMult(Noise15PercDF, myTitle="ROC 15% Noise");
write.table(Res15Noise, "Res10Noise.txt", sep="\t", row.names=T);
print("done 15 noise");

#Add 40 % noise
simExprs40 <- round(abs(addNoise(.4)+simExprs)); 
simObj40Noise <- list(simExprs40, simMeta);
Res40Noise <- CreateMatrix(simObj40Noise);
Res40Noise[,3] <- as.numeric(as.character(Res40Noise[,3]));
Res40Noise[,4] <- as.numeric(as.character(Res40Noise[,4]));
Res40Noise[,5] <- as.numeric(as.character(Res40Noise[,5]));
Res40Noise[,2] <- as.numeric(as.character(Res40Noise[,2]));
Noise40PercDF <- rbind(data.frame(createROCFrame(Res40Noise, 3, posControlList), method="Cox Regression"),
data.frame(createROCFrame(Res40Noise, 2, posControlList), method="K-Means"),
data.frame(createROCFrame(Res40Noise, 5, posControlList), method="Quantile 25th-75th"),
data.frame(createROCFrame(Res40Noise, 4, posControlList), method="Median"))
Noise40PercROC <- roconMult(Noise40PercDF, myTitle="ROC 15% Noise");
write.table(Res40Noise, "Res40Noise.txt", sep="\t", row.names=T);
print("done 40 noise");



noNoiseROC[[3]]
Noise05PercROC[[3]]
Noise075PercROC[[3]]
Noise10PercROC[[3]]
Noise15PercROC[[3]]



#Add 5 % noise
#simExprs05 <- round(abs(addNoise(.05)+simExprs)); 
#simObj05Noise <- list(simExprs05, simMeta);
#Res05Noise <- CreateMatrix(simObj05Noise);
#write.table(Res05Noise, "Res05Noise.txt", sep="\t", row.names=T);
#print("done 5 noise");


#Add 20 % noise
#simExprs20 <- round(abs(addNoise(.2)+simExprs)); 
#simObj20Noise <- list(simExprs20, simMeta);
#Res20Noise <- CreateMatrix(simObj20Noise);
#write.table(Res20Noise, "Res20Noise.txt", sep="\t", row.names=T);
#print("done 20 noise");




#save.image("../data/ALLDATAFin2.RData");




