##############################################
#Code to simulate RNA-Seq data and see
#how well methods can deal with noise 
#in the system
#
##############################################

#Call libraries
library("SimSeq");

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
simDat <- SimData(counts = round(as.matrix(exprs_hn)), treatment=annot_hn[,"eventVar"], n.genes=5000, n.diff=0, k.ind=100, sort.method="unpaired");

#Counts
simExprs <- data.frame(simDat$counts);
simExprs <- simExprs[,1:150];

#metaData
simMeta <- data.frame(colnames(simExprs), simDat$treatment[1:150]);
set.seed(100);

survTimes1 <- sort(annot_hn[annot_hn["eventVar"]==0, "TimeVar"])[(nrow(annot_hn[annot_hn["eventVar"]==0,])-99): nrow(annot_hn[annot_hn["eventVar"]==0,])];
survTimes2 <- sort(annot_hn[annot_hn["eventVar"]==1, "TimeVar"])[c(1:50)];
survTimes <- c(survTimes1, survTimes2);

simMeta[,"TimeVar"] <- survTimes;
colnames(simMeta)[1:2] <- c("Sample", "eventVar");

#Now Let's generate positive controls
genPosControlMat <- function(x)
{
tmpMean <- mean(x)
tmpSrv <- survTimes[101:150];
survTimesNorm <- c(rep(0,100),(max(tmpSrv)-tmpSrv)/(max(tmpSrv)-min(tmpSrv)))
multiplier <- round(runif(1, min=2, max=6))
xN <- round(multiplier*x*survTimesNorm);
}
#Get genes with only values > 100
posControlGenes <- rownames(simExprs)%in%sample(rownames(simExprs)[apply(simExprs, FUN=mean, MARGIN=1)>100], 250)

simExprs[posControlGenes,] <- simExprs[posControlGenes,]+t(apply(simExprs[posControlGenes,], FUN=genPosControlMat, MARGIN=1));


simObjNoNoise <- list(simExprs, simMeta);

addNoise <- function(myPerc)
{
myMeans <-as.numeric(apply(simExprs, FUN=mean, MARGIN=1))*myPerc;
mySD <-as.numeric(apply(simExprs, FUN=sd, MARGIN=1))*myPerc;

myDat <- rnorm(ncol(simExprs), mean=myMeans[1], sd=mySD[1])
for(i in 2:length(myMeans))
{
tmpMat <- rnorm(ncol(simExprs), mean=myMeans[i], sd=mySD[i]);
myDat <- rbind(myDat, tmpMat);
}
return(as.matrix(myDat));
}


CreateMatrix <- function(simObj)
{
kmeans_sim <- sapply(rownames(simExprs), FUN= kmeansSA, simObj, tVar="TimeVar", eVar="eventVar");
kmeans_sim <- data.frame(t(data.frame(kmeans_sim)));
colnames(kmeans_sim) <- c("Gene", "P.Value");

coxReg_sim <- sapply(rownames(simExprs), FUN= coxReg, simObj);
coxReg_sim <- data.frame(t(data.frame(coxReg_sim)));
colnames(coxReg_sim) <- c("Gene", "P.Value");

qCut50_sim <- sapply(rownames(simExprs), FUN=quantCutSA, simObj, F, quantLow=.50,  quantHigh=.50, tVar="TimeVar", eVar="eventVar");
qCut50_sim <- data.frame(t(data.frame(qCut50_sim)));
colnames(qCut50_sim) <- c("Gene", "P.Value");

qCut2575_sim <- sapply(rownames(simExprs), FUN=quantCutSA, simObj, F, quantLow=.25,  quantHigh=.75, tVar="TimeVar", eVar="eventVar");
qCut2575_sim <- data.frame(t(data.frame(qCut2575_sim)));
colnames(qCut2575_sim) <- c("Gene", "P.Value");

kmScan_sim <- sapply(rownames(simExprs), FUN=kapmPlot, simObj, F, tVar="TimeVar", eVar="eventVar");
kmScan_sim <- data.frame(t(data.frame(kmScan_sim)));
colnames(kmScan_sim) <- c("Gene", "P.Value", "Adj.P.Value");

coxReg_sim <- coxReg_sim[rownames(kmeans_sim),];
qCut50_sim <- qCut50_sim[rownames(kmeans_sim),];
qCut2575_sim <- qCut2575_sim[rownames(kmeans_sim),];
kmScan_sim <- kmScan_sim[rownames(kmeans_sim),];

myDF <- data.frame(kmeans_sim, coxReg_sim[2], qCut50_sim[2], qCut2575_sim[2], kmScan_sim[2:3]);
colnames(myDF) <- c("gene", "kmeans", "coxreg", "qcut50", "qcut2575", "kmScanP", "kmScanQ"); 
return(myDF)
}

#Run no noise
ResNoNoise <- CreateMatrix(simObjNoNoise);
write.table(ResNoNoise, "ResNoNoise.txt", sep="\t", row.names=T);

#Add 5 % noise
simExprs05 <- round(abs(addNoise(.05)+simExprs)); 
simObj05Noise <- list(simExprs05, simMeta);
Res05Noise <- CreateMatrix(simObj05Noise);
write.table(Res05Noise, "Res05Noise.txt", sep="\t", row.names=T);

#Add 10 % noise
simExprs10 <- round(abs(addNoise(.1)+simExprs)); 
simObj10Noise <- list(simExprs10, simMeta);
Res10Noise <- CreateMatrix(simObj10Noise);
write.table(Res10Noise, "Res10Noise.txt", sep="\t", row.names=T);

#Add 20 % noise
simExprs20 <- round(abs(addNoise(.2)+simExprs)); 
simObj20Noise <- list(simExprs20, simMeta);
Res20Noise <- CreateMatrix(simObj20Noise);
write.table(Res20Noise, "Res20Noise.txt", sep="\t", row.names=T);

#Add 30 % noise
simExprs30 <- round(abs(addNoise(.3)+simExprs)); 
simObj30Noise <- list(simExprs30, simMeta);
Res30Noise <- CreateMatrix(simObj30Noise);
write.table(Res30Noise, "Res30Noise.txt", sep="\t", row.names=T);

#Add 40 % noise
simExprs40 <- round(abs(addNoise(.4)+simExprs)); 
simObj40Noise <- list(simExprs40, simMeta);
Res40Noise <- CreateMatrix(simObj40Noise);
write.table(Res40Noise, "Res40Noise.txt", sep="\t", row.names=T);

#Add 50 % noise
simExprs50 <- round(abs(addNoise(.5)+simExprs)); 
simObj50Noise <- list(simExprs50, simMeta);
Res50Noise <- CreateMatrix(simObj50Noise);
write.table(Res50Noise, "Res50Noise.txt", sep="\t", row.names=T);

#Add 60 % noise
simExprs60 <- round(abs(addNoise(.6)+simExprs)); 
simObj60Noise <- list(simExprs60, simMeta);
Res60Noise <- CreateMatrix(simObj60Noise);
write.table(Res60Noise, "Res60Noise.txt", sep="\t", row.names=T);

#Add 70 % noise
simExprs70 <- round(abs(addNoise(.7)+simExprs)); 
simObj70Noise <- list(simExprs70, simMeta);
Res70Noise <- CreateMatrix(simObj70Noise);
write.table(Res70Noise, "Res70Noise.txt", sep="\t", row.names=T);

#Add 80 % noise
simExprs80 <- round(abs(addNoise(.8)+simExprs)); 
simObj80Noise <- list(simExprs80, simMeta);
Res80Noise <- CreateMatrix(simObj80Noise);
write.table(Res80Noise, "Res80Noise.txt", sep="\t", row.names=T);

#Add 90 % noise
simExprs90 <- round(abs(addNoise(.9)+simExprs)); 
simObj90Noise <- list(simExprs90, simMeta);
Res90Noise <- CreateMatrix(simObj90Noise);
write.table(Res90Noise, "Res90Noise.txt", sep="\t", row.names=T);

#Add 100 % noise
simExprs100 <- round(abs(addNoise(1)+simExprs)); 
simObj100Noise <- list(simExprs100, simMeta);
Res100Noise <- CreateMatrix(simObj100Noise);
write.table(Res100Noise, "Res100Noise.txt", sep="\t", row.names=T);

#Add 200 % noise
simExprs200 <- round(abs(addNoise(2)+simExprs)); 
simObj200Noise <- list(simExprs200, simMeta);
Res200Noise <- CreateMatrix(simObj200Noise);
write.table(Res200Noise, "Res200Noise.txt", sep="\t", row.names=T);

#Add 300 % noise
simExprs300 <- round(abs(addNoise(3)+simExprs)); 
simObj300Noise <- list(simExprs300, simMeta);
Res300Noise <- CreateMatrix(simObj300Noise);
write.table(Res300Noise, "Res300Noise.txt", sep="\t", row.names=T);


save.image("../data/ALLDATAFin2.RData");




