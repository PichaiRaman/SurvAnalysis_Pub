##############################
# Author: Pichai Raman 
# Function: Code to run PCA Analysis
# Date: 12/3/2015
##############################

#########################################
#Read in data DONE
print("Starting Read of Data");

exprs_ov <- read.delim("../data/raw/ovca/exprs.txt")
colnames(exprs_ov) <- gsub("\\.", "-", colnames(exprs_ov))
rownames(exprs_ov) <- gsub(" ", "", rownames(exprs_ov));

#prostate
exprs_pr <- read.delim("../data/raw/prca/exprs.txt")
colnames(exprs_pr) <- gsub("\\.", "-", colnames(exprs_pr));
rownames(exprs_pr) <- gsub(" ", "", rownames(exprs_pr));

#head and neck
exprs_hn <- read.delim("../data/raw/hnca/exprs.txt")
colnames(exprs_hn) <- gsub("\\.", "-", colnames(exprs_hn))
rownames(exprs_hn) <- gsub(" ", "", rownames(exprs_hn));

#Kidney
exprs_ki <- read.delim("../data/raw/kica/exprs.txt")
colnames(exprs_ki) <- gsub("\\.", "-", colnames(exprs_ki))
rownames(exprs_ki) <- gsub(" ", "", rownames(exprs_ki));

print("Finished Reading Data");
##########################################

#Sort all row anmes
exprs_ov <- exprs_ov[rownames(exprs_ov),];
exprs_pr <- exprs_pr[rownames(exprs_ov),];
exprs_hn <- exprs_hn[rownames(exprs_ov),];
exprs_ki <- exprs_ki[rownames(exprs_ov),];
bigDF <- cbind(exprs_ov, exprs_pr, exprs_hn, exprs_ki);

#Let's remove all rows with mean count < 100
myRowMean <- rowMeans(bigDF);
bigDF <- bigDF[myRowMean>100,]
bigDF <- log2(bigDF+1);
bigDF <- data.frame(bigDF);
cancCol <- c(rep("OV",ncol(exprs_ov)), rep("PR",ncol(exprs_pr)), rep("HN",ncol(exprs_hn)), rep("KI",ncol(exprs_ki)));

bigAnnot <- data.frame(colnames(bigDF), cancCol);
outPCA <- prcomp(bigDF)$rotation;

write.table(outPCA, "../results/PCA.txt", sep="\t", row.names=T);

write.table(bigAnnot, "../results/PCA_Annot.txt", sep="\t", row.names=T);




