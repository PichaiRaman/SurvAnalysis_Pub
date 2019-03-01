####################################
#Code to format and clean TCGA Data
#RNA-Seq-V2 in a folder called expr
#Clinical annotation in a folder called annot
#
#Pichai Raman
#9/17/2015
###################################




#1. First let's get all the metadata and format

clinFile <- list.files("annot/Clinical/Biotab/")[grep("org_clinical_patient", list.files("annot/Clinical/Biotab/"))];
annDat <- read.delim(paste("annot/Clinical/Biotab/", clinFile, sep=""));
annDat <- annDat[,c("bcr_patient_barcode", "last_contact_days_to", "death_days_to", "vital_status")];
annDat <- annDat[-1:-2,];


########################
#Format days column
########################
annDat[,"death_days_to"] <- as.character(annDat[,"death_days_to"]);
annDat[annDat[,"death_days_to"]=="[Not Applicable]","death_days_to"]<- "0";
annDat[,"death_days_to"] <- as.numeric(annDat[,"death_days_to"]);
annDat[,"last_contact_days_to"] <- as.character(annDat[,"last_contact_days_to"]);
annDat[annDat[,"last_contact_days_to"]=="[Not Available]","last_contact_days_to"]<- "0";
annDat[,"last_contact_days_to"] <- as.numeric(annDat[,"last_contact_days_to"]);
annDat[,"TimeVar"] <- ifelse(annDat[,"death_days_to"]>annDat[,"last_contact_days_to"], annDat[,"death_days_to"], annDat[,"last_contact_days_to"]);


########################
#Format vital_status column
########################
annDat[,"eventVar"] <- ifelse(annDat[,"vital_status"]=="Alive", 0, 1);


#Finally remove all NA
annDat <- na.omit(annDat);
rownames(annDat) <- annDat[,1];





#2. Now let's create RNA-Seq Matrix

########################
#List all files
########################
rPath <- "expr/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/";
rnaFiles <- list.files(rPath);

########################
#Pull out just the normalized gene values
########################
rnaFiles <- rnaFiles[grep("rsem.genes.normalized_results", rnaFiles)]
print(length(rnaFiles));

########################
#Read in files and create RNA-Seq Matrix
########################
allDat <- read.delim(paste(rPath, rnaFiles[1], sep="")); 
for(i in 2:length(rnaFiles))
{
tmpDat <- read.delim(paste(rPath, rnaFiles[i], sep="")); 
allDat <- cbind(allDat, tmpDat[2]);
}
rownames(allDat) <- allDat[,1];
allDat <- allDat[-1];
colnames(allDat) <- rnaFiles;

########################
#Replace colnames with TCGA Identifier
########################
mapFile <- read.delim("expr/FILE_SAMPLE_MAP.txt");
rownames(mapFile) <- mapFile[,1];
mapFile <- mapFile[colnames(allDat),];
colnames(allDat) <- substring(as.character(mapFile[,2]), 1, 12);


#3. Finally make sure we have same samples in both

intersectSamps <- intersect(colnames(allDat), rownames(annDat));

allDat <- allDat[,intersectSamps];
annDat <- annDat[intersectSamps,];


write.table(annDat, "annot.txt", sep="\t", row.names=T);
write.table(allDat, "exprs.txt", sep="\t", row.names=T);











