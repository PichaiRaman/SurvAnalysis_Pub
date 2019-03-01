##########################################################
<<<<<<< HEAD
#CHOP Network sucks so created
#a function to download GSE's
#one GSM at a time. You just need to enter the GSE
#name and all the GSM's will be downloaded
#as well as a file with the annotation
#it will download to the directory you're in currently ...
#Pichai Raman
##########################################################

#Call libraries
=======
#Code to download all data and uncompress
#GSE425, GSE1159, GSE14468, GSE30258, GSE30285, and GSE37642
##########################################################



#Call libraries
library("GEOquery");

#Main function that will download all GSM's into a directory
#given a particular GSE
DownloadGSE <- function(gseName=null, myWD=".")
{
    setwd(myWD);
    system(paste("mkdir ", gseName, sep=""));
    setwd(paste(gseName));
    #Read in files to download
    myGSE <- getGEO(gseName,GSEMatrix=T);
    sampleAnnotation <- pData(phenoData(myGSE[[1]]));
    sampleLoc <- as.character(sampleAnnotation[,"supplementary_file"]);

    #Download File
    for(i in 1:length(sampleLoc))
    {
    system(paste("wget ", sampleLoc[i], sep=""));
    }
    write.table(sampleAnnotation, "metaDataSample.txt", sep="\t", row.names=F);

    #Okay now let's check how many GSM's there are
    numFiles <- length(list.files())-1;
    numSamples <- dim(sampleAnnotation)[1];
    check="Failure"
    if(numFiles==numSamples)
    {
        check="Success";
    }
    print(paste("There are ", numFiles, " files downloaded and ", numSamples, " samples in your annotation file so your download was a ", check, sep=""));
}

DownloadGSE("GSE1159",myWD="../data/"); #Done

DownloadGSE("GSE37642",myWD="../data/"); #Done

DownloadGSE("GSE14468", myWD="/data/SurvAnalysis_Pub/data");

DownloadGSE("GSE30258", myWD="/data/SurvAnalysis_Pub/data");

DownloadGSE("GSE30285", myWD="/data/SurvAnalysis_Pub/data");








>>>>>>> 9257864f2a990c53e086a60e1c59a2b6a56cdbc5





<<<<<<< HEAD
list <- read.delim("filelist.txt");


lname <- list[list[,1]=="File",];
lname <- as.character(list[,"Name"]);
lname <- lname[-grep("GSE", lname)];
lname <- gsub(".CEL.gz", "", lname);
lname2 <- substring(lname, 1, 9);


dFile <- function(fName)
{
    fnnn <- gsub("GSM", "", fName);
    fnnn <- substr(fnnn, 1, nchar(fnnn)-3);
    fnnn <- paste("GSM", fnnn, "nnn", sep="");
    myDir <- paste("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/", fnnn, "/", fName, "/suppl/", sep="");
    
    fName2 <-strsplit(getURL(myDir), "\n")[[1]]
    fName2 <- substr(fName2 , str_locate(fName2, "GSM")[[1]], nchar(fName2));
    dld <- paste(myDir, fName2, sep="");
    print(dld);
    download.file(dld, paste(fName, ".CEL.gz",sep=""));
    
}


lapply(lname2, FUN=dFile);
=======






>>>>>>> 9257864f2a990c53e086a60e1c59a2b6a56cdbc5
