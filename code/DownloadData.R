##########################################################
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

#DownloadGSE("GSE1159",myWD="../data/"); #Done

#DownloadGSE("GSE37642",myWD="../data/"); #Done

DownloadGSE("GSE14468", myWD="/data/SurvAnalysis_Pub/data");

DownloadGSE("GSE30258", myWD="/data/SurvAnalysis_Pub/data");

DownloadGSE("GSE30285", myWD="/data/SurvAnalysis_Pub/data");

#DownloadGSE("GSE425",myWD="../data/");



















