##########################################################
#CHOP Network sucks so created
#a function to download GSE's
#one GSM at a time. You just need to enter the GSE
#name and all the GSM's will be downloaded
#as well as a file with the annotation
#it will download to the directory you're in currently ...
#Pichai Raman
##########################################################

#Call libraries





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