######################################################
# Author: Code to merge data files into one data frame
# Function: Pichai Raman
# Date: 01/30/2018
######################################################

# folder with primary results
primRes <- "../data/SurvAnalysisResults"

# folder with disc/validation results
reproRes <- "../data/SurvAnalysisResultsRepro"

# list all data files
myFiles <- list.files(primRes)

AllRes <- read.delim(paste(primRes, "/", myFiles[1], sep = ""))
for(i in 2:length(myFiles)) {
  tmpFile <- read.delim(paste(primRes, "/", myFiles[i], sep = ""))
  if(grepl("kmScan", myFiles[i])) {
    tmpFile <- tmpFile[,c(1,3)]
  }
  colnames(tmpFile)[1] <- c("Gene")
  AllRes <- merge(AllRes, tmpFile, by.x = "Gene", by.y = "Gene")
}
rownames(AllRes) <- AllRes[,1]
AllRes <- AllRes[-1]
colnames(AllRes) <- gsub(".txt", "", myFiles)
write.table(AllRes, "../data/MergedData/resultsPrimary.txt", sep = "\t", row.names = T)

###########################################################
# Now do it for all samples with discovery and validation
###########################################################

# list all data files
myFiles <- list.files(reproRes)

AllRes <- read.delim(paste(reproRes, "/", myFiles[1], sep = ""))
for(i in 2:length(myFiles)) {
  tmpFile <- read.delim(paste(reproRes, "/", myFiles[i], sep = ""))
  if((dim(tmpFile)[2]) == 3) {
    tmpFile <- tmpFile[,c(1,3)]
  }
  colnames(tmpFile)[1] <- c("Gene")
  AllRes <- merge(AllRes, tmpFile, by.x="Gene", by.y="Gene")
}
rownames(AllRes) <- AllRes[,1]
AllRes <- AllRes[-1]
colnames(AllRes) <- gsub(".txt", "", myFiles)
write.table(AllRes, "../data/MergedData/resultsDiscValid.txt", sep = "\t", row.names = T)


