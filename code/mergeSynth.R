############################################
# Author: Pichai Raman
# Function: code to merge in silico results
############################################

# 0
data <- read.delim("../data/SyntheticData/Res_0Noise.txt")
dataJM <- read.delim("../data/SurvResultsJM/survivalData/sim.0.txt")
data <- merge(data, dataJM[1:2], by.x = "gene", by.y = "genes", all.x = T)
data[is.na(data)] <- 1
row.names(data) <- data[,1]
colnames(data)[10] <- "dmod"
write.table(data, "../data/SyntheticDataMerge/Res_0Noise.txt", sep = "\t", row.names = T)

# .1
data <- read.delim("../data/SyntheticData/Res_0.1Noise.txt")
dataJM <- read.delim("../data/SurvResultsJM/survivalData/sim.0.1.txt")
data <- merge(data, dataJM[1:2], by.x = "gene", by.y = "genes", all.x = T)
data[is.na(data)] <- 1
row.names(data) <- data[,1]
colnames(data)[10] <- "dmod"
write.table(data, "../data/SyntheticDataMerge/Res_0.1Noise.txt", sep = "\t", row.names = T)

# .25
data <- read.delim("../data/SyntheticData/Res_0.25Noise.txt")
dataJM <- read.delim("../data/SurvResultsJM/survivalData/sim.0.25.txt")
data <- merge(data, dataJM[1:2], by.x = "gene", by.y = "genes", all.x = T)
data[is.na(data)] <- 1
row.names(data) <- data[,1]
colnames(data)[10] <- "dmod"
write.table(data, "../data/SyntheticDataMerge/Res_0.25Noise.txt", sep = "\t", row.names = T)

# .5
data <- read.delim("../data/SyntheticData/Res_0.5Noise.txt")
dataJM <- read.delim("../data/SurvResultsJM/survivalData/sim.0.5.txt")
data <- merge(data, dataJM[1:2], by.x = "gene", by.y = "genes", all.x = T)
data[is.na(data)] <- 1
row.names(data) <- data[,1]
colnames(data)[10] <- "dmod"
write.table(data, "../data/SyntheticDataMerge/Res_0.5Noise.txt", sep = "\t", row.names = T)

# .75
data <- read.delim("../data/SyntheticData/Res_0.75Noise.txt")
dataJM <- read.delim("../data/SurvResultsJM/survivalData/sim.0.75.txt")
data <- merge(data, dataJM[1:2], by.x = "gene", by.y = "genes", all.x = T)
data[is.na(data)] <- 1
row.names(data) <- data[,1]
colnames(data)[10] <- "dmod"
write.table(data, "../data/SyntheticDataMerge/Res_0.75Noise.txt", sep = "\t", row.names = T)

# 1
data <- read.delim("../data/SyntheticData/Res_1Noise.txt")
dataJM <- read.delim("../data/SurvResultsJM/survivalData/sim.1.txt")
data <- merge(data, dataJM[1:2], by.x = "gene", by.y = "genes", all.x = T)
data[is.na(data)] <- 1
row.names(data) <- data[,1]
colnames(data)[10] <- "dmod"
write.table(data, "../data/SyntheticDataMerge/Res_1Noise.txt", sep = "\t", row.names = T)

# 1.5
data <- read.delim("../data/SyntheticData/Res_1.5Noise.txt")
dataJM <- read.delim("../data/SurvResultsJM/survivalData/sim.1.5.txt")
data <- merge(data, dataJM[1:2], by.x = "gene", by.y = "genes", all.x = T)
data[is.na(data)] <- 1
row.names(data) <- data[,1]
colnames(data)[10] <- "dmod"
write.table(data, "../data/SyntheticDataMerge/Res_1.5Noise.txt", sep = "\t", row.names = T)

