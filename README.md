# SurvAnalysis_Pub
Determining the best method for survival analysis in the case of continuous variables.

In order to do this we will look at a few different methods

1. Quantile based cut-off at the median
2. Quantile based cut-off at 25th and 75th percentile
3. Scanning approach - finding optimal cut-off
4. Cox Regression
5. K-means
6. Distrubtion Specific cut-off (Jess Method)

### Specifics on how to run the code & what the files are

* KaplanScan.R - Contains code to run KaplanScan
* quantileCut.R - Contains code to run Quantile based cutoff at Median and 25th/75th
* coxReg.R - Contains code to run Cox Regression model
* kmeansSA.R - Contains code to run the K-means model
* **Please note that the Distribution specific cut-off was supplied by Dr. Mar**
* driver.R - Runs all the methods on the 4 cancer datasets
* driver_split.R - Runs all methods on the split data
* insilicoDataMethCompare.R - Runs all methods on the in silico data
* readAllData.R - Can be substituted for the first part of the the driver script, reads in all the data. 

### These functions are for creating final figures and tables 

* FormatJMResults.R - This will format Jess Survival analysis for use in the analysis 
* mergeResults.R - This will merge all methods into two data frames, one for all data and one for the split discovery / validation
* compareResults.R - This creates some of the images an plots from the merged data 
* kidneyExplainCreateDF.R - This is a standalone script to explain why all methods seem to predict kidney survival markers well
* mergeSynth.R - Merges Jess's results and mine on the synthetic data
* synthDataSection.R - Generates Synthetic data plot
* checkMeanQuantileCutoff.R - Shows cut-offs are sensitive to composition of datasets

#### These files can be ignored KR

* DownloadData.R - Utility to download data
* createFiles.R - Utility to format data
* readAllData.R - Utility to read formatted data and save as R object
* createPCA.R - Creates PCA plot









you can run the code by going into the code directory and typing 


```Rscript driver.R ```




