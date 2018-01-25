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

 KaplanScan.R - Contains code to run KaplanScan
 quantileCut.R - Contains code to run Quantile based cutoff at Median and 25th/75th
 coxReg.R - Contains code to run Cox Regression model
 kmeansSA.R - Contains code to run the K-means model
 Rscript driver.R - Runs all the methods on the 4 cancer datasets




you can run the code by going into the code directory and typing 


```Rscript driver.R ```




