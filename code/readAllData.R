######################################
# Author: Pichai Raman 
# Function: Code to read all data and 
# put it in an R object for 
# further processing
# Date: 12/11/2015
######################################

# ovarian
annot_ov <- read.delim("../data/raw/ovca/annot.txt");
exprs_ov <- read.delim("../data/raw/ovca/exprs.txt")
ov <- list(exprs_ov, annot_ov);

# prostate
annot_pr <- read.delim("../data/raw/prca/annot.txt");
exprs_pr <- read.delim("../data/raw/prca/exprs.txt")
pr <- list(exprs_pr, annot_pr);

# head and neck
annot_hn <- read.delim("../data/raw/hnca/annot.txt");
exprs_hn <- read.delim("../data/raw/hnca/exprs.txt")
hn <- list(exprs_hn, annot_hn);

# kidney
annot_ki <- read.delim("../data/raw/kica/annot.txt");
exprs_ki <- read.delim("../data/raw/kica/exprs.txt")
ki <- list(exprs_ki, annot_ki);

save.image("../data/ALLDATA.RData");


















