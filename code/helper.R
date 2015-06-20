##########################################
#Objective : Function to help with code
#
##########################################

#Turn factor into numeric
factToNum <- function(x)
{
    x <- as.numeric(as.character(x));
}


#Z-score
myCV <- function(x) {sd(x)/mean(x)}
