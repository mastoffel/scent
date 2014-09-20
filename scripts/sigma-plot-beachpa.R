# get nMDS plot data in the right format for sigma plot

# raw data, just two columns (x, y)
beachdata <- read.csv("beachpa.csv", header = F, row.names = NULL)
names(beachdata) <- c("x", "y")

row.names(beachdata) <- row.names(factors)

#add Family factor from factors
beachdata$family <- factors$Family

# sort by family
library(dplyr)
beach <- arrange(beachdata, family)

beach <- beach[, -3]

# short for loop
data <- matrix(NA, nrow=2, ncol=82)
beach <- as.matrix(beach)

for (i in seq(1,(ncol(data)-1), by = 2)) {
        data[1:2, i:(i+1)] <- beach[i:(i+1), 1:2]
}