# gene files
library(hierfstat)
library(diveRsity)
library(StAMPP)
library(pegas)
# 41 loci genotype file
genotypes <- read.table(".\\txt\\raw_41loci_ordered.txt", na.strings = "NA", row.names = 1)

# factors
factors <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                          "MSc.Behaviour\\Research\\Seal Scent\\",
                          "R code\\data\\csv_files\\",
                          "factors.csv", sep = ""),
                          row.names=1) 

data(potato)
data(Big_data)
