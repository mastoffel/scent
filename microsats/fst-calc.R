# gene files

library(StAMPP)
library(hierfstat)
library(pegas)

# 41 loci genotype file
genotypes <- read.table(".\\txt\\raw_41loci_ordered.txt", na.strings = "NA", row.names = 1)

# factor
factors <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                          "projects\\sealscent\\data_files\\",
                          "Rdata\\csv_files\\",
                          "factors.csv", sep = ""),
                        row.names=1) 

gene <- cbind(factors$Beach, genotypes)
names(gene)[1] <- "Beach"
library(dplyr)

gene <- gene[order(gene$Beach), ]
