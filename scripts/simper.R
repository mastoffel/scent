# simper analysis for each beach separately for mum-pup similarity

library(ggplot2)
library(vegan)
# mum pup 
scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "MSc.Behaviour\\Research\\Seal Scent\\R code\\",
                                                  "data\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))
coord <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                        "MSc.Behaviour\\Research\\Seal Scent\\R code\\",
                        "data\\csv_files\\coordinates.csv",
                        sep = ""),row.names=1) 

factors <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                          "MSc.Behaviour\\Research\\Seal Scent\\",
                          "R code\\data\\csv_files\\",
                          "factors.csv", sep = ""),
                        row.names=1) 

heterozygosity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                 "MSc.Behaviour\\Research\\Seal Scent\\",
                                 "R code\\data\\csv_files\\",
                                 "heterozygosity_41loci.csv", sep = ""),
                           row.names=1) 


# get simper results in
source("extract_simper.R")
simp_one <- extract_simper("beach_one_simper.csv")
simp_two <- extract_simper("beach_two_simper.csv")

simp <- extract_simper("allsimper.csv")

# best two substances from every pair per beach---------------------------------
toptwo_b1 <- sapply(simp_one, function(x) x$comp[1:2])
table(toptwo_b1)
toptwo_b2 <- sapply(simp_two, function(x) x$comp[1:2])
table(toptwo_b2)

# best substances from all samples----------------------------------------------
topfive <- sapply(simp, function(x) x$comp[1:5])

hist(table(topfive), breaks=c(1:41))

comps <- sort(table(topfive), decreasing = TRUE)

# getting top compounds
topcomp <- names(sort(table(topfive), decreasing = TRUE)[1:9])
top_mp_ind <- which(names(scent_abundance) %in% topcomp)

# top comps
top_mp_ind <- c(58, 68,  86,  90, 106, 107, 164)

