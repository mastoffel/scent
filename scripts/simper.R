# simper analysis for each beach separately for mum-pup similarity

library(ggplot2)
library(vegan)

# mum pup 
scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "projects\\sealscent\\data_files\\",
                                                  "Rdata\\csv_files\\scent abundances.csv", 
                                                   sep = ""), row.names=1)))

# relatedness matrix (old: relatedness_41loci.csv)
relatedness <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                              "projects\\sealscent\\data_files\\",
                              "Rdata\\csv_files\\",
                              "relatednessnew.csv", sep = ""),
                               row.names=1)

## heterozygosity SH
heterozygosity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                 "projects\\sealscent\\data_files\\",
                                 "Rdata\\csv_files\\",
                                 "heterozygosity_41loci.csv", sep = ""),
                                  row.names=1) 

# beach and family factor
factors <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                          "projects\\sealscent\\data_files\\",
                          "Rdata\\csv_files\\",
                          "factors.csv", sep = ""),
                           row.names=1) 


# beach simper results
options(scipen=8)

beach <- as.factor(factors$Beach)
simper_beach <- simper(scent_abundance, beach)

simper_beach_names <- rownames(summary(simper_beach)[[1]])
simper_beach_ind <- which(names(scent_abundance) %in% simper_beach_names[1:15])
contribution <- summary(simper_beach)[[1]]$contr[1:15]


source("extract_simper.R")

# mum-pup simper results from primer
simp <- extract_simper("allsimper.csv")

# best substances from all samples----------------------------------------------
# new idea: sorting substances by the average variance explained

toptwo <- sapply(simp, function(x) x$comp[1:2])
topcom <-  sapply(simp, function(x) x$comp[1:length(comp)])
vars <- sapply(simp, function(x) x$var[1:length(var)])

library(plyr)
library(dplyr)

# all best simper substances plus explained variances in one data.frame
allsimp <- ldply(simp, data.frame)

allsimp$comp <- as.factor(allsimp$comp)
allsimp$var <- as.numeric(allsimp$comp)

simpsum <- allsimp %>%
                group_by(comp) %>%
                summarise(
                meanvar = mean(var, na.rm = TRUE),
                  meansd = sd(var, na.rm = TRUE))


hist(table(toptwo), breaks=c(1:41))

comps <- sort(table(toptwo), decreasing = TRUE)

# getting top compounds
topcomp <- names(sort(table(toptwo), decreasing = TRUE)[1:9])
top_mp_ind <- which(names(scent_abundance) %in% topcomp)

# top comps
top_mp_ind <- c(58, 68,  86,  90, 106, 107, 164)

