# Hauptskript Datananalyse Seal Scent Project  PUPS

## packages
library("ggplot2")
library("grid")
library("vegan")
library("MASS") #required for vegan
## functions
source("PCAScores.R")
source("PCADiff.R")
source("SubsetAll.R")
source("ResultsRelatedness.R")
source("resultsHet.R")
source("multiplot.R")


## loading data

# already standardized and transformed, transposed abundance matrix
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 

# 6 different diversity indices (from primer)
scent.diversity <- read.csv(".\\csv_files\\scent diversity.csv",row.names=1) 

# relatedness matrix
relatedness <- read.csv(".\\csv_files\\relatedness_41loci.csv",row.names=1)

# heterozygosity vector SH
heterozygosity <- read.csv(".\\csv_files\\heterozygosity_41loci.csv", row.names=1) 

# beach and family vector
factors <- read.csv(".\\csv_files\\factors.csv",row.names=1) 


## computing PCA with vegan package (10 Principal Components) //

# pca on the whole dataset and then subsetted
pca.scores <- PCAScores(scent.abundance)


## subsetting all data objects

Sub <- SubsetAll("Pups",scent.abundance,scent.diversity,relatedness,heterozygosity,pca.scores) # out:list 

# from list back to original formats
abund.pups <- Sub[[1]]
div.pups <- Sub[[2]]
relate.pups <- Sub[[3]]
het.pups <- Sub[[4]]
pca.scores.pups <- Sub[[5]]


## get pairwise pca differences and  relatedness as vectors or diff.pc matrix

# diff.pc.mat.mums <- PCADiff(relate.mums,pca.scores.mums,num.pca=6, df=TRUE)
diff.relate.pups <- PCADiff(relate.pups,pca.scores.pups,num.pca=6)

# pairwise pc differences as vectors // bad code here (just transforming around)
relate.vec.pups <- subset(diff.relate.pups,select = relatedness)
genetic.dist.vec.pups <- 1-relate.vec.pups # introduce distance instead relatedness
diff.pca.vec.pups <- subset(diff.relate.pups,select = pca)

# create data.frame for relatedness plotting in ggplot //Main data.fame for relatedness here
relate.pups.df <- data.frame("genetic.distance" = genetic.dist.vec.pups[, 1],"pca.diff"=diff.pca.vec.pups[, 1])

## Mantel tests (pearson) for all pcs
results.relatedness.pups <- ResultsRelatedness(relate.pups,pca.scores.pups,abund.pups)

# create data.frame for heterozygosity plotting in ggplot //main data.frame for heterozygosity here
het.pups.df <- data.frame("het"= het.pups[,1],row.names=rownames(het.pups))
het.pups.df$NumComps <- div.pups$S
het.pups.df[, 3:12] <- pca.scores.pups[, 1:10]

# PC vs. Heterozygosity
results.het.pups <- resultsHet(het.pups.df)

# create pca.scores data.frame with factor beaches
# pca.scores$beach <- factor(factors$Beach,labels=c("beach 1","beach 2")) # bloody hell, this is cool!!!

# delete unnecessary objects
rm(scent.abundance,scent.diversity,relatedness,heterozygosity,Sub)
rm(abund.pups, div.pups,factors, het.pups,pca.scores,relate.pups,pca.scores.pups)
rm(diff.relate.pups,relate.vec.pups,genetic.dist.vec.pups,diff.pca.vec.pups)