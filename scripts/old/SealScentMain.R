# Hauptskript Datananalyse Seal Scent Project  MUMS

## packages
library("ggplot2")
library("grid")
library("vegan")
library("MASS") #required for vegan
## functions
source("PCAScores.R")
source("FAScores.R")
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
Sub <- SubsetAll("Mums",scent.abundance,scent.diversity,relatedness,heterozygosity,pca.scores) # out:list 

# from list back to original formats
abund.mums <- Sub[[1]]
div.mums <- Sub[[2]]
relate.mums <- Sub[[3]]
het.mums <- Sub[[4]]
pca.scores.mums <- Sub[[5]]


## get pairwise pca differences and  relatedness as vectors or diff.pc matrix

# diff.pc.mat.mums <- PCADiff(relate.mums,pca.scores.mums,num.pca=6, df=TRUE)
diff.relate.mums <- PCADiff(relate.mums,pca.scores.mums,num.pca=1)

        # pairwise pc differences as vectors // bad code here (just transforming around)
        relate.vec.mums <- subset(diff.relate.mums,select = relatedness)
        genetic.dist.vec.mums <- 1-relate.vec.mums # introduce distance instead relatedness
        diff.pca.vec.mums <- subset(diff.relate.mums,select = pca)

# create data.frame for relatedness plotting in ggplot //Main data.fame for relatedness here
relate.mums.df <- data.frame("genetic.distance" = genetic.dist.vec.mums[, 1],"pca.diff"=diff.pca.vec.mums[, 1])

## Mantel tests (pearson) for all pcs and bray curtis matrix of overall scent abundance(PC0)
results.relatedness.mums <- ResultsRelatedness(relate.mums,pca.scores.mums, abund.mums)

# create data.frame for heterozygosity plotting in ggplot //main data.frame for heterozygosity here
het.mums.df <- data.frame("het"= het.mums[,1],row.names=rownames(het.mums))
het.mums.df$NumComps <- div.mums$S
het.mums.df[, 3:12] <- pca.scores.mums[, 1:10]

# PC vs. Heterozygosity
results.het.mums <- resultsHet(het.mums.df)

# create pca.scores data.frame with factor beaches
# pca.scores$beach <- factor(factors$Beach,labels=c("beach 1","beach 2")) # bloody hell, this is cool!!!

# delete unnecessary objects
rm(scent.abundance,scent.diversity,relatedness,heterozygosity,Sub)
rm(abund.mums, div.mums,factors, het.mums,pca.scores,relate.mums,pca.scores.mums)
rm(diff.relate.mums,relate.vec.mums,genetic.dist.vec.mums,diff.pca.vec.mums)