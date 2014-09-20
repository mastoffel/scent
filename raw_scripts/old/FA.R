
library(HDMD) # provides factor analysis when D>>N
scent.fa <- factor.pa.ginv(scent.abundance, nfactors = 6, prerotate=T, rotate = "promax", scores = T)
row.names(scent.fa$scores) = row.names(scent.abundance)
#order loadings of first factor
scent.fa$loadings[order(scent.fa$loadings[,2]),]
scent.fa$scores

VSS.scree(scent.fa)

## loading data

# already standardized and transformed, transposed abundance matrix
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 
scent.diversity <- read.csv(".\\csv_files\\scent diversity.csv",row.names=1) 
relatedness <- read.csv(".\\csv_files\\relatedness_41loci.csv",row.names=1)
heterozygosity <- read.csv(".\\csv_files\\heterozygosity_41loci.csv", row.names=1) 
factors <- read.csv(".\\csv_files\\factors.csv",row.names=1) 

pca.scores <- PCAScores(scent.abundance)
pca.scores[, 1:6] <- scent.fa$scores[, 1:6]

## subsetting all data objects

Sub <- SubsetAll("Mums",scent.abundance,scent.diversity,relatedness,heterozygosity,pca.scores) # out:list 

# from list back to original formats
abund.mums <- Sub[[1]]
div.mums <- Sub[[2]]
relate.mums <- Sub[[3]]
het.mums <- Sub[[4]]
pca.scores.mums <- Sub[[5]]


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

rm(scent.diversity,relatedness,heterozygosity,Sub)
rm(abund.mums, div.mums,factors, het.mums,pca.scores,relate.mums,pca.scores.mums)
rm(diff.relate.mums,relate.vec.mums,genetic.dist.vec.mums,diff.pca.vec.mums)


