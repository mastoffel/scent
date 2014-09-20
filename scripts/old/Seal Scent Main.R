# Hauptskript Datananalyse Seal Scent Project 

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
        relatedness <- read.csv(".\\csv_files\\relatedness.csv",row.names=1)
        
        # heterozygosity vector SH
        heterozygosity <- read.csv(".\\csv_files\\heterozygosity.csv", row.names=1) 
        
        # beach and family vector
        factors <- read.csv(".\\csv_files\\factors.csv",row.names=1) 


## computing PCA with vegan package (10 Principal Components) //
  
        # pca on the whole dataset and then subsetted
        pca.scores <- PCAScores(scent.abundance)
      
        
## subsetting all data objects
        
        Sub <- SubsetAll("All",scent.abundance,scent.diversity,relatedness,heterozygosity,pca.scores) # out:list 
        
        # from list back to original formats
        abund <- Sub[[1]]
        div <- Sub[[2]]
        relate <- Sub[[3]]
        het <- Sub[[4]]
        pca.scores <- Sub[[5]]
        
        #rm(scent.abundance,scent.diversity,relatedness,heterozygosity)
        
## get pairwise pca differences and  relatedness as vectors or diff.pc matrix
        
        diff.pc.mat <- PCADiff(relate,pca.scores,num.pca=1, df=TRUE)
        diff.relate <- PCADiff(relate,pca.scores,num.pca=2)
        
## glm´s for all pc´s 
        results.relatedness.pups <- ResultsRelatedness(relate,pca.scores)
        
        
# pairwise pc differences as vectors
        relate.vec <- subset(diff.relate,select = relatedness)
        genetic.dist.vec <- 1-relate.vec # introduce distance instead relatedness
        diff.pca.vec <- subset(diff.relate,select = pca)
        
        # rm(diff.relate)
        
        All.Pairs <- data.frame("genetic.distance" = genetic.dist.vec[, 1],"pca.diff"=diff.pca.vec[, 1])
        
        # relate.diff.plot <- qplot(relate,pca.diff,data=All.Pairs) +
                # geom_smooth(method="lm")
        
        # fit <- lm(All.Pairs$relate ~ All.Pairs$pca.diff)
        # summary(fit)
        
        
# create data.frame with all vars for heterozygosity analysis in
        All <- data.frame("het"= het[,1],row.names=rownames(het))
        All$NumComps <- div$S
        All[, 3:12] <- pca.scores[, 1:10]
       
        # All$Shannon <- div[,5]
        # All$Simpson <- div[,6]

# PC vs. Heterozygosity
        results.het <- resultsHet(All)
        
# create pca.scores data.frame with factor beaches
        pca.scores$beach <- factor(factors$Beach,labels=c("beach 1","beach 2")) # bloody hell, this is cool!!!
        

        


