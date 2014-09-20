ScentResultsNoScores <- function(subgroup, type, plotPC) {
        ## gives a summary of the main results
        ## subgroup = "mums","pups", "all" 
        ## type = "pca", "fa"
        ## plotPC = factornumber to compute for plotting
        
        ## packages
        library("ggplot2")
        library("grid")
        library("vegan")
        library("MASS") #required for vegan
        ## functions
        source("PCAScores.R")
        #source("FAScores.R")
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
        if (type == "pca"){
                scores <- PCAScores(scent.abundance)
        } else if (type == "fa") {
                scores <- FAScores(scent.abundance)
        }
        
        NumFactors <- ncol(scores)
        
        ## subsetting all data objects
        Sub <- SubsetAll(subgroup,scent.abundance,scent.diversity,relatedness,
                         heterozygosity,scores,factors) # out:list 
        
        # from list back to original formats
        abund <- Sub[[1]]
        div <- Sub[[2]]
        relate <- Sub[[3]]
        het <- Sub[[4]]
        scores <- Sub[[5]]
        factors <- Sub[[6]]
        
        
        ## get pairwise pca differences and  relatedness as vectors or diff.pc matrix
        # diff.pc.mat.mums <- PCADiff(relate.mums,pca.scores.mums,num.pca=6, df=TRUE)
        diff.relate.vecs <- PCADiff(relate,scores, df=F)
        diff.relate.dfs <- PCADiff(relate,scores, df=T)
        ## diff.relate.vecs$genetic.distance <- 1 - diff.relate.vecs$genetic.distance
        relate.df <- diff.relate.vecs
        relate.list <- diff.relate.dfs
        
        ## Mantel tests (pearson) for all pcs and bray curtis matrix of overall scent abundance(PC0)
        results.relatedness <- ResultsRelatedness(relate,scores, abund, NumFactors)
        
        # create data.frame for heterozygosity plotting in ggplot //main data.frame for heterozygosity here
        het.df <- data.frame("het"= het[,1],row.names=rownames(het))
        het.df$NumComps <- div$S
        het.df[, 3:(2+NumFactors)] <- scores[, 1:NumFactors]
        
        # PC vs. Heterozygosity
        results.het <- resultsHet(het.df, NumFactors)
        
        # create pca.scores data.frame with factor beaches
        scores$beach <- factor(factors$Beach,labels=c("beach 1","beach 2")) # bloody hell, this is cool!!!
        
        ScentResults <- list(relate.df=relate.df,
                             results.relatedness=results.relatedness,
                             relate.list=relate.list,
                             relatedness=relate,
                             het.df=het.df,
                             results.het=results.het,
                             scores=scores,
                             factors=factors)
}