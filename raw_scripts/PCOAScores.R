PCOAScores <- function(abundance.matrix){
        # computes pca.scores for first two PC´s from relative abundance matrix
     
        library(vegan) # provides factor analysis when D>>N
        PCOAscores <- as.data.frame(cmdscale(vegdist(abundance.matrix, method = "euclidian"), k = 2))   
        #scent.fa <- principal(abundance.matrix, nfactors = 3, rotate = "simplimax")
        
        # scree plot --> 2 PC seem to be optimal
        # plot(scent.pca,type="lines")
}