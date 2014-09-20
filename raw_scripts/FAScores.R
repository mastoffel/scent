FAScores <- function(abundance.matrix){
        # computes pca.scores for first two PC´s from relative abundance matrix
        library(psych)
        library(HDMD) # provides factor analysis when D>>N
        scent.fa <- factor.pa.ginv(abundance.matrix, nfactors = 4, prerotate=T,
                                   rotate = "promax", scores = T, m=3)
        #scent.fa <- principal(abundance.matrix, nfactors = 3, rotate = "simplimax")
        FAScores <- as.data.frame(scent.fa$scores)
        # scree plot --> 2 PC seem to be optimal
        # plot(scent.pca,type="lines")
}