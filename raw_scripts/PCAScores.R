PCAScores <- function(abundance.matrix){
        # computes pca.scores for first two PC´s from relative abundance matrix
        scent.pca <- prcomp(abundance.matrix)
        pc1.score <- scores(scent.pca, choices=1)
        pc2.score <- scores(scent.pca, choices=2)
        pc3.score <- scores(scent.pca, choices=3)
        pc4.score <- scores(scent.pca, choices=4)
        pc5.score <- scores(scent.pca, choices=5)
        pc6.score <- scores(scent.pca, choices=6)
        pc7.score <- scores(scent.pca, choices=7)
        pc8.score <- scores(scent.pca, choices=8)
        pc9.score <- scores(scent.pca, choices=9)
        pc10.score <- scores(scent.pca, choices=10)
        pca.scores <- data.frame("pc1"=pc1.score,"pc2"=pc2.score,"pc3"=pc3.score,
                                 "pc4"=pc4.score,"pc5"=pc5.score,"pc6"=pc6.score,
                                 "pc7"=pc7.score,"pc8"=pc8.score,"pc9"=pc9.score,
                                 "pc10"=pc10.score)
        # scree plot --> 2 PC seem to be optimal
        # plot(scent.pca,type="lines")
}