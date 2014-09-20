library(vegan)

testMDS <- metaMDS(scent.abundance, k = 4, plot=TRUE, autotransform=FALSE, wascores=TRUE)
mdsscores <- scores(testMDS)
mumscores <- mdsscores[1:41, ]
scores[, 1:4] <- mumscores[, 1:4]
