## factor analysis approach
library(psych)
## determining number of components
# parallel analysis --> 10
library(nFactors)
ev <- eigen(cor(scent.abundance))
ap <- parallel(subject=nrow(scent.abundance), var=ncol(scent.abundance) ,rep=100, cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

## diverse tests, vss ,map
vss(scent.abundance, n = 15, rotate = "varimax", fm="pc")


# tried the same on the prcomp way, didn´t work
scent.pca <- prcomp(scent.abundance)
loads <- scent.pca$rotation
newloads <- varimax(loads)$loadings[, 1:5]
newscores <- as.matrix(scent.abundance) %*% newloads


#fa.parallel(scent.abundance)
scores <- as.data.frame(scent.pca$scores)

pca.scores[, 1:5] <- scores[, 1:5]


cor(scent.fa$scores, scent.pca$scores)

all <- vector()
for (i in 1:10) {
        scent.fa <- factor.pa.ginv(scent.abundance, nfactors = i, prerotate=FALSE, rotate = "promax", scores = T)
        all[i] <- scent.fa$fit
}