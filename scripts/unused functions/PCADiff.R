PCADiff <- function(relate.df,pca.df,df=F) {
        # creates data.frame with  1) pairwise mums relatedness 2) pairwise differences in principal components
        # input should be: relatedness data frame (lower triangular), data frame with factor/pc scores in columns
        # if df=TRUE, pca.diff will return a list of pca dataframes (for each component/factor) with pairwise pc-differences. 
   

# copy similarity matrix and clear
PCA.mat <- relate.df
PCA.mat[,] <- NA

# get relatedness vector
#relate <- as.vector(as.matrix(relate.df))
#pairnames <- names(unlist(relatedness))
#rownames(relate) <- pairnames

##get vector of pair-rownames
allnames <- vector()
for (i in 1:ncol(relate.df)) {
        for (k in 1:nrow(relate.df)) {
                nametemp <- paste(names(relate.df)[i], row.names(relate.df)[k], sep = "")
                allnames <- append(allnames, nametemp)
        }
}


relate <- unlist(relate.df)
names(relate) <- allnames
relate <- relate[!is.na(relate)]
pairnamessub <- names(relate)
FacDiffVecs <- data.frame("relatedness"= relate)

# construct similarity matrix out of pairwise differences in principal components

names <- rownames(relate.df)
FacDiffMats <- list()

for (z in 1:ncol(pca.df)) {
        for (i in names) {
                for (k in names) {
                        if (!(is.na(relate.df[i,k]))) {
                                diff.PCA <- abs(pca.df[i,z] - pca.df[k,z])
                                PCA.mat[i,k] <- diff.PCA
                        }
                }
        }
        
        # create list of data frames, containing difference matrices per Factor
        FacDiffMats <- c(FacDiffMats, list(PCA.mat))
        
        # turn into vector
        pca <- as.vector(as.matrix(PCA.mat))
        pca <- pca[!is.na(pca)]
        
        FacDiffVecs <- cbind(FacDiffVecs, pca)
}

## check argument for what to return
        if (df == T) {
                names(FacDiffMats) <- names(pca.df)
                return(FacDiffMats)
        } else if (df == F) {
                names(FacDiffVecs) <- c("relatedness",names(pca.df))
                row.names(FacDiffVecs) <- pairnamessub
                return(FacDiffVecs)
        }


}