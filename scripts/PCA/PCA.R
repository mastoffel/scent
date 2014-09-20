# load mum relatedness matrix
mums.relate <- read.csv(".\\csv_files\\Relatedness_Mums.csv",row.names=1)

# load PCS Mums
mums.PCA <- read.csv(".\\csv_files\\PCA_Mums.csv",row.names=1)

# copy similarity matrix and clear
mums.PCA.mat <- mums.relate
mums.PCA.mat[,] <- NA

# construct similarity matrix out of pairwise differences in principal components

names <- rownames(mums.relate)

for (i in names) {
        for (k in names) {
                if (!(is.na(mums.relate[i,k]))) {
                        diff.PCA <- abs(mums.PCA[i,1]-mums.PCA[k,1])
                        mums.PCA.mat[i,k] <- diff.PCA
                }
        }
}

#turn into vector

mums.pca <- as.vector(as.matrix(mums.PCA.mat))
mums.pca <- mums.pca[!is.na(mums.pca)]
mums.relate <- as.vector(as.matrix(mums.relate))
mums.relate <- mums.relate[!is.na(mums.relate)]

lm <- lm(mums.pca ~ mums.relate)
cor <- cor(mums.pca,mums.relate)
# plot
library(ggplot2)
library(ecodist)
S <- qplot(mums.relate,mums.pca) +
        geom_point(colour = "black", size = 2) +
        geom_smooth(method="lm",size = 1 ,col="black",alpha=0.5) +
        geom_text(size=8,aes(0.3,6.5, label="r = 0.07, p = 0.03")) +
        theme(axis.title.x = element_text(vjust=0.1,size = 16),
              axis.title.y = element_text(vjust=0.1,size = 16))

# mantel(mums.relate~mums.pca)
