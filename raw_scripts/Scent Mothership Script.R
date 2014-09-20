## Main Scripts
# 1 = special study beach, 2 = freshwater beach
# loading raw data
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 

source("ScentResults.R")
res <- ScentResults("mums","fa", beach=0)

# subset results list output
relate.df <- res[[1]]
results.relatedness <- res[[2]]
relate.list <- res[[3]]
relatedness <- res[[4]]
het.df <- res[[5]]
results.het <- res[[6]]
scores <- res[[7]]
factors <- res[[8]]

# get factor solution
scent.fa <- factor.pa.ginv(scent.abundance, nfactors = 4, prerotate=T,
                           rotate = "promax", scores = T, m=4)


load    <- scent.fa$loadings
sorted.loadings <- load[order(load[, 1],decreasing=TRUE), 1] # change both numbers for PC change
#sorted.loadings.1 <- load[order(load[, 1]), 1]
loaddf <- as.data.frame(sorted.loadings)
myTitle <- "Loadings Plot for PC1" 
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")
