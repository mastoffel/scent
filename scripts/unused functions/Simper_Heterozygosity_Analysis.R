## Most Explaining Substances for Heterozygosity Analysis


## Highest loadings on PC2 ##

scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 

# very nice plot: Variable loadings vs. PC
scent.pca <- prcomp(scent.abundance) 
load <- scent.pca$rotation

# getting absolute value of all loadings
load[, 2] <- abs(load[, 2])
sorted.loadings.PC2 <- load[order(load[, 2], decreasing=TRUE), 2] # change both numbers for PC change
top10 <- sorted.loadings.PC2[1:10]

#sorted.loadings.1 <- load[order(load[, 1]), 1]
myTitle <- "Loadings Plot for PC2" 
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")



## Simper ##
