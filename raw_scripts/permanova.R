## adonis (PERMANOVA) approach
## Analysis of variance using distance matrices - for partitioning distance matrices among sources of variation and 
## fitting linear models (e.g., factors, polynomial regression) to distance matrices; uses a permutation test with pseudo-F ratios.

library(vegan)
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 
scent.mum <- scent.abundance[1:41, ]
source("ScentResults.R")
res <- ScentResults("mums","fa",1)

# subset results list output
relate.df <- res[[1]]
results.relatedness <- res[[2]]
relate.list <- res[[3]]
relatedness <- res[[4]]
het.df <- res[[5]]
results.het <- res[[6]]
scores <- res[[7]]
factors <- res[[8]]

## all relevant similarity matrices
relatedness <- as.matrix(relatedness)
F1 <- as.matrix(relate.list[[1]])
F2 <- as.matrix(relate.list[[2]])
F3 <- as.matrix(relate.list[[3]])
F4 <- as.matrix(relate.list[[4]])

## make matrices symmetric
for (c in 1:(ncol(F1))) {
        for (r in 1:(nrow(F1))) {
                if (!is.na(F1[c,r])) {
                        relatedness[r,c] <- relatedness[c,r]
                        F1[r,c] <- F1[c,r]
                        F2[r,c] <- F2[c,r]
                        F3[r,c] <- F3[c,r]
                        F4[r,c] <- F4[c,r]
                }
        }
}

relatedness <- relatedness - min(relatedness[!(is.na(relatedness))]) ## add lowest value to make all values positive
## fill diagonals
diag(relatedness) <- 0
diag(F1) <- 0
diag(F2) <- 0
diag(F3) <- 0
diag(F4) <- 0

## create dist objects
rel
rel.dist <- as.dist(relatedness)
F1.dist <- as.dist(F1)
F2.dist <- as.dist(F2)
F3.dist <- as.dist(F3)
F4.dist <- as.dist(F4)


model <- adonis(rel.dist ~ scores$F1, method="euclidian" )

# Multiple Regression on distance Matrices
library(ecodist)
model <- MRM(rel.dist ~ F1.dist + F2.dist + F3.dist + F4.dist, nperm=100000)

# partial mantel
manteltest <- mantel(relatedness ~  F1 + F3+ F2 + F4, data = relate.df, nperm=10000)






