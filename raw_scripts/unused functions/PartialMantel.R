## Partial Mantel

## get df
source("ScentResults.R")
res <- ScentResults("pups","fa",1)

# subset results list output
relate.df <- res[[1]]

## vegan and ecodist´s mantel interfere, thus attach library after getting df´s
library(ecodist)

## y ~ x + z1 + z2 + z3 will do a partial Mantel test of the relationship 
## between x and y given z1, z2, z3. All variables can be either a distance matrix of 
## class dist or vectors of dissimilarities.

manteltest <- mantel(relatedness ~  F1 + F3+ F2 + F4, data = relate.df, nperm=10000)


