## factor analysis approach
library(psych)

## determining number of components

# parallel analysis --> 10

library(nFactors)
library(HDMD)
## from HDMD tutorial
VSS.scree(scent_abundance)
sc <- scree(scent_abundance)
# 
scent.fa <- factor.pa.ginv(scent_abundance, nfactors = 4, prerotate=T,
                           rotate = "promax", scores = T, m=4)

## diverse tests, vss ,map, scree, parallel
crits <- VSS(scent_abundance, n = 10, rotate = "promax", fm="pa")
fa.parallel(scent_abundance, fm = "promax")


nfactors(scent_abundance, n=10, rotate = "promax", fm = "pa")


scent.fa$scores

## nice! screeplot final solution
plot(scent.fa$values[1:10], type="lines")
