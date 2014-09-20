## factor analysis approach
library(psych)
## determining number of components
# parallel analysis --> 10
library(nFactors)
library(HDMD)
## from HDMD tutorial
VSS.scree(scent.abundance)
sc <- scree(scent.abundance)

# 
scent.fa <- factor.pa.ginv(scent.abundance, nfactors = 4, prerotate=T,
                           rotate = "promax", scores = T, m=4)
# loading raw data
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 

## diverse tests, vss ,map, scree, parallel
crits <- VSS(scent.abundance, n = 10, rotate = "promax", fm="pa")
fa.parallel(scent.abundance, fm = "promax")


nfactors(scent.abundance, n=10, rotate = "promax", fm = "pa")


scent.fa$scores

## nice! screeplot final solution
plot(scent.fa$values[1:10], type="lines")
