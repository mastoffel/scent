## factor structure investigation
## idea: plot --> deviance explained by each model for relatedness and heterozygosity
## taking the model with best deviance explanation
source("ParsModel.R")
source("ScentResults.R")



## 1 factor
res <- ScentResults("mums","fa",1)

relate.df <- res[[1]]
results.relatedness <- res[[2]]
relate.list <- res[[3]]
relatedness <- res[[4]]
het.df <- res[[5]]
results.het <- res[[6]]
scores <- res[[7]]