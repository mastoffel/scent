## Mum Pup seperate Analysis

# creating dist objects
scent.mum.dist <- as.dist(scent.mum)
scent.pup.dist <- as.dist(scent.pup)
rel.mum.dist <- as.dist(rel.mum)
rel.pup.dist <- as.dist(rel.pup)
scent.pa.mum.dist <- as.dist(scent.pa.mum)
scent.pa.pup.dist <- as.dist(scent.pa.pup)


# plots
nfrow = c(3,2)
plot(rel.mum.dist, scent.mum.dist)
plot(rel.pup.dist, scent.pup.dist)
plot(rel.mum.dist, scent.pa.mum.dist)
plot(rel.pup.dist, scent.pa.pup.dist)

library(ecodist)
nonranked.mantel  <- mantel(rel.mum.dist ~ scent.mum.dist) 
nonranked.mantel  <- mantel(rel.pup.dist ~ scent.pup.dist)