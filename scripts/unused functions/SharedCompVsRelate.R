## idea: number of shared substances vs. relatedness

scent.mum <- scent.abundance[1:41, ]

current <- 1

## count function
CompShare <- function(x) {
        if (x[1]>0 & x[2]>0) {
                CompShare <- 1
        } else {
                CompShare <- 0
        }
}

numInd <- nrow(scent.mum)
## define df for shared comps
ShareSim <- relatedness 

# count and fill
for (i in c(1:numInd)) {
        for (k in c(1:numInd)[-current]) {
                name1 <- rownames(scent.mum)[i]
                name2 <- rownames(scent.mum)[k]
                
                shared <- apply(scent.mum[c(i,k), ],2,CompShare) ##get rows
                shared.sum <- sum(shared)
                
                if (!is.na(relatedness[name1,name2])) {
                        ShareSim[name1,name2] <- shared.sum
                } else {
                        ShareSim[name2, name1] <- shared.sum
                }
                
        }
        current <- current +1 
}


SharedVec <- as.vector(as.matrix(ShareSim))
RelVec <- 1-relate.df$genetic.distance
Shared <- SharedVec[!is.na(SharedVec)]

plot(RelVec,Shared)
mantel(ShareSim,relatedness)
