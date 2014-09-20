# building the bootstrapping test
library(vegan)
source("bio.env.R")
library(stringr)
#library(dplyr)
#source("bv.step.R")
#source("bioenvmod.R")
# create empty vector to count elements


all_best <- vector()

for (i in 1:10) {
        
        index_seals <- sort(sample(1:41, size = 20, replace = F))
        
        # subset relate and abund
        reltemp <- 1-as.dist(relate[index_seals, index_seals])
        abundtemp <- abund[index_seals, ]
        
        for (i in 1:10) {
                
        # sample 10 compounds
                index_comps <- sort(sample(1:183, size = 10, replace = F))
                abundtemp_sub <- abundtemp[, index_comps]
                
                # get vector with 0 for null-column and 1 for non-null column
                
                nullcomps <- apply(abundtemp_sub, 2, function(x) sum(x>0))
                
                abundtemp_sub <- subset(abundtemp_sub, 
                                        subset = c(rep(TRUE,nrow(abundtemp_sub))), 
                                        select = (nullcomps >= 2))
                
                # new iteration if too less substances left
                if (ncol(abundtemp_sub) <= 2) next
                
                results <- bio.env(reltemp, abundtemp_sub, 
                                var.dist.method = "bray", 
                                scale.fix = F, scale.var = F)
        
                mods <- results$best.model.vars
                best <- unlist(str_split(mods, ","))
                all_best <- append(all_best, best)
        }
}



                 