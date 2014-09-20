# resampling test for finding out the substances associated with relatedness
# parallel computing using 40 cores


library(vegan)
library(stringr)
library(dplyr)
library(snow)
library(snowfall)
source("bio.env.R")

# create empty vector to count elements

scent_abundance <- as.data.frame(t(read.csv(paste(".\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))
relatedness <- read.csv(paste(".\\csv_files\\relatednessnew.csv", sep = ""),
                                row.names=1) 

# subset
age <- c(rep(TRUE,41),rep(FALSE,41))
abund <- filter(scent_abundance, age)
relate <-  relatedness[age, age]

rm(scent_abundance, relatedness)

# bootstrap
all_best <- vector()

# initialise cluster
sfInit(parallel=TRUE, cpus=2, type="SOCK")

# libraries
sfSource("bio.env.R")
sfLibrary(vegan)
sfLibrary(stringr)
sfLibrary(dplyr)

bootstrap <- function(iter_comp) { # main resampling function
        
        for (i in 1:500) {
                
                index_seals <- sort(sample(1:41, size = 20, replace = F))
                
                # subset relate and abund
                reltemp <- 1-as.dist(relate[index_seals, index_seals])
                abundtemp <- abund[index_seals, ]
                
                for (i in iter_comp) {
                        
                # sample 10 compounds
                        index_comps <- sort(sample(1:213, size = 10, replace = F))
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
                        # write(best, file = "best.txt", append = TRUE, sep = " ")
                }
        }
        
        return(all_best)
}

# export objects
sfExportAll(except = NULL, debug = FALSE)
sfClusterEval(ls())

# calculations

# create list of values for different cores
vals <- list()
for (i in 1:40) {
        vals[[i]] <- 1:500
}

best <- sfLapply(vals, bootstrap)

# stop cluster
sfStop()

# bring all results together and save as txt
results <- unlist(best)
write(results, file = "best.txt", append = TRUE, sep = " ")


