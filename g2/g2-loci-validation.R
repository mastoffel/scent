# script for resampling test for g2
source("g2.R")
num_iter <- 1000
genotypes <- read.table("raw_41loci_ordered.txt", row.names = 1)


results <- data.frame(matrix(nrow = num_iter, ncol = (ncol(genotypes)-1)/2 ))

for (i in seq(from = 2, to = (ncol(genotypes))/2)) {
        
        for (k in seq_along(1:num_iter)) {
                
                ## create genotype matrix with subset of loci (defined by loop)
                
                # create index from which to sample loci from genotypes
                ind <- seq(from = 1, to = ncol(genotypes), by = 2)
                
                # loci to take
                s <- sample(ind, i, replace = FALSE)
                s <- sort(s)
                
                # initialise subset data.frame
                genotypes.sub <- data.frame(matrix(nrow = nrow(genotypes),ncol = length(s)*2))
            
                
                # finalize subset data.frame with specific number of loci
                count <- 1
                for (m in s) {
                        genotypes.sub[, c(count,count+1)] <- as.matrix(subset(genotypes, select= c(m,m+1)) )
                        count <- count+2
                }
                
                # compute het data
                g_val <- g2(genotypes.sub)
                
                results[k, i] <- g_val
        }
}


library(ggplot2)
source("SumResults.R")
allresults <- SumResults(results)
ggplot(allresults, aes(x = locnum, y = cormean)) +
        geom_errorbar(aes(ymin = cormean-se, ymax = cormean+se),
                      width=0.8, alpha=0.7, colour="blue") +
        geom_point() +
        geom_line(colour = "blue", size=1) 

