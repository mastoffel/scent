bootstrap_loci <- function(genotypes, tocorrelateDF, tocorrelate, numIter) {
# "genotypes" is a matrix with ID´s in the first column and the genotypes in the others (one locus in two columns)
# "tocorrelateDF" is a data frame with one or more variables to correlate heterozygosity with
# tocorrelate is a number indicating the column to choose in the df
# numIter is the number of resampling per loci adding

# get mums from correlation DF
tocorrelateDF <- tocorrelateDF[1:41, ]

# initialize results df
results <- data.frame(matrix(nrow = numIter, ncol = (ncol(genotypes)-1)/2 ))
        
        for (i in seq(from = 1, to = (ncol(genotypes)-1)/2 )) {
                
                for (k in seq_along(1:numIter)) {
                        
                        ## create genotype matrix with subset of loci (defined by loop)
                        
                        # create index from which to sample loci from genotypes
                        ind <- seq(from = 2, to = ncol(genotypes), by = 2)
                        
                        # loci to take
                        s <- sample(ind, i, replace = FALSE)
                        s <- sort(s)
                        
                        # initialise subset data.frame
                        genotypes.sub <- data.frame(matrix(nrow = nrow(genotypes),ncol = length(s)*2 + 1))
                        genotypes.sub[, 1] <- genotypes[ ,1]
                        
                        # finalize subset data.frame with specific number of loci
                        count <- 2
                        for (m in s) {
                                genotypes.sub[, c(count,count+1)] <- as.matrix(subset(genotypes, select= c(m,m+1)) )
                                count <- count+2
                        }
                        
                        # compute het data
                        het <- mlhWS(genotypes.sub,"NA",4)
                        
                        # index het data for Mums
                        het.mums <- het[1:41, ]
                        
                        results[k, i] <- cor(as.vector(het.mums$SH),tocorrelateDF[, tocorrelate])
                }
        }
results
}