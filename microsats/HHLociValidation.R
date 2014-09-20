HHLociValidation <- function(genotypes, numIter) {
        # "genotypes" is a matrix with ID´s in the first column and the genotypes in the others (one loci in two columns)
        # "tocorrelateDF" is a data frame with one or more variables to correlate heterozygosity with
        # tocorrelate is a number indicating the column to choose in the df
        # numIter is the number of resampling per loci adding
        
        genotypes <- genotypes[, -1]
        # initialize results df
        results <- data.frame(matrix(nrow = numIter, ncol = (ncol(genotypes))/2 ))
        
        for (i in seq(from = 2, to = (ncol(genotypes))/2 )) {
                
                for (k in seq_along(1:numIter)) {
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
                        
                        # compute hethet correlation
                        hethet <- hh(genotypes.sub,n = 100,method = "sh")
                        hetmean <- mean(hethet)
                        results[k, i] <- hetmean
                }
        }
        results
}