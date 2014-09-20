# script for plot with increasing loci and correlation with relatedness for best substances
library(vegan)
library(Demerelate)
library(snow)
library(snowfall)
# load genotypes file with all 46 loci, which is ordered according to ID
genotypes <- read.table(".\\data\\raw_41loci_ordered.txt", na.strings = "NA", row.names = 1)
gen_mum <- genotypes[1:41, ]
# relatedness matrix (old: relatedness_41loci.csv)
scent_abundance <- as.data.frame(t(read.csv(".\\data\\scent abundances.csv", row.names=1)))

# best comps
comp_ind_m <- c(36,52,86,88,96,103,110,203,206)

# bc mat for mums
scent <- scent_abundance[1:41, ]
scent_bc <- vegdist(scent)
scent_mat <- as.matrix(scent_bc)
scent_mat[upper.tri(scent_mat)] <- 0

# initialise cluster
sfInit(parallel=TRUE, cpus=2, type="SOCK")
sfLibrary(vegan)
sfLibrary(Demerelate)

# get bc from best comps
numIter <- 1000

# every core gets a certain number of loci from 1:40

resample_loci <- function(loci) {
        
        out <- matrix(nrow = numIter, ncol = 1)
                
                for (k in seq_along(1:numIter)) {
                        
                        ## create genotype matrix with subset of loci (defined by loop)
                        
                        # create index from which to sample loci from genotypes
                        ind <- seq(from = 1, to = ncol(genotypes), by = 2)
                        
                        # loci to take
                        s <- sample(ind, loci, replace = FALSE)
                        s <- sort(s)
                        
                        # initialise subset data.frame
                        genotypes.sub <- data.frame(matrix(nrow = nrow(genotypes),ncol = length(s)*2))
                        
                        # finalize subset data.frame with specific number of loci
                        count <- 1
                        for (m in s) {
                                genotypes.sub[, c(count,count+1)] <- as.matrix(subset(genotypes, select= c(m,m+1)) )
                                count <- count+2
                        }
                        
                        genotypes_sub <- cbind(row.names(genotypes), as.factor(rep("AG", 82)), genotypes.sub)
                        
                        # compute relatedness data
                        rel <- emp.calc(genotypes_sub, value="rxy", ref.pop = genotypes_sub)
                        
                        # index het data for Mums
                        rel.mums <- rel[1:41, 1:41]
                        rel.mums[upper.tri(rel.mums)] <- 0
                        rel.m <- 1-rel.mums
                        mant <- mantel(rel.m, scent_mat, method = "spearman", na.rm = TRUE)
                        out[k, 1] <- mant$stat
                        
                
        }
        out
}

# export objects
sfExportAll(except = NULL, debug = FALSE)
sfClusterEval(ls())

locinum <- list()
for (i in 1:40) locinum[[i]] <- i

# let it run
results <- sfLapply(locinum, resample_loci)

# stop cluster
sfStop()

results_df <- data.frame(matrix(nrow = numIter, ncol = 40))

for (k in 1:40) {
        results_df[, k] <- results[[k]]
}

write.csv(results_df, file = "locirelate.csv")
