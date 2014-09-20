## Simper analysis

# packages
library(vegan)

# load

# already standardized and transformed, transposed abundance matrix
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 

# 6 different diversity indices (from primer)
scent.diversity <- read.csv(".\\csv_files\\scent diversity.csv",row.names=1) 

# relatedness matrix
relatedness <- read.csv(".\\csv_files\\relatedness_41loci.csv",row.names=1)

# heterozygosity vector SH
heterozygosity <- read.csv(".\\csv_files\\heterozygosity_41loci.csv", row.names=1) 

# beach and family vector
factors <- read.csv(".\\csv_files\\factors.csv",row.names=1) 

# get beach factor for simper
beach <- as.factor(factors$Beach)

## average dissimilarity
simp <- simper(scent.abundance, beach)
best30 <- summary(simp)[[1]]


source("element_indices.R")
## element indices
elemind <- element_indices(best30)
# element names
elemnames <- rownames(best30)[1:30]
# write df
#scent.abund.nobeach <- scent.abundance[, -which(names(scent.abundance) %in% elemnames)]
#write.csv(t(scent.abund.nobeach), "scent abundance nobeach.csv")

library(HDMD)
scent.fa <- factor.pa.ginv(scent.abundance, nfactors = 4, prerotate=T,
                           rotate = "promax", scores = T, m=4)
load    <- scent.fa$loadings
sorted.loadings <- load[order(load[, 1],decreasing=FALSE), 1] # change both numbers for PC change
#sorted.loadings.1 <- load[order(load[, 1]), 1]
loaddf <- as.data.frame(sorted.loadings)
loaddf$simper <- 0
loaddf$simper[rownames(loaddf) %in% elemnames] <- 1
loaddf$simper[c(5,6,86,96,110,206)] <- 1
loaddf$simper <- as.factor(loaddf$simper)

## check loadings from best analysis substances

loadbest <- as.data.frame(load[, 1:4])
loadbest$best <- 0
loadbest$best[c(5,6,86,96,110,206)] <- 1


ggplot(data=loaddf, aes(y = 1:213, x = sorted.loadings, alpha = simper)) +
        geom_point(size = 4) +
        theme_minimal(base_size=16)
