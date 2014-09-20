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


## element indices

elemind <- element.indices()
