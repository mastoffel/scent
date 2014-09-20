# coordinates analysis: TO DO: relatedness!

library(ggplot2)
library(vegan)
# creating clean xy data.frame (just Beach 1)

# load coord and het data
scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "MSc.Behaviour\\Research\\Seal Scent\\R code\\",
                                                  "data\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))
coord <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                        "MSc.Behaviour\\Research\\Seal Scent\\R code\\",
                        "data\\csv_files\\coordinates.csv",
                        sep = ""),row.names=1) 

factors <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                          "MSc.Behaviour\\Research\\Seal Scent\\",
                          "R code\\data\\csv_files\\",
                          "factors.csv", sep = ""),
                    row.names=1) 

heterozygosity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                 "MSc.Behaviour\\Research\\Seal Scent\\",
                                 "R code\\data\\csv_files\\",
                                 "heterozygosity_41loci.csv", sep = ""),
                           row.names=1) 

# create subset with relevant individuals (40)
sub <- row.names(coord) %in% row.names(heterozygosity) 
coordinates <- subset(coord, sub)
coord <- subset(coordinates, factors$Beach==1)
coord["M48",1] <- "0618"
coord["P48",1] <- "0618"
coord$XY <- as.character(coord$XY)

# split coord
X <- apply(coord, 1, function(x) substring(x,1,2))
Y <- apply(coord, 1, function(x) substring(x,3,4))
coord$X <- X
coord$Y <- Y
coord$X <- as.numeric(coord$X)
coord$Y <- as.numeric(coord$Y)

# plot Seal distribution
ggplot(coord, aes(X,Y)) +
        geom_point(colour="blue", size = 4, alpha = 0.5) +
        theme_bw(base_size=18)

# create dist mat
dist.mat <- as.matrix(dist(coord[, 2:3]))
distances <- as.data.frame(dist.mat)
distances <- as.matrix(distances)
distances[upper.tri(distances, diag=TRUE)] <- NA
distances <- as.data.frame(distances)
#write.xlsx(distances,"distances.xlsx")

# getting bray curtis matrix from beach one
abund <- scent_abundance[factors$Beach==1, ]

# substances explaining mum pup similarity
simper_mp_ind <- c(58, 68,  86,  90, 106, 107, 164)
abund_sub <- abund[, simper_mp_ind]

bcSim <- as.matrix(vegdist(abund_sub, method="bray"))        
bcSim[upper.tri(bcSim, diag=TRUE)] <- NA   
bcSim <- as.data.frame(bcSim)


# seperate analysis for mums, pups
bcSimMum <- bcSim[1:20,1:20]
distMum <- distances[1:20,1:20]

bcSimPup <- bcSim[21:40,21:40]
distPup <- distances[21:40,21:40]

# mantel test
# mums, pups
mantmum <- mantel(bcSimMum,distMum, method = "spearman")
mantpup <- mantel(bcSimPup,distPup, method = "spearman")

# all
mantall <- mantel(bcSim, distances, method = "spearman")

# plotting
bcSim.vec <- as.vector(as.matrix(bcSim))
distances.vec <- as.vector(as.matrix(distances))
plot(bcSim.vec,distances.vec)
fit <- lm(distances.vec~bcSim.vec)
abline(fit)
