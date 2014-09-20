## Models
# Scent Mothership script has to be run

## creating stuff for mum-pup model
require(vegan)
source("ParsModel.R")
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 
source("ScentResults.R")

res <- ScentResults("all","fa",1)

# subset results list output
relate.df <- res[[1]]
results.relatedness <- res[[2]]
relate.list <- res[[3]]
relatedness <- res[[4]]
het.df <- res[[5]]
results.het <- res[[6]]
scores <- res[[7]]
factors <- res[[8]]

#### mum-pup model ####

# get bray curtis dissimilarities from R 
BC <- as.data.frame(as.matrix(vegdist(scent.abundance, method = "bray",na.rm=TRUE)))

# get Mum-Pup pair values for Scent Dissimilarity
BCrel <- BC[42:82,1:41]
MP <- diag(as.matrix(BCrel))

# load scores with Scent Mothership Script

FaDiff <- matrix(nrow=41,ncol=4)
FaDiff[1:41, 1:4] <- abs(as.matrix(scores[1:41, 1:4]) - as.matrix(scores[42:82, 1:4]))
FaDiff <- as.data.frame(FaDiff)
names(FaDiff) <- c("F1","F2","F3","F4")

#create df with all
MPdf <- cbind(MP,FaDiff)
MumPupModel <- ParsModel(MPdf, test = "Chisq")



#### beach model ####

beachdata <- het.df[, -2]
beachdata$het <- as.factor(factors$Beach)
names(beachdata)[1] <- "beach"

BeachModel <- ParsModel(beachdata, family="binomial")



