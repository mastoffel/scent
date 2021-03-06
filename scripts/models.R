## Models

require(vegan)
source("ParsModel.R")

scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "projects\\sealscent\\data_files\\",
                                                  "Rdata\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))


# Scent Mothership script has to be run
source("all_scent.R")
res <- all_scent("all", beach = 0)

# subset results list output
relate.df <- res[[1]]
results.relatedness <- res[[2]]
relate.list <- res[[3]]
relatedness <- res[[4]]
het.df <- res[[5]]
results.het <- res[[6]]
scores <- res[[7]]
factors <- res[[8]]


#### Heterozygosity model ####

hetmodel <- ParsModel(het.df[, -2])



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

mod <- glm(beachdata$beach ~. , data = beachdata[, 4:5], family = binomial)

# relatedness models wioth ecodist mantel
source("all_scent.R")
res <- all_scent("mums", beach=0)
relate.df.mums <- res[[1]]

source("all_scent.R")
res <- all_scent("pups", beach=0)
relate.df.pups <- res[[1]]

res <- ScentResults("mums","fa", beach=1)
relate.df.mums.1 <- res[[1]]

res <- ScentResults("mums","fa", beach=2)
relate.df.mums.2 <- res[[1]]

library(ecodist)

mantel(relatedness ~ F1 + F2 + F3 + F4, data = relate.df.mums, mrank = T)
mantel(relatedness ~ F2 + F1 + F3 + F4, data = relate.df.mums, mrank = T)
mantel(relatedness ~ F3 + F2 + F3 + F4, data = relate.df.mums, mrank = T)
mantel(relatedness ~ F4 + F1 + F3 + F4, data = relate.df.mums, mrank = T)

mantel(relatedness ~ F1 + F2 + F3 + F4, data = relate.df.pups, mrank = T)
mantel(relatedness ~ F2 + F1 + F3 + F4, data = relate.df.pups, mrank = T)
mantel(relatedness ~ F3 + F2 + F3 + F4, data = relate.df.pups, mrank = T)
mantel(relatedness ~ F4 + F1 + F3 + F4, data = relate.df.pups, mrank = T)

detach("package:ecodist", unload=TRUE)
# taken out
mantel_mums_B1 <- mantel(relatedness ~ F1 + F2 + F3 + F4, data = relate.df.mums.1)
mantel_mums_B2 <- mantel(relatedness ~ F1 + F2 + F3 + F4, data = relate.df.mums.2)
