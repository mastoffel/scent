#load data for mums

scent <- read.csv(".\\csv_files\\ScentSim_Mums.csv",row.names=1)
rel <- read.csv(".\\csv_files\\Relatedness_Mums.csv",row.names=1)

# from matrix to vector
rel.vec <- as.vector(as.matrix(rel))
rel.vec <- rel.vec[!is.na(rel.vec)]

scent.vec <- as.vector(as.matrix(scent))
scent.vec <- scent.vec[!is.na(scent.vec)]

plot(scent.vec ~ rel.vec)
lm <- lm(scent.vec ~ rel.vec)
abline(lm,col="red")

