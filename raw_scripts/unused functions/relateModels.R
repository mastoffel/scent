## PERMANOVA models

F1 <- as.dist(relate.list[[1]])
F2 <- relate.list[[2]]
F3 <- relate.list[[3]]
F4 <- relate.list[[4]]
F5 <- relate.list[[5]]

relate <- as.dist(abs(relatedness))
adonis(relate ~ F1)

gendist <- relate.df$genetic.distance
F1v <- as.vector(as.matrix(F1))
F2v <- as.vector(as.matrix(F2))
F3v <- as.vector(as.matrix(F3))
F4v <- as.vector(as.matrix(F4))
F5v <- as.vector(as.matrix(F5))

F1v <- F1v[!is.na(F1v)]
F2v <- F2v[!is.na(F2v)]
F3v <- F3v[!is.na(F3v)]
F4v <- F4v[!is.na(F4v)]
F5v <- F5v[!is.na(F5v)]