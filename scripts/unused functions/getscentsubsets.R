## try: get scent abundances just for substantial substances
library(xlsx)
SA <- apply(scent.abundance, 2, function(x) length(x[x>0]) )

Scent <- as.data.frame(SA)

scent.abund <- scent.abundance[, which(Scent>8)]

write.xlsx(t(scent.abund),"scent.abund.xlsx")

scent.abund <- scent.abundance[, c(86,96,107,164,189)]

sa.dist <- vegdist(scent.abund)
sa.dist <- as.data.frame(as.matrix(sa.dist))


mantel(sa.dist,1-relatedness,na.rm=TRUE)
