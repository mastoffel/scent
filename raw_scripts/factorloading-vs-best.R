## 
# already standardized and transformed, transposed abundance matrix
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 

library(HDMD)
scent.fa <- factor.pa.ginv(scent.abundance[1:41], nfactors = 4, prerotate=T,
                           rotate = "promax", scores = T, m=4)

load    <- scent.fa$loadings
loaddf <- as.data.frame(load[, 1:4])
row.names(loaddf) <- names(scent.abundance)

# factor: 1 for best results compound (--> relatedness)
library(dplyr)
loaddf$best <- 0
loaddf$beach <- 0
loaddf[elemnames, 6] <- 1
loaddf$best[intersect.pa] <- 1
loaddfsort <- arrange(loaddf, F2)

# adding rank factor for number of occurences
occ <- apply(scent.abundance, 1, function(x) out <- sum(x > 0))
loaddf$occ <- occ

library(ggplot2)
ggplot(data=loaddfsort, aes(y = 1:213, x = F2)) +
        geom_point(aes(size = occ, color = as.factor(best))) +
        theme_minimal(base_size=16)



scent <- as.matrix((scent.abundance))
importance <-  loadF1 %*%  as.matrix(t(scent.abundance))

sorted.loadings <- load[order(load[, 1],decreasing=FALSE), 1] 

importance <- vector()
for (i in 1:213) {
        importance[i] <- sum(loadF1[i] * scent.abundance[, i])
}
names(importance) <- names(loadF1)

final <- sort(importance)

impdf <- as.data.frame(final)

impdf$simper <- 0

impdf$simper[c(2,3,6,15,20,29,33,50,60,86,96,104,110,114,118,126,130,133,134,141,157,186,206)] <- 1
impdf$simper <- as.factor(impdf$simper)

library(ggplot2)
ggplot(data=impdf, aes(y = 1:213, x = final, alpha = simper)) +
        geom_point(size = 4) +
        theme_minimal(base_size=16)