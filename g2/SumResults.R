SumResults <- function(results) {
# sums results matrix into a matrix with means, sd and se

mean.cor <- apply(results,2,mean, na.rm=T)
sd.cor <- apply(results,2,sd, na.rm=T)
allresults <- data.frame(locnum = 1:ncol(results), cormean = mean.cor, corsd = sd.cor)
allresults$se <- allresults$corsd/(sqrt(length(allresults$corsd)))

allresults
}