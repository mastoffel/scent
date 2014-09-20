ResultsRelatedness <- function(relate,pca.scores,abund, NumFactors) {
# inputs: relatedness matrix, table of pc.scores. Output: data.frame with all necessary outcomes for the possible Mantels
# caution: genetic DISTANCE is used instad of relatedness

        
results <- data.frame("PC"=0:NumFactors)
# gen.dist <- relate ## 1- ?
rel <- 1-relate
# create bray curtis matrix from scent abundance
bcSim <- as.matrix(vegdist(abund, method="bray"))        
bcSim[upper.tri(bcSim, diag=TRUE)] <- NA   

        for (i in 0:NumFactors) {
                if (i == 0) {
                        model<- mantel(rel,bcSim,method="pearson",na.rm=TRUE,permutation=9999)
                        results$mantelR[i+1] <- model$statistic
                        results$p[i+1] <- model$signif
                }
                if (i>0) {
                pc.mat <- PCADiff(rel,pca.scores, df=TRUE)
                model <- mantel(rel,pc.mat[[i]],method="pearson",na.rm=TRUE,permutation=9999)
                results$mantelR[i+1] <- model$statistic
                results$p[i+1] <- model$signif
                }
        }

results$fdr_corrected <- p.adjust(results$p,method="fdr")
results$bonferroni <- p.adjust(results$p,method="bonferroni")
results$holm <- p.adjust(results$p,method="bonferroni")
return(results)
}