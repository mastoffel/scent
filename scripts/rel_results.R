rel_results <- function(relate, scores, abund, num_factors) {
# inputs: relatedness matrix, table of pc.scores. Output: data.frame with all
# necessary outcomes for the possible Mantels
# caution: genetic DISTANCE is used instad of relatedness
        
        # create data frame
        results <- data.frame("FA" = 0:num_factors)
        
        # gen.dist <- relate ## 1- ?
        
        rel <- 1 - relate
        
        # create bray curtis matrix from scent abundance and put na in upper tri
        bcSim <- as.matrix(vegdist(abund, method="bray"))        
        bcSim[upper.tri(bcSim, diag=TRUE)] <- NA   
        
        
        for (i in 0:num_factors) {
                if (i == 0) {
                        model<- mantel(rel,bcSim,method="pearson",na.rm=TRUE,permutation=9999)
                        results$mantelR[i+1] <- model$statistic
                        results$p[i+1] <- model$signif
                }
                if (i>0) {
                        facdiff_list <- get_pairdiff(rel, scores, df=TRUE)
                        model <- mantel(rel, facdiff_list[[i]],method="pearson",na.rm=TRUE,permutation=9999)
                        results$mantelR[i+1] <- model$statistic
                        results$p[i+1] <- model$signif
                }
        }
        
        results$fdr_corrected <- p.adjust(results$p,method="fdr")
        results$bonferroni <- p.adjust(results$p,method="bonferroni")
        results$holm <- p.adjust(results$p,method="bonferroni")
        return(results)
}