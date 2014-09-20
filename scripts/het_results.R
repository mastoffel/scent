het_results <- function(het_df, num_factors) {
        # all is a data.frame containing het values and pc scores
        # FA 0 is the label used for comparison between het and all compounds 
        
        results <- data.frame("FA" = 0:num_factors)
        
        for (i in 2:(2 + num_factors)) {
                
                model <- glm(het_df$het ~ het_df[, i])
                model.ano <- anova(model,test="F")
                c <- cor(het_df$het, het_df[, i])
                r.squared <- c^2
                results$r[i-1] <- c
                results$rsquared[i-1] <- r.squared
                results$f[i-1] <- model.ano[2, 5]
                results$p[i-1] <- model.ano[2, 6]
                
        }
        
        results$fdr_corrected <- p.adjust(results$p, method="fdr")
        return(results)
}