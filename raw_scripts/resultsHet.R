resultsHet <- function(All, NumFactors) {
# all is a data.frame containing het values and pc scores
# PC 0 is the label used for comparison between het and all compounds 
results <- data.frame(data.frame("F"=0:NumFactors))

for (i in 2:(2+NumFactors)) {
        model <- glm(All$het ~ All[, i])
        model.ano <- anova(model,test="F")
        c <- cor(All$het,All[, i])
        r.squared <- c^2
        results$r[i-1] <- c
        results$rsquared[i-1] <- r.squared
        results$f[i-1] <- model.ano[2, 5]
        results$p[i-1] <- model.ano[2, 6]
}

results$fdr_corrected <- p.adjust(results$p, method="fdr")
return(results)
}