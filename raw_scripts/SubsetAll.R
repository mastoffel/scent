SubsetAll <- function(subset,scent.abundance,scent.diversity,relatedness,heterozygosity,pca.scores, factors) {
# subsets all input data objects. the first input should be "All", "Mums", "Pups" for specific subsetting
        
        #create subset factor
        if (subset == "all") factor <- rep(TRUE,82)
        if (subset == "mums") factor <- c(rep(TRUE,41),rep(FALSE,41))
        if (subset == "pups") factor <- c(rep(FALSE,41),rep(TRUE,41))
        
        # subset all data by logical value (factor)
        scent.abundance <- scent.abundance[factor, ]
        scent.diversity <- scent.diversity[factor, ]
        relatedness <- relatedness[factor,factor]
        heterozygosity <- subset(heterozygosity, factor)
        pca.scores <- subset(pca.scores,factor)
        factors <- factors[factor, ]
        SubsetAll <- list(scent.abundance=scent.abundance,
                          scent.diversity=scent.diversity,
                          relatedness=relatedness,
                          heterozygosity=heterozygosity,
                          pca.scores=pca.scores,
                          factors=factors)
        return(SubsetAll)
}
