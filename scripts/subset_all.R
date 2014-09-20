subset_all <- function(subset, beach, scent_abundance,scent_diversity,
                       relatedness, heterozygosity, allscores, factors) {
# subsets all input data objects. the first input should be "All", "Mums", "Pups" for specific subsetting
        require(dplyr)
        require(magrittr)
        
        # create subset factor
        if (subset == "all") age <- rep(TRUE,82)
        if (subset == "mums") age <- c(rep(TRUE,41),rep(FALSE,41))
        if (subset == "pups") age <- c(rep(FALSE,41),rep(TRUE,41))
        
        # create non-subset if beach = 0
        if (!beach == (1 | 2)) beach <- factors$Beach
        
        # subset all data by age and beach
        scent_abund <- filter(scent_abundance, age & (factors$Beach == beach))
        scent_div <- filter(scent_diversity, age & factors$Beach == beach)
        rel <-  relatedness[age & factors$Beach == beach, age & factors$Beach == beach]
        het <- filter(heterozygosity, age & factors$Beach == beach)
        allscores <- filter(allscores, age & factors$Beach == beach)
        facs <- filter(factors, age & factors$Beach == beach)
        
        
        SubsetAll <- list(scent.abund = scent_abund,
                          scent.div = scent_div,
                          rel = rel,
                          heterozygosity=het,
                          allscores=allscores,
                          factors=facs)
        return(SubsetAll)
}
