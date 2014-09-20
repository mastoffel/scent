get_scores <- function(abundances, method = "fa", num_dim = 4, 
                       rotation = "varimax") {
        library(HDMD)
        library(vegan)
        library(magrittr)
        
        if  (method == "pca") {
                scent_ord <- prcomp(abundances)
                allscores <- as.data.frame(scores(scent_ord, choices = 1:num_dim))
                
        } else if (method == "fa") {
                scent_ord <- factor.pa.ginv(abundances, nfactors = num_dim, 
                                            prerotate = T,rotate = rotation, 
                                            scores = T, m = 3)
                allscores <- as.data.frame(scent_ord$scores)
                
        } else if (method == "pcoa") {
                allscores <- vegdist(abundances, method ="euclidian") %>%
                cmdscale(k = num_dim) %>%
                as.data.frame()
                
        }
        allscores
}





