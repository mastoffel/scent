best_vec <- names(sort(table(all_best), decreasing = TRUE))
best_vec_nob <- names(sort(table(all_best_nob), decreasing = TRUE))


abund_df <- as.data.frame(abund)
gendiff <- 1-as.dist(relate)

stat <- vector()
sig <- vector()

for (i in 2:100) {
        bc_dist <- vegdist(abund_df[, best_vec_nob[1:i]], method = "bray")
        mod <- mantel(gendiff, bc_dist, na.rm = T)
        stat <- append(stat, mod$statistic)
        sig <- append(sig, mod$sig)
}


# get indices
elemnames <- which(colnames(scent_abundance) %in% best_vec[1:15])


library(ggplot2)