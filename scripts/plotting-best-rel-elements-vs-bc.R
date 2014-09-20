# best bootstrap relatedness mums
best <- c(36,52,86,88,96,103,110,203,206)
abund_best <- abund[, best]

library(vegan)

bc <- as.vector(vegdist(abund_best))
rel <- as.vector(as.dist(1-relate))

bc[bc==1] <- NA
rel[rel==1] <- NA

plot(rel, bc)
abline(lm((bc ~ rel),col="red"))

glm(bc ~ rel)

df <- as.data.frame(cbind(bc, rel))

ggplot(df, aes(x = rel, y = bc)) +
        geom_point() +
        geom_smooth(method = "lm") +
        theme_minimal(base_size = 20)

