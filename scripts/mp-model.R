# multinomial logistic regression

library(mlogit)

allscores$family <- factors$Family

mldata <- mlogit.data(allscores, shape = "wide" )

library(MASS)
fit <- lda(family ~ F1 + F2 + F3 + F4, data = allscores, CV = TRUE)

scores <- allscores[-5]

score_diff <- abs(scores[1:41, ] - scores[42:82, ])




diff_sum <- colSums(score_diff)
diff_sd <- apply(score_diff, 2, sd)
diff_sum
diff_sd

bcsim <- as.matrix(vegdist(scent_abundance))
bc <- bcsim[1:41, 42:82]
bcmp <- diag(bc)

df <- cbind(bcmp, score_diff)

