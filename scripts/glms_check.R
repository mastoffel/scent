# trying to redo the models

bcSim <- as.data.frame(as.matrix(vegdist(scent_abundance[factors$Beach == 1, ])))
bcSim_pairs <- as.matrix(bcSim[1:20, 21:40])
scores_b1 <- scores[factors$Beach == 1, ]
# mum pup similarity
mpsim <- diag(bcSim_pairs)

names <- c("F1", "F2", "F3", "F4")
for (i in 1:4){
        assign(names[i], abs(scores_b1[1:20, i] - scores_b1[21:40, i]))
}

model <- glm(mpsim ~ F1 + F2 + F3 + F4)



bcSim <- as.data.frame(as.matrix(vegdist(scent_abundance[factors$Beach == 2, ])))
bcSim_pairs <- as.matrix(bcSim[1:21, 22:42])
scores_b2 <- scores[factors$Beach == 2, ]

# mum pup similarity
mpsim <- diag(bcSim_pairs)

names <- c("F1", "F2", "F3", "F4")
for (i in 1:4){
        assign(names[i], abs(scores_b2[1:21, i] - scores_b2[22:42, i]))
}

model <- glm(mpsim ~ F1 + F2 + F3 + F4)

