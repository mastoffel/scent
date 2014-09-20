library(data.table)

data <- fread("best_pups.txt")
all <- sort(table(data), decreasing = TRUE)
all <- as.data.frame(all)
names(all) <- c( "occurences")
write.csv(all, "bootstrap_pups.csv")
