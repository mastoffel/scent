best.simper.list <- function(simperfile) {



# reading in simper results as csv (first row has to be "Group X")
# output is list (length: number of groups) with simper results (first column: component, second column: explained variance)

simper <- read.csv(simperfile, header = FALSE)
simper2 <- simper[,-c(2,3,5)]


# extracting group names
groups <- vector()
count <- 1
for (i in 1:nrow(simper2)) {
  if (pmatch("Group",simper2$V1[i], nomatch = 0) == 1) {
    groups[count] <- as.character(simper2$V1[i])
    count <- count + 1
  }
}
rm(count,i)

# extract simper substances of each group and put them into list of indix vectors

library(stringr)

extract.list <- list() #initialize
count <- 1
mat.start <- 5 # starts at 5 because its easier to extract then..
mat.stop <- 0
for (i in (5:nrow(simper2))) {
  if (str_detect((simper2[i,1]), "Species") == TRUE) {
    mat.start <- i+1
  }
  if (str_detect((simper2[i,1]), "Group") == TRUE) {
    mat.stop <- i-2
    extract.list[[count]] <- c(mat.start:mat.stop)
    count <- count + 1
  } else if (i <- nrow(simper2)) {
    mat.stop <- i
    extract.list[[count]] <- c(mat.start:mat.stop)
  }
}

# getting data frame of simper values in columns // creating list of data frames, first row component, second row
#simper.data <- data.frame(t(rep(NA,length(groups))))
simper.data <- list()

for (i in seq_along(extract.list)) {
  simper.data[[i]] <- simper2[extract.list[[i]],]
  names(simper.data[[i]]) <- c("comp","var")
}











