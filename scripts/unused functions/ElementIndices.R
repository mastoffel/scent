# computing indices of simper results for primer 

element.indices <- function(best.simper) {

RT <- read.csv(".\\csv_files\\RT_vector.csv",header=FALSE)
elements <- as.numeric(rownames(best.simper))

indices <- vector()
for (i in 1:length(elements)) {
  temp <- which(RT == elements[i])
  indices <- append(indices, temp)
}

ind.string <- as.character(indices)
element.indices <- paste(ind.string,collapse=",")
}
