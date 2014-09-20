## analyses on Simper
# input: list of data.frames with 2 columns each(element + explained variance)
# output: data frame with elements from simper as rows, groups as columns and explained variance as values

best.simper <- function(simper.data, number.extract) {

extract <- function(x) return(x[1:number.extract,1])

# extract alsways first two elements
to.extract <- sapply(simper.data,extract)

to.extract <- as.numeric(as.vector(to.extract))

# delete duplicates and get a vector of unique elements
uni <- unique(to.extract)
uni <- sort(uni)


# make data.frame with unique as column names, individuals as rownames and explained variance in matrix
# create empty frame
simper.frame <- as.data.frame(matrix(NA, ncol = length(groups), nrow = length(uni)))
names(simper.frame) <- groups
rownames(simper.frame) <- as.character(uni)

## find explained variances to fill matrix

for (i in seq_along(simper.data)) {
  temp.log <- uni %in% simper.data[[i]][,1]
  for (k in seq_along(temp.log)) {
    if (temp.log[k] == TRUE) {
      simper.frame[k,i] <- as.numeric(as.vector(simper.data[[i]][k,2]))
    }
  }
}

final.data <- simper.frame


}

