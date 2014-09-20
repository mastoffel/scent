## analyses on Simper
# input: list of( list of data.frames with 2 columns each(element + explained variance) + list of group labels )
# output: data frame with elements from simper as rows, groups as columns and explained variance as values

best.simper.df <- function(simper.groups.list, number.extract) {

simper.data <- simper.groups.list[[1]] # first part of input list is simper data
groups <- simper.groups.list[[2]] # second part is character vector with group names

extract <- function(x) return(x[1:number.extract,1]) # extract given number of components from simper analysis

# extract alsways first (two, three, four etc.) elements
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
  temp.log <-  simper.data[[i]][,1] %in% uni #get vector of logicals for fitting indices
  
  for (k in seq_along(temp.log)) {
    if (temp.log[k] == TRUE) {
      ind <- which(simper.data[[i]][k,1] == uni)
      simper.frame[ind,i] <- as.numeric(as.vector(simper.data[[i]][k,2]))
    }
  }
}

final.data <- simper.frame


}

