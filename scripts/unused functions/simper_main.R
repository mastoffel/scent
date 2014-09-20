# analysing simper

rm(list = ls())

# when time use assign() to make function out of here and asssign names

# loading functions
source("best.simper.list.R")
source("best.simper.df.R")
source('element.indices.R')

for (i in 1:5) {
# getting list of simpers

simper.list <- best.simper.list("simper.csv")

#getting data.frame with best simper results

best.simper <- best.simper.df(simper.list,i)
name <- paste("best.simper", i, sep=".")
assign(name,best.simper)

#get indices

ind <- element.indices(best.simper)
name.ind <- paste("ind", i, sep = ".")
assign(name.ind,ind)
}

# arranging simper outputs for mum-pup- pairs in a way that sorts the elements by number of occurences in mum-pup-pair simpers
count.occurences <- apply(best.simper, 1, function(x) sum(!is.na(x)))
sorted.count.occurences <- sort(count.occurences, decreasing = TRUE)
occurences <- as.matrix(sorted.count.occurences)

# same for explained variance
count.var <- apply(best.simper, 1, function(x) sum(x, na.rm = TRUE))
sorted.count.var <- sort(count.var, decreasing = TRUE)
expl.var <- as.matrix(sorted.count.var)

# write data frame to xlsx
library(xlsx)
write.xlsx(x = best.simper.one, file = "simper_data_1.xlsx", sheetName = "best_one", row.names = TRUE)
write.xlsx(x = best.simper.two, file = "simper_data_2.xlsx", sheetName = "best_two", row.names = TRUE)
write.xlsx(x = best.simper.three, file = "simper_data_3.xlsx", sheetName = "best_three", row.names = TRUE)
write.xlsx(x = best.simper.four, file = "simper_data_4.xlsx", sheetName = "best_four", row.names = TRUE)
write.xlsx(x = best.simper.five, file = "simper_data_5.xlsx", sheetName = "best_five", row.names = TRUE)


