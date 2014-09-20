# load new relatedness values in strange vector form from athinas program
relatedness.raw <- read.csv(".\\csv_files\\relatedness_raw.csv")

# "copy" (old) relatedness matrix
relate_new <- relatedness
relate_new[, ] <- NA
# loop through matrix and fill non-NA values with new ones 
for (i in row.names(relatedness)) {
        for (k in names(relatedness)) {
                if (!is.na(relatedness[i,k])) {
                        rowind <- which(relatedness.raw$Ind1 == i & relatedness.raw$Ind2 == k)
                        # may be otherway rounf
                        if (length(rowind) == 0) {
                                rowind <- which(relatedness.raw$Ind1 == k & relatedness.raw$Ind2 == i)    
                        }
                        relateval <- relatedness.raw[rowind,3] # third column contains values
                        relate_new[i,k] <- relateval
                }
        }
}


library(rJava)
library(xlsx)

write.xlsx(relate_new,"relatedness_new.xlsx",row.names=TRUE,col.names=TRUE)
