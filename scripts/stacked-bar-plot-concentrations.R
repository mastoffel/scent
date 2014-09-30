# goal: stacked bar plots to show differences in relative concentrations between colonies

scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "projects\\sealscent\\data_files\\",
                                                  "Rdata\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))

# beach and family factor
factors <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                          "projects\\sealscent\\data_files\\",
                          "Rdata\\csv_files\\",
                          "factors.csv", sep = ""),
                                 row.names=1) 

# elements
simper_mp_ind <- c(58, 68,  86,  90, 106, 107, 164)
comp_ind_m <- c(36,52,86,88,96,103,110,203,206,207)
# add beach factor
scent <- cbind(scent_abundance, factors$Beach)
names(scent)[length(names(scent))] <- "beach"

# bar plots with juxtaposed mums and pups
library(gplots)


par(mfrow=c(9,2), mai = c(0.1,0.4,0.3,0.1))
# alternative: k in simper_mp_ind
for (k in comp_ind_m){
        for (i in 1:2) {
                
                comp <- scent[scent$beach==i, k]
                len <- length(comp)
                comp_mat <- matrix(rep(NA, len), ncol = len/2)
                
                comp_mat[1, ] <- comp[1:(len/2)]
                comp_mat[2, ] <- comp[(len/2+1):len]
                
                
                barplot2(comp_mat, beside = TRUE, col = c("grey12", "grey82"),
                         legend = c("mums", "pups"), space = c(0, 2),
                         main = names(scent_abundance[k]), ylim = c(0,5))
                
        }
}


