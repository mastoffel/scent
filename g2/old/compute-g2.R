
# selfing rate (g2) from microsat genotype data---------------------------------
data <- read.table("raw_41loci_ordered.txt", row.names = 1)

# turn data into 0 (homozygote), 1 (heterozygote) or -1 (NA on one locus)
s1 <- seq(1, ncol(data), 2)

checkhet <- function(x) {
        
        s1 <- seq(1, length(x), 2)
        newx <- as.vector(rep(NA, length(x)/2))
        count <- 1
        
        for(i in s1){
                
                if (is.na(x[i] | x[i + 1])) {
                        
                        newx[count] = -1
                        count = count + 1
                        
                } else if (x[i] == x[i + 1]) {
                        newx[count] = 0
                        count = count + 1
                        
                } else if (x[i] != x[i + 1]) {
                        newx[count] = 1
                        count = count + 1
                }
        }
        newx 
}

# transform in
het <- t(apply(data, 1, checkhet))


# simple variant for no missing data--------------------------------------------

# missing vals are treated as homozygous
het[-1] <- 0

allsum <- 0

for (i in 1:ncol(het)){
        loc_vec <- as.vector(het[, i])
        loc_mat <- as.matrix(het[, -i])
        sumvm <- loc_vec * loc_mat
        allsum <- allsum + (sum(colSums(sumvm)))
}

numerator <- allsum

denominator <- 0
est <- 1/(nrow(het) - 1)

for (i in 1:ncol(het)){
        
        loc_vec <- as.vector(het[, i])
        to_loop <- c(1:ncol(het))[-i]
        
        sum_j <- 0
        
        for (k in to_loop) {
                
                loc_vec_comp <- as.vector(het[, k])
                
                sum_inds <- 0
                
                for (n in 1:length(loc_vec)) {
                        
                        ind_sum <- sum(loc_vec[n] * loc_vec_comp[-n])
                        
                        sum_inds <- sum_inds + ind_sum
                        }
                
                sum_j <- sum_j + (est * sum_inds)      
                }
        
         denominator <- denominator + sum_j           
        }
     


g2hat <- numerator/denominator - 1

selfrate <- (1 + (5 * g2hat) - (sqrt(1 + 10*g2hat + 9*g2hat))) / (2*g2hat)




