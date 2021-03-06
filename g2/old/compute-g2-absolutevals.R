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

# original full data matrix
origin <- (apply(data, 1, checkhet))

# H matrix with 0 for -1
h <- origin
h[h==-1] <- 0

# define matrix with 1 for missing and 0 for all others
m <- origin
m[m==1] <- 0
m[m==-1] <- 1

# mij: proportion of individuals missing on i and j �s locus
n <- ncol(h)
l <- nrow(h)
m_ij <- m %*% t(m)


# vector with rowsums
m_loc <- apply(m, 1, sum)


# numerator --------------------------------------------------------------------
# pij entry says total amount of individuals that het locus i and locus j
p <- h %*% t(h) 
missmat <- 1/(n - m_ij)
diag(missmat) <- 0

numerator_mat <- p * missmat

numerator <- sum(numerator_mat)

# denominator-------------------------------------------------------------------

missmat2 <- matrix(rep(0, l*l), ncol = l)

# sum over loci
for (i in seq(1:nrow(h))){
        for (j in seq(1:nrow(h))[-i]){
                
                sum_miss  <- 1/((n * (n - 1) - m_loc[i]*m_loc[j] + m_ij[i, j]))
                                  
                # sum over individuals
                sum_ind <- 0
                hk1 <- h[i, ]
                hk2 <- h[j, ]
                mat <- hk1 %*% t(hk2)
                diag(mat)  <- 0
                sum_ind <- sum(mat)
                
                missmat2[i,j] <- sum_miss * sum_ind
                
        }
}

denominator_mat <- missmat2

denominator <- sum(denominator_mat)

g2 <- (numerator / denominator) - 1













