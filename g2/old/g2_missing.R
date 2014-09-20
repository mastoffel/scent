# selfing rate (g2) from microsat genotype data---------------------------------
data <- read.table("raw_41loci_ordered.txt", row.names = 1)

        
# turn data into 0 (homozygote), 1 (heterozygote) or -1 (NA on one locus)

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
origin <- apply(data, 1, checkhet)

# define matrix with 1 for missing and 0 for all others
m <- origin
m[m==1] <- 0
m[m==-1] <- 1

# H matrix with 0 for -1
h <- origin
h[h==-1] <- 0

n <- ncol(origin) # number of individuals
l <- nrow(origin) # number of loci

# mij: proportion of individuals missing on i and j ´s locus
mtemp <- m %*% t(m)
m_ij <- mtemp/n

# vector with rowsums
m_loc_temp <- apply(m, 1, sum)
m_loc <- m_loc_temp / n

# numerator --------------------------------------------------------------------
# pij entry says total amount of individuals that het locus i and locus j
p <- h %*% t(h) 
numerator <- 0
missmat <- matrix(rep(0, l*l), ncol = l)

for (i in seq(1:nrow(h))){
        for (j in seq(1:nrow(h))[-i]){
                missmat[i,j]  <- 1/ (n * (1 - m_loc[i] - m_loc[j] + m_ij[i,j])) 
        }
}

numerator_mat <- p * missmat
numerator <- sum(numerator_mat)

# denominator-------------------------------------------------------------------
denominator_mat <- matrix(rep(0, l*l), ncol = l)

# sum over loci
for (i in seq(1:nrow(h))){
        for (j in seq(1:nrow(h))[-i]){
                sum_miss  <- 1/(((n * (n - 1) * (1 - m_loc[i] - m_loc[j] + m_loc[i] * m_loc[j]))) - 
                                        (n * (m_ij[i, j] - m_loc[i] * m_loc[j])))
                
                # sum over individuals
                sum_ind <- 0
                hk1 <- h[i, ]
                hk2 <- h[j, ]
                mat <- hk1 %*% t(hk2)
                diag(mat)  <- 0
                sum_ind <- sum(mat)
                
                denominator_mat[i,j] <- sum_miss * sum_ind
                
        }
}

denominator <- sum(denominator_mat)

g2 <- (numerator / denominator) - 1












