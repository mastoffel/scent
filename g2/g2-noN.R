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
m_ij <- mtemp

# vector with rowsums
m_loc_temp <- apply(m, 1, sum)
m_loc <- m_loc_temp 

# numerator --------------------------------------------------------------------
# pij entry says total amount of individuals that het locus i and locus j
p <- h %*% t(h) 
numerator <- 0
missmat <- matrix(rep(0, l*l), ncol = l)

for (i in seq(1:nrow(h))){
        for (j in seq(1:nrow(h))[-i]){
                missmat[i,j]  <- 1/ (n - m_loc[i] - m_loc[j] + m_ij[i,j]) 
        }
}

numerator_mat <- p * missmat
numerator <- sum(numerator_mat)

# denominator-------------------------------------------------------------------
denominator_mat <- matrix(rep(0, l*l), ncol = l)
missmat2 <- matrix(rep(0, l*l), ncol = l)

# sum over loci
for (i in seq(1:nrow(h))){
        for (j in seq(1:nrow(h))[-i]){
                missmat2[i,j] <- 1/((((n - 1) * (n - m_loc[i] - m_loc[j] + 1/n * m_loc[i] * m_loc[j]))) - 
                                            ((m_ij[i, j] - 1/n * m_loc[i] * m_loc[j])))
                
        }
}

nullmat <- matrix(rep(1, n*n), ncol=n)
diag(nullmat) <- 0
q <- h %*% (nullmat %*% t(h))

denominator_mat <- missmat2 * q
denominator <- sum(denominator_mat)

g2 <- (numerator / denominator) - 1