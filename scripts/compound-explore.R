# analyzing mum-pup correlation of mum-pup simper substances and xy location and
# beach substances

library(ggplot2)

# mum pup 
scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "MSc.Behaviour\\Research\\Seal Scent\\R code\\",
                                                  "data\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))
coord <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                  "MSc.Behaviour\\Research\\Seal Scent\\R code\\",
                  "data\\csv_files\\coordinates.csv",
                  sep = ""),row.names=1) 

factors <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                          "MSc.Behaviour\\Research\\Seal Scent\\",
                          "R code\\data\\csv_files\\",
                          "factors.csv", sep = ""),
                        row.names=1) 

heterozygosity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                 "MSc.Behaviour\\Research\\Seal Scent\\",
                                 "R code\\data\\csv_files\\",
                                 "heterozygosity_41loci.csv", sep = ""),
                                row.names=1) 

# subset beach
scent_beach1 <- scent_abundance[factors$Beach==2, ]

# 
abund_mums <- scent_beach1[1:21, ]
abund_pups <- scent_beach1[22:42, ]

# simper_mp_ind <- c(58,60,68,74,86,90,96,107,164,181,189,209)
# new version
simper_mp_ind <- c(58, 68,  86,  90, 98, 106, 107, 164, 181)
names(scent_abundance)[simper_mp_ind]

# plotting mum vs pup concentration in best mum-pup substances

plotM <- function(l){
        
        mums_temp <- abund_mums[, l]
        #mums_temp[mums_temp==0] <- NA
        pups_temp <- abund_pups[, l]
        #pups_temp[pups_temp==0] <- NA
        
        # check if the contain at least two pairs
        df <- cbind(mums_temp, pups_temp)
        count <- apply(df, 1, function(x) !(is.na(x[1] & x[2])))
       
        if (sum(count) >= 3) {
        plot(mums_temp, pups_temp, xlab="mums",ylab="pups",
             main=paste("element:", toString(names(scent_abundance)[l])))
        abline(lm((mums_temp ~ pups_temp)))
        
        
        }
}
# plot all mum vs. pup compound concentration for best simper mum-pup compounds
par(mfrow = c(3,4))

for (i in simper_mp_ind){
        plotM(i)
}


# coordinates-------------------------------------------------------------------

# create subset with relevant individuals (40)
sub <- row.names(coord) %in% row.names(heterozygosity) 
coordinates <- subset(coord, sub)
coord <- subset(coordinates, factors$Beach==1)
coord["M48",1] <- "0618"
coord["P48",1] <- "0618"
coord$XY <- as.character(coord$XY)
# split coord
X <- apply(coord, 1, function(x) substring(x,1,2))
Y <- apply(coord, 1, function(x) substring(x,3,4))
coord$X <- X
coord$Y <- Y
coord$X <- as.numeric(coord$X)
coord$Y <- as.numeric(coord$Y)
#-------------------------------------------------------------------------------

# plot Seal distribution
ggplot(coord, aes(X,Y)) +
        geom_point(colour="blue", size = 4, alpha = 0.5) +
        theme_bw(base_size=18)

# create scent_abundance subset with individuals with xy coordinates------------
scent <- scent_abundance[factors$Beach==1, ]

# get simper from beach
library(vegan)
beach <- as.factor(factors$Beach)
simper_beach <- simper(scent_abundance, beach)
simper_beach_names <- rownames(summary(simper_beach)[[1]])
simper_beach_ind <- which(names(scent_abundance) %in% simper_beach_names[1:15])

# add mum pup substances to coords data
coord[, 4:10] <- scent[ ,simper_mp_ind]

# plot
library(scatterplot3d)

plotxy <- function(l){
        coord2 <- coord
        # coord2[coord2==0] <- NA # without missing data
        plot <- scatterplot3d(coord2[, 2], coord2[, 3], coord2[, l], 
                              xlab="x", ylab="y",
                              zlab="concentration",
                      type = "h", highlight.3d=TRUE,
                      main=paste("element:", toString(names(coord)[l])))
             #main=paste("element:", toString(names(scent_abundance)[l])))
        fit <- lm(coord2[, l] ~ coord2[, 2] + coord2[, 3])
        plot$plane3d(fit)
}

# plot mum-pup substance concentration vs. xy coords

par(mfrow = c(3,4))

for (i in 4:10){
        plotxy(i)
}

# add beach substances to coords data
coord[, 11:25] <- scent[ ,simper_beach_ind]

par(mfrow = c(3,5))

for (i in 11:25){
        plotxy(i)
}


# plot: mum-pup r2 vs spatial model deviance------------------------------------

# extract mum vs. pup concentration r2 for best mum-pup substances
options(scipen = 8)
df_lm_mp <- data.frame("r-squared" = NA, "p-val f test" = NA)

get_stats_mp <- function(l){
        mums_temp <- abund_mums[, l]
        #mums_temp[mums_temp==0] <- NA
        pups_temp <- abund_pups[, l]
        #pups_temp[pups_temp==0] <- NA
        
        df <- cbind(mums_temp, pups_temp)
        count <- apply(df, 1, function(x) !(is.na(x[1] & x[2])))
        
        if (sum(count) >= 2) {
        model <- summary(lm(mums_temp ~ pups_temp))
        rsqu <- model$r.squared
        f <-  model$fstatistic
        p_val <- pf(f[1],f[2],f[3],lower.tail=F)
        df_lm_mp <- rbind(df_lm_mp, c(rsqu, p_val))
        
        } else {
                df_lm_mp <- rbind(df_lm_mp, c(NA, NA))
        }
        
        df_lm_mp
}

for (i in simper_mp_ind){
        df_lm_mp <- get_stats_mp(i)
}

df_lm_mp <- df_lm_mp[-1, ]
row.names(df_lm_mp) <- names(coord[, 4:10])


# extract spatial model deviance for best mum-pup substances

df_lm_xy <- data.frame("rsquared" = NA, "pval_f" = NA, "pval_x" = NA, "pval_y" = NA)

get_stats_xy <- function(l){
        coord_pres <- coord
        coord_pres[coord_pres==0] <- NA # without missing data
        model <- summary(lm(coord_pres[, l] ~ coord_pres[, 2] + coord_pres[, 3]))
        rsqu <- model$r.squared
        # pvals <- summary(lm(coord[, l] ~ coord[, 2] + coord[, 3]))$coefficients[,4][2:3]
        f <-  model$fstatistic
        xy_pval <- model$coefficients[, 4][2:3]
        p_val <- pf(f[1],f[2],f[3],lower.tail=F)
        
        df_lm_xy <- rbind(df_lm_xy, c(rsqu, p_val, xy_pval))
        df_lm_xy
}

for (i in 4:10){
        df_lm_xy <- get_stats_xy(i)
}

df_lm_xy <- df_lm_xy[-1, ]
row.names(df_lm_xy) <- names(scent[ ,simper_mp_ind])


# plotting them against each other
inds <- which(row.names(df_lm_xy) %in% row.names(df_lm_mp))
df_xy <- df_lm_xy[inds, ]
inds2 <- which(row.names(df_lm_mp) %in% row.names(df_xy))
df_mp <- df_lm_mp[inds2, ]
plot(df_xy$rsquared, df_mp$r.squared)
abline(lm(df_mp$r.squared ~ df_xy$rsquared))
summary(lm(df_mp$r.squared ~ df_xy$rsquared))


# repeat xy - conc models just for mums

xy_mums <- data.frame("rsquared" = NA, "pval_f" = NA, "pval_x" = NA, "pval_y" = NA)

get_stats_xy <- function(l){
        
        model <- summary(lm(coord[1:20, l] ~ coord[1:20, 2] + coord[1:20, 3]))
        rsqu <- model$r.squared
        # pvals <- summary(lm(coord[, l] ~ coord[, 2] + coord[, 3]))$coefficients[,4][2:3]
        f <-  model$fstatistic
        p_val <- pf(f[1],f[2],f[3],lower.tail=F)
        xy_pval <- model$coefficients[, 4][2:3]
        xy_mums <- rbind(xy_mums, c(rsqu, p_val, c(xy_pval)))
        xy_mums
}

for (i in 16:30){
        xy_mums <- get_stats_xy(i)
}

xy_mums <- xy_mums[-1, ]
row.names(xy_mums) <- names(scent[ ,simper_beach_ind])