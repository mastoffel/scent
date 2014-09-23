# plotting most important substances vs. factor loadings

# analyzing best substances from bootstrap analysis

library(ggplot2)
require(dplyr)
require(magrittr)
library(vegan)
source("multiplot.R")
library(reshape2)

# load compounds from bootstrap!
best_mums <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                            "projects\\sealscent\\data_files\\",
                            "Rdata\\csv_files\\",
                            "bootstrap_mums.csv", sep = ""),
                      row.names=1)



# significance development------------------------------------------------------

scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "projects\\sealscent\\data_files\\",
                                                  "Rdata\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))

# diversity measures
scent_diversity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                  "projects\\sealscent\\data_files\\",
                                  "Rdata\\csv_files\\",
                                  "scent diversity.csv", sep = ""),
                            row.names=1)

# relatedness matrix (old: relatedness_41loci.csv)
relatedness <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                              "projects\\sealscent\\data_files\\",
                              "Rdata\\csv_files\\",
                              "relatednessnew.csv", sep = ""),
                        row.names=1)

## heterozygosity SH
heterozygosity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                 "projects\\sealscent\\data_files\\",
                                 "Rdata\\csv_files\\",
                                 "heterozygosity_41loci.csv", sep = ""),
                           row.names=1) 

# beach and family factor
factors <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                          "projects\\sealscent\\data_files\\",
                          "Rdata\\csv_files\\",
                          "factors.csv", sep = ""),
                    row.names=1) 

# subset------------------------------------------------------------------------

mums <- c(rep(TRUE,41),rep(FALSE,41)) # mums


# gendiff (1-rel) because bray curtis dist in R is similarity (not dissimilarity)
abundm <- scent_abundance[mums, ]
relm <-  1-relatedness[mums, mums]

# best vs. factorloading---------------------------------------------------------

sub_names_mums <- row.names(best_mums)

# get indices of compounds within scent_abundance
comp_ind_m <- which(colnames(scent_abundance) %in% sub_names_mums[1:9])

# PRIMER ready
paste(comp_ind_m, collapse = ",")

# correlation matrix
# pairs(scent_abundance[, comp_ind_m])
# pairs(scent_abundance[, comp_ind_p])

# factor analysis
library(HDMD)
scent_fa <- factor.pa.ginv(scent_abundance, nfactors = 4, prerotate=T,
                           rotate = "promax", scores = T, m=4)

# extract loadings to new matrix
load    <- scent_fa$loadings
loaddf <- as.data.frame(load[, 1:4])
row.names(loaddf) <- names(scent_abundance)

# add a factor for marking best results in plot
loaddf$best <- 0
loaddf$best[comp_ind_m] <- 1

loaddf$best <- as.factor(loaddf$best)
levels(loaddf$best) <- c("unrelevant", "best mums")

# add factor for number of occurences
occ <- apply(scent_abundance, 2, function(x) out <- sum(x > 0))
loaddf$occ <- occ

# get beach and mp factors for simper
beach <- as.factor(factors$Beach)

simper_beach <- simper(scent_abundance, beach)
simper_beach_names <- rownames(summary(simper_beach)[[1]])

simper_beach_ind <- which(names(scent_abundance) %in% simper_beach_names[1:15])
# top 5 from each pair sorted by occurences, than 7 comps occured in more than 10 pairs
simper_mp_ind <- c(58, 68,  86,  90, 106, 107, 164) 

# add beach and mp best substances factor
loaddf$bestbeach <- 0
loaddf$bestbeach[simper_beach_ind] <- 1
loaddf$bestbeach <- as.factor(loaddf$bestbeach)
levels(loaddf$bestbeach) <- c("unrelevant", "beach")

loaddf$bestmp <- 0
loaddf$bestmp[simper_mp_ind] <- 1
loaddf$bestmp <- as.factor(loaddf$bestmp)
levels(loaddf$bestmp) <- c("unrelevant", "mum-pup")

# plotting best relatedness comps vs factors------------------------------------
allfig <- list()

for (i in c("F1", "F2", "F3", "F4")) {
        
        loaddfsort <- loaddf[order(loaddf[, i]), ]
        
        guide = FALSE
        if (i == "F4") guide="legend"
        
        ggplot(data=loaddfsort, aes_string(x = i)) +
                geom_point(size = 4, aes(y = 1:213, colour = bestmp, alpha = bestmp)) +
                
                scale_colour_manual(values = c("grey", "red", "blue",
                                               "green"),
                                    guide = guide) +
                scale_alpha_manual(values = c(0.1, 0.6, 0.6, 0.6),
                                   guide = FALSE) +
                theme_classic(base_size=14) +
                scale_x_continuous(limits = c(-1, 1),
                                   breaks = seq(-0.8,0.8, 0.2)) +
                xlab(paste(i," - loadings", sep="")) +
                ylab("substances") +
                theme(legend.position=c(0.8,0.38)) +
                ggtitle(paste("Intersect", i, "- best relatedness",
                              sep = " "))
        
        allfig[[length(allfig) + 1]] <- last_plot()        
}

source("multiplot.R")
multiplot(allfig[[1]], allfig[[2]], allfig[[3]], allfig[[4]], cols = 2)

# plotting simper mum pup comps vs. factors---------------------------------------

allfig2 <- list()

for (i in c("F1", "F2", "F3", "F4")) {
        
        loaddfsort <- loaddf[order(loaddf[, i]), ]
        
        guide = FALSE
        if (i == "F4") guide="legend"
        
        ggplot(data=loaddfsort, aes_string(x = i)) +
                geom_point(size = 4, aes(y = 1:213, colour = bestenv, alpha = bestenv)) +
                scale_colour_manual(values = c("grey", "red", "blue",
                                               "green"),
                                    guide = guide) +
                scale_alpha_manual(values = c(0.1, 0.6, 0.6, 0.6),
                                   guide = FALSE) +
                theme_classic(base_size=14) +
                scale_x_continuous(limits = c(-1, 1),
                                   breaks = seq(-0.8,0.8, 0.2)) +
                xlab(paste(i," - loadings", sep="")) +
                ylab("substances") +
                theme(legend.position=c(0.8,0.38)) +
                ggtitle(paste("Intersect", i, "- best relatedness",
                              sep = " "))
        
        allfig2[[length(allfig2) + 1]] <- last_plot()  
        
}

source("multiplot.R")
multiplot(allfig2[[1]], allfig2[[2]], allfig2[[3]], allfig2[[4]], cols = 2)


# extracting loess for sigma plot
mod <- loess(lf_mantel$value[1:99] ~ lf_mantel$x[1:99], span = 0.2)
plot(lf_mantel$value[1:99] ~ lf_mantel$x[1:99])
plot(mod)
lines(predict(mod))
xfit <- predict(mod, newdata = seq(0,100,0.1), band="local")
write.csv(xfit, "loess.csv")
