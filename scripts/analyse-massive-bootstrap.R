# analyzing best substances from bootstrap analysis

library(ggplot2)
require(dplyr)
require(magrittr)
library(vegan)
source("multiplot.R")
library(reshape2)

# load compounds from bootstrap!
best_mums <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                  "MSc.Behaviour\\Research\\Seal Scent\\",
                                  "R code\\data\\csv_files\\",
                                  "bootstrap_mums.csv", sep = ""),
                                   row.names=1)
best_pups <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                            "MSc.Behaviour\\Research\\Seal Scent\\",
                            "R code\\data\\csv_files\\",
                            "bootstrap_pups.csv", sep = ""),
                             row.names=1)



# check out the intersect
intersect(row.names(best_mums)[1:30], row.names(best_pups)[1:30]) # "20.361875","11.84727273","8.307142857" "20.568" ,"21.12214286", "8.516153846", "19.86137931"
best_df <- rbind(best_mums, best_pups)
best_df$age <- c(rep(1, nrow(best_mums)), rep(2, nrow(best_pups)))
best_df$x <- c(1:nrow(best_df))

# occurences plotted
ggplot(best_df, aes(x = x, y = occurences, colour = age)) +
        geom_point() +
        theme_bw(base_size = 20)

# significance development------------------------------------------------------

# loading
scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "MSc.Behaviour\\Research\\Seal Scent\\R code\\",
                                                  "data\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))
# relatedness matrix (old: relatedness_41loci.csv)
relatedness <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                              "MSc.Behaviour\\Research\\Seal Scent\\",
                              "R code\\data\\csv_files\\",
                              "relatednessnew.csv", sep = ""),
                              row.names=1)

# diversity measures
scent_diversity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                  "MSc.Behaviour\\Research\\Seal Scent\\",
                                  "R code\\data\\csv_files\\",
                                  "scent diversity.csv", sep = ""),
                            row.names=1)

## heterozygosity SH
heterozygosity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                 "MSc.Behaviour\\Research\\Seal Scent\\",
                                 "R code\\data\\csv_files\\",
                                 "heterozygosity_41loci.csv", sep = ""),
                           row.names=1) 

# beach and family factor
factors <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                          "MSc.Behaviour\\Research\\Seal Scent\\",
                          "R code\\data\\csv_files\\",
                          "factors.csv", sep = ""),
                          row.names=1) 

# subset------------------------------------------------------------------------

mums <- c(rep(TRUE,41),rep(FALSE,41)) # mums
pups <- c(rep(FALSE,41),rep(TRUE,41)) # pups

# gendiff (1-rel) because bray curtis dist in R is similarity (not dissimilarity)
abundm <- scent_abundance[mums, ]
relm <-  1-relatedness[mums, mums]

abundp <- scent_abundance[pups, ]
relp <-  1-relatedness[pups, pups]


# get vectors of best substance names
sub_names_mums <- row.names(best_mums)
sub_names_pups <- row.names(best_pups)

statm <- vector()
sigm <- vector()
statp <- vector()
sigp <- vector()

for (i in 2:100) {
        bc_dist <- vegdist(abundm[, sub_names_mums[1:i]], method = "bray")
        mod <- mantel(relm, bc_dist, na.rm = T, method = "spearman")
        statm <- append(statm, mod$statistic)
        sigm <- append(sigm, mod$sig)
}

for (i in 2:100) {
        bc_dist <- vegdist(abundp[, sub_names_pups[1:i]], method = "bray")
        mod <- mantel(relp, bc_dist, na.rm = T, method = "spearman")
        statp <- append(statp, mod$statistic)
        sigp <- append(sigp, mod$sig)
}

mantelRs <- data.frame("statm" = statm, "statp" = statp)
mantelsigs <- data.frame("sigm" = sigm, "sigp" = sigp)       


# longformat allstat

lf_mantel <- melt(mantelRs)
lf_sig <- melt(mantelsigs)

lf_mantel$variable <- as.factor(lf_mantel$variable)
lf_mantel$x <- rep(1:99, each = 1)
lf_sig$variable <- as.factor(lf_sig$variable)
lf_sig$x <- rep(1:99, each = 1)

library(grid)
# simple plot
ggplot(lf_mantel[1:99, ], aes(x = x, y = value, colour = variable)) +
        geom_line(size = 1.5) +
        theme_classic(base_size = 16) +
        theme(legend.title = element_blank()) +
        xlab("best substances from bootstrap bioenv") +
        ylab("mantelR") +
        scale_colour_discrete(labels = c("mothers", "pups"))

# plot for presentation
ggplot(lf_mantel[1:99, ], aes(x = x, y = value, colour = variable)) +
        geom_point(colour = "black", size = 3) +
        theme_minimal(base_size = 26) +
        geom_line(colour = "green") +
        geom_smooth(method = "loess", size = 0.5) +
        theme(strip.text.x = element_text(vjust=1,size = 18),
              axis.title.x = element_text(vjust= -2 ,size = 28),
              axis.title.y = element_text(vjust=3,size = 28),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = (unit(c(.5, .5, 2, 2), "cm"))) +
       # scale_x_continuous(breaks=c(seq(from = 0.8, to = 1.20, by = 0.1))) +
        #geom_text(aes(0.85,80, label="(a) r = 0.34, p = 0.027"),size=4) +
        xlab("cumulative substances from bootstrap") +
        ylab("mantelR") 

# plot for paper with smoothed line---------------------------------------------
# make df with mantel r´s and 

#span maybe 0.13
ggplot(lf_mantel[1:99, ], aes(x = x, y = value)) +
        stat_smooth(se = FALSE, span = 0.13, size = 1.5, method = "loess") +
        geom_point(colour = "black", size = 3, alpha = 0.4) +
        theme_minimal(base_size = 26) +
        theme(strip.text.x = element_text(vjust=1,size = 18),
              axis.title.x = element_text(vjust= -2 ,size = 28),
              axis.title.y = element_text(vjust=3,size = 28),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = (unit(c(.5, .5, 2, 2), "cm"))) +
        # scale_x_continuous(breaks=c(seq(from = 0.8, to = 1.20, by = 0.1))) +
        #geom_text(aes(0.85,80, label="(a) r = 0.34, p = 0.027"),size=4) +
        xlab("cumulative substances from bootstrap") +
        ylab("mantelR") 

# best vs. factorloading---------------------------------------------------------

# get indices of compounds within scent_abundance
comp_ind_m <- which(colnames(scent_abundance) %in% sub_names_mums[1:9])
comp_ind_p <- which(colnames(scent_abundance) %in% sub_names_pups[1:20])

# PRIMER ready
paste(comp_ind_m, collapse = ",")

# correlation matrix
# pairs(scent_abundance[, comp_ind_m])
# pairs(scent_abundance[, comp_ind_p])

# factor analysis
library(HDMD)
scent_fa <- factor.pa.ginv(scent_abundance, nfactors = 4, prerotate=F,
                           rotate = "varimax", scores = T, m=4)

# extract loadings to new matrix
load    <- scent_fa$loadings
loaddf <- as.data.frame(load[, 1:4])
row.names(loaddf) <- names(scent_abundance)

# add a factor for marking best results in plot
loaddf$best <- 0
loaddf$best[comp_ind_m] <- 1
loaddf$best[comp_ind_p] <- 2
loaddf$best[intersect(comp_ind_m, comp_ind_p)] <- 3
loaddf$best <- as.factor(loaddf$best)
levels(loaddf$best) <- c("unrelevant", "best mums", "best pups", "mums and pups")

# add factor for number of occurences
occ <- apply(scent_abundance, 2, function(x) out <- sum(x > 0))
loaddf$occ <- occ

# get beach and mp factors for simper
beach <- as.factor(factors$Beach)

simper_beach <- simper(scent_abundance, beach)
simper_beach_names <- rownames(summary(simper_beach)[[1]])

simper_beach_ind <- which(names(scent_abundance) %in% simper_beach_names[1:15])
simper_mp_ind <- c(58,60,68,74,86,90,96,107,164,181,189,209) # best 2 from each pair

# add beach and mp best substances factor
loaddf$bestenv <- 0
loaddf$bestenv[simper_beach_ind] <- 1
loaddf$bestenv[simper_mp_ind] <- 2
loaddf$bestenv[intersect(simper_beach_ind, simper_mp_ind)] <- 3
loaddf$bestenv <- as.factor(loaddf$bestenv)
levels(loaddf$bestenv) <- c("irrelevant", "beach", "mumpup", "Both")




# plotting best relatedness comps vs factors------------------------------------
allfig <- list()

for (i in c("F1", "F2", "F3", "F4")) {
         
         loaddfsort <- loaddf[order(loaddf[, i]), ]
         
         guide = FALSE
         if (i == "F4") guide="legend"

        ggplot(data=loaddfsort, aes_string(x = i)) +
                        geom_point(size = 4, aes(y = 1:213, colour = best, alpha = best)) +
                        
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
