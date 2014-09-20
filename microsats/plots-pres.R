# plot stuff for presentation

load("C:/Users/Martin/Studium/MSc.Behaviour/Research/",
     "Seal Scent/R code/Microsatellites/old/Results4Factor.RData")

library(grid)
library(ggplot2)
library(Rhh)
library(plyr)
library(ggthemes)

# load own function , which is different from mlh() since it 
# loads from workspace, doesn´t save as text file and 
# just calculates SH
source("mlhWS.R")
source("LociValidation.R")
source("SumResults.R")
source("HHLociValidation.R")
source("LociValidationPups.R")

res1 <- allresults41
res2 <- allresults41.FA
res3 <- allresults41.Pups
res4 <- allresults41.FA.Pups

all_comp_het <- rbind(res1,res3)
all_comp_het$fac <- as.factor(c(rep(1, 41),rep(2, 41)))
levels(all_comp_het$fac) <- c("Mums", "Pups")

ggplot(all_comp_het, aes(x = locnum, y = cormean, group = fac, colour = fac)) +
        geom_errorbar(aes(ymin = cormean-se, ymax = cormean+se),
                      width=0.8, alpha=0.7, size = 0.8) +
        geom_point(size = 3, shape = 16) +
        geom_line(size = 0.8) +
        theme_minimal(base_size = 20) +
        theme(axis.title.x = element_text(vjust= -2 ,size = 22),
              axis.title.y = element_text(vjust=3,size = 22),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = (unit(c(.5, .5, 2, 2), "cm"))) +
        #geom_hline(yintercept=0.305) +
        ylab(expression(r [mean+-se])) +
        xlab("number of loci") +
        theme(legend.title=element_blank()) +
        scale_colour_brewer(type="qual",palette="Set2") +
        scale_shape_manual(values=c(1,19,0,15)) 
        # ggtitle("correlation between heterozygosity and compound richness")
ggsave(file="hetall.pdf", width = 6, height = 4)

# plotting the same with factors

fa_het <- rbind(res1,res3,res2,res4)
fa_het$fac <- as.factor(c(rep(1, 41),rep(2, 41),rep(3, 41),rep(4, 41)))
levels(fa_het$fac) <- c("Mums overall profile", "Pups overall profile",
                        "Mums factor 1 + factor 2", "Pups factor 1 + factor 2")

ggplot(fa_het, aes(x = locnum, y = cormean, group = fac, colour = fac)) +
        geom_errorbar(aes(ymin = cormean-se, ymax = cormean+se),
                      width=0.8, alpha=0.7, size = 0.8) +
        geom_point(size = 3, shape = 16) +
        geom_line(size = 0.8) +
        theme_minimal(base_size = 24) +
        theme(
                axis.title.x = element_text(vjust= -2 ,size = 26),
                axis.title.y = element_text(vjust=3,size = 26),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                plot.margin = (unit(c(.5, .5, 2, 2), "cm"))) +
        #geom_hline(yintercept=0.305, linetype = "") +
        ylab(expression(r [mean+-se])) +
        xlab("number of loci") +
        theme(legend.title=element_blank(),
              legend.position="none") +
        scale_colour_brewer(type="qual",palette="Set2") +
        scale_shape_manual(values=c(1,19,0,15)) 