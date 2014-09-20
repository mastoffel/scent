# Aim of this script: Graph that shows correlation with heterozygosity for increasing number of loci

library(Rhh)
library(plyr)
library(ggplot2)
library(ggthemes)

# load own function , which is different from mlh() since it 
# loads from workspace, doesn´t save as text file and 
# just calculates SH
source("mlhWS.R")
source("LociValidation.R")
source("SumResults.R")
source("HHLociValidation.R")
source("LociValidationPups.R")

# load PC2 and NumbComp data
CorVars <- read.csv("het_df.csv",row.names=1) 
CorVars <- het.df


# load genotypes file with all 46 loci, which is ordered according to ID
genotypes <- read.table(".\\txt\\raw_41loci_ordered.txt", na.strings = "NA")


# slightly changed procedure for factor analysis
het.df$F1F2 <- het.df$F1 + het.df$F2
CorVars <- het.df
# run algorithm to get results (2 for NumberCompounds, 4 for PC 2)
results41 <- LociValidation(genotypes,CorVars,2,100) 
results41.FA <- LociValidation(genotypes,CorVars,7,100) 

results41.Pups <- LociValidationPups(genotypes,CorVars,2,100) 
results41.FA.Pups <- LociValidationPups(genotypes,CorVars,7,100) 

resultsHH41 <- HHLociValidation(genotypes, numIter=100)

# ddply for summarizing results

#labels <- colnames(results)[2:dim(results)[2]]
# summarizing in another df
allresults41 <- SumResults(results41)
allresults41.FA <- SumResults(results41.FA)
allresults41.Pups <- SumResults(results41.Pups)
allresults41.FA.Pups <- SumResults(results41.FA.Pups)

allresultsHH41 <- SumResults(resultsHH41)



# allresults$locnum <- as.factor(allresults$locnum)

# plotting stuff
ggplot(allresults41, aes(x = locnum, y = cormean)) +
        geom_errorbar(aes(ymin = cormean-corsd, ymax = cormean+corsd),
                      width=0.8, alpha=0.7, colour="blue") +
        geom_point() +
        geom_line(colour = "blue", size=1) 


# create dataframe with both old and new genytype data
res1 <- allresults41
res2 <- allresults41.FA
res3 <- allresults41.Pups
res4 <- allresults41.FA.Pups

resall <- rbind(res1,res2,res3,res4)
resall$fac <- as.factor(c(rep(1, 41),rep(2, 41),rep(3, 41),rep(4, 41)))
levels(resall$fac) <- c("Compound richness - Mothers", "F1F2 Score - Mothers",
                        "Compound richness - Pups", "F1F2 Score - Pups")
row.names(resall) <- NULL

ggplot(resall, aes(x = locnum, y = cormean, group = fac, colour = fac,shape = fac)) +
        geom_errorbar(aes(ymin = cormean-se, ymax = cormean+se),
                      width=0.8, alpha=0.7) +
        geom_point(size = 2) +
        geom_line() +
        theme_minimal() +
        ylab(expression(r [mean +- se])) +
        xlab("number of loci") +
        theme(legend.title=element_blank()) +
        scale_colour_brewer(type="qual",palette="Set1") +
        scale_shape_manual(values=c(1,19,0,15)) +
        ggtitle("Heterozygosity correlation with compound richness / F1F2 \n for an increasing number of loci\n")


# plotting het-het

# plot
ggplot(allresultsHH41, aes(x = locnum, y = cormean)) +
        geom_errorbar(aes(ymin = cormean-se, ymax = cormean+se),
                      width=0.8, alpha=0.7) +
        geom_point(size = 3) +
        geom_line() +
        theme.paper +
        ylab(expression(r [mean+-se])) +
        xlab("number of loci") +
        theme(legend.title=element_blank()) +
        ggtitle("het-het correlation \n(100 loci-resamplings per data point) \n ")
