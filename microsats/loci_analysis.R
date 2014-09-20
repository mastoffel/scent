# Aim of this script: Graph that shows correlation with heterozygosity for increasing number of loci

library(Rhh)
library(plyr)
library(ggplot2)
library(ggthemes)
# load PC2 and NumbComp data
CorVars <- read.csv("het_df.csv",row.names=1) 

# load own function , which is different from mlh() since it loads from workspace, doesn´t save as text file and 
# just calculates SH
source("mlhWS.R")
source("LociValidation.R")
source("SumResults.R")
source("HHLociValidation.R")
# load genotypes file with all 46 loci, which is ordered according to ID
genotypes <- read.table("raw_41loci_ordered.txt", na.strings = "NA")
genotypes2 <- read.table("raw_46loci_ordered.txt", na.strings = "NA")

# run algorithm to get results (2 for NumberCompounds, 4 for PC 2)
results41 <- LociValidation(genotypes,CorVars,2,100) 
results46 <- LociValidation(genotypes2,CorVars,2,100)
results41.PC <- LociValidation(genotypes,CorVars,4,100) 
results46.PC <- LociValidation(genotypes2,CorVars,4,100)

resultsHH41 <- HHLociValidation(genotypes, numIter=500)
resultsHH46 <- HHLociValidation(genotypes2, numIter=500)

# ddply for summarizing results

#labels <- colnames(results)[2:dim(results)[2]]
# summarizing in another df
allresults41 <- SumResults(results41)
allresults46 <- SumResults(results46)
allresults41.PC <- SumResults(results41.PC)
allresults46.PC <- SumResults(results46.PC)

allresultsHH41 <- SumResults(resultsHH41)
allresultsHH46 <- SumResults(resultsHH46)



# allresults$locnum <- as.factor(allresults$locnum)

# plotting stuff
ggplot(allresults41, aes(x = locnum, y = cormean)) +
        geom_errorbar(aes(ymin = cormean-corsd, ymax = cormean+corsd),
                      width=0.8, alpha=0.7, colour="blue") +
        geom_point() +
        geom_line(colour = "blue", size=1) 
        

# create dataframe with both old and new genytype data
res1 <- allresults41
res2 <- allresults46[1:41, ]
res3 <- allresults41.PC
res4 <- allresults46.PC[1:41, ]
resall <- rbind(res1,res2,res3,res4)
resall$fac <- as.factor(c(rep(1, 41),rep(2, 41),rep(3,41),rep(4,41)))
levels(resall$fac) <- c("41 loci","46 loci","41 loci PC2","46 loci PC2")
row.names(resall) <- NULL

ggplot(resall, aes(x = locnum, y = cormean, group = fac, colour = fac)) +
        geom_errorbar(aes(ymin = cormean-se, ymax = cormean+se),
                      width=0.8, alpha=0.7) +
        geom_point() +
        geom_line() +
        theme_minimal(base_size=18) +
        ylab("r [mean+-se]") +
        xlab("number of loci") +
        theme(legend.title=element_blank()) +
        ggtitle("correlation between heterozygosity and number of compounds \n(100 random resamplings per data point)")

# create dataframe for hethet data
res1HH <- allresultsHH41
res2HH <- allresultsHH46[1:41, ]
resallHH <- rbind(res1HH,res2HH)
resallHH$fac <- as.factor(c(rep(1, 41),rep(2, 41)))
levels(resallHH$fac) <- c("41 loci","46 loci")

# plot
ggplot(resallHH, aes(x = locnum, y = cormean, group = fac, colour = fac)) +
        geom_errorbar(aes(ymin = cormean-se, ymax = cormean+se),
                      width=0.8, alpha=0.7) +
        geom_point() +
        geom_line() +
        theme_minimal(base_size=18) +
        ylab("r [mean+-se]") +
        xlab("number of loci") +
        theme(legend.title=element_blank()) +
        ggtitle("het-het correlation \n(1000 random resamplings per data point)")
