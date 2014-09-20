DimExplore <- function(y, x, method = "fa", distance = "bray", subset = "all") {
## function to explore different dimension reduction / ordination methods for similarity matrix models
## 
## libraries
library(psych)
library("ggplot2")
library("grid")
library("vegan")
library("MASS") #required for vegan
library("HDMD")
## functions
source("PCADiff.R")
source("SubsetAll.R")
source("MinMod.R")
#source("ResultsRelatedness.R")
#source("resultsHet.R")
#source("multiplot.R")


## predefine
x <- as.data.frame(x)
AbundMat <- x
results <- vector()

## subset

if (subset == "mums") {ind <- 1:41}
if (subset == "pups") {ind <- 42:82}
if (subset == "all") {ind <- 1:82}


for (numDim in c(1:10)) {
        
        ## define function fo the chosen ordination method with increasing factor number
        if (method == "pcoa") {
                ordiscores <- function(AbundMat,numDim){
                        score <- as.data.frame(cmdscale(vegdist(AbundMat, method = distance), k = numDim))   
                }
                
        } else if (method == "fa") {
                ordiscores <- function(AbundMat, numDim){
                        model <- factor.pa.ginv(AbundMat, nfactors = numDim, prerotate=T,
                                                   rotate = "promax", scores = T, m=4)
                        score <- as.data.frame(model$scores) 
                }
        } else if (method == "pca") {
                ordiscores <- function(AbundMat, numDim){
                        model <- prcomp(AbundMat)
                        score <- as.data.frame(scores(model, choices = 1:numDim))               
                }
                
        } else if (method == "mds") {
                ordiscores <- function(AbundMat, numDim){
                        model <- metaMDS(AbundMat, distance = "bray", k = numDim, trymax = 50)
                        score <- as.data.frame(scores(model))
                }
        }
        
        
        # get factor scores
        tempscores <- ordiscores(AbundMat, numDim)
        
        # create data frame
        alldf <- as.data.frame(cbind(y, tempscores))
        
        ## get chosen subset
        alldf <- alldf[ind, ]
        bestmod <- MinMod(alldf)[[2]]
        devExpl <- (bestmod$null.deviance - bestmod$deviance)/bestmod$null.deviance
        results[numDim] <- devExpl
        
}

results

}