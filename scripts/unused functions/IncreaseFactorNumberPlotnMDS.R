## calculating explained deviance for mum-pup, beach, het, relate for factor solutions 1-10
## important: comment(disable) loading of FAScores function in script "ScentResults.R"

## maybe problem with ParsModel
library(psych)

################################################
## get Heterzygosity and Relatedness for MUMS ##
################################################
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 
source("MinMod.R")
allrelresults <- list()
allhetdf <- list()
allreldf <- list()
for (numFac in c(1:10)) {
        ## change number of extracted factors
        FAScores <- function(abundance.matrix){
                # computes pca.scores for first two PC´s from relative abundance matrix
                scent.pca <- as.data.frame(cmdscale(dist(abundance.matrix), k = numFac))
               # PCAScores <- as.data.frame(scent.pca$points)
                # PCAScores
                # scree plot --> 2 PC seem to be optimal
                # plot(scent.pca,type="lines")
        }
        
        source("ScentResultsNoScores.R")
        res <- ScentResultsNoScores("mums","fa",1)
        # subset results list output
        relate.df <- res[[1]]
        results.relatedness <- res[[2]]
        #relate.list <- res[[3]]
        #relatedness <- res[[4]]
        het.df <- res[[5]]
        #results.het <- res[[6]]
        #scores <- res[[7]]
        factors <- res[[8]]
        
        ## list of relatedness results
        allrelresults <- append(allrelresults, list(results.relatedness))
        
        ## list of relatedness data frames for partial mantel
        allreldf <- append(allreldf, list(relate.df))
        
        ## list of heterozygosity data frames (het in first column + factors in other columns)
        allhetdf <- append(allhetdf, list(het.df[, -2]))
        
}

#names(allhetdf[[1]])[2] <- "F1"

## small function to extract best model and deviance explained
devExpl <- function(df) {
        df <- as.data.frame(df)
        ## get deviance from het model
        model <- MinMod(df)
        model <- model[[2]]
        # get deviance explained --> (Null dev - Res dev) / Null dev
        exp <- (model$null.deviance - model$deviance)/model$null.deviance
        exp
}

##### get deviance explained for het models ######
hetDevExpl <- sapply(allhetdf, devExpl) # alternative:lapply



## get MantelR for partial Mantel-test
## find highest (most negative) MantelR
findmin <- function(rellist) {
        rel <- as.data.frame(rellist)
        ind <- which.min(rel$mantelR)
        ind
}

## vector containing the column-number of the independent factor for partial mantel
library(ecodist)
MantelMax <- sapply(allrelresults, findmin)

getPartialMantel <- function(reldf, indVar) {
        reldf <- as.data.frame(reldf)
        newOrderdf <- cbind(reldf[,1], reldf[, indVar]) # getting independent and dependent vars on the right position
        newOrderdf <- as.data.frame(cbind(newOrderdf, reldf[, -c(1,indVar)]))
        names(newOrderdf)[1] <- c("rel")
        PartialMantel <- mantel(rel ~. ,data=newOrderdf)
        PartialMantel[1]
}

allMantel <- vector()
for (i in 1:length(MantelMax)) {
        allMantel[i] <- getPartialMantel(allreldf[[i]],MantelMax[i])
}

##### all Partial Mantel values #####
allMantel

detach(package:ecodist)




#####################################################
##get Beach and Mum Pup deviances for ALL samples ##
####################################################

allscores <- list()

for (numFac in c(1:10)) {
        ## change number of extracted factors
        FAScores <- function(abundance.matrix){
                # computes pca.scores for first two PC´s from relative abundance matrix
                scent.pca <- as.data.frame(cmdscale(dist(abundance.matrix), k = numFac))
                #PCAScores <- as.data.frame(scent.pca$points)
               # PCAScores
                # scree plot --> 2 PC seem to be optimal
                # plot(scent.pca,type="lines")
        }
        
        source("ScentResultsNoScores.R")
        res <- ScentResultsNoScores("all","fa",1)
        
        scores <- res[[7]]
        factors <- res[[8]]
        
        ## list of relatedness results
        allscores <- append(allscores, list(scores))  
        
}


## small function to extract best model and deviance explained
devExpl <- function(df) {
        df <- as.data.frame(df)
        ## get deviance from het model
        model <- MinMod(df)
        model <- model[[2]]
        # get deviance explained --> (Null dev - Res dev) / Null dev
        devExpl <- (model$null.deviance - model$deviance)/model$null.deviance
        devExpl
}

##### get deviance explained for beach models ######
## change df format to have beach on column 1
changeCols <- function(df) {
        dfnew <- df
        dfnew[, 1] <- df[, ncol(df)]
        dfnew[, 2:ncol(dfnew)] <- df[, 1:(ncol(df)-1)]
        names(dfnew)[1] <- names(df)[ncol(df)]
        names(dfnew)[2:ncol(dfnew)] <- names(df)[1:(ncol(df)-1)]
        dfnew
}
newallscores <- lapply(allscores, changeCols)

## new function for binomial
devExpl2 <- function(df) {
        df <- as.data.frame(df)
        ## get deviance from het model
        model <- ParsModel(df, family = binomial)
        # get deviance explained --> (Null dev - Res dev) / Null dev
        devExpl <- (model$null.deviance - model$deviance)/model$null.deviance
        devExpl
}

beachDevExpl <- sapply(newallscores, devExpl2) # alternative:lapply

##### get deviance explained for mum pup models ######
## get family instead of beach

changetoMP <- function(df, factor) {
        df[, 1] <- factor
        df
}

allscoresfamily <- lapply(newallscores, changetoMP, factor=factors$Family)

MumPupDevExpl <- sapply(allscoresfamily, devExpl) # alternative:lapply



## construct ALLDATA frame
devExplained <- data.frame("HetMums" <- hetDevExpl, "RelMums" <- allMantel, "Beach" <- beachDevExpl, "MumPup" <- MumPupDevExpl)
names(devExplained) <- c("HetMums","RelMums","Beach","MumPup" )
devExplained$RelMums <- abs(devExplained$RelMums)

##add BIC
bic <- crits$vss.stats$BIC
devExplained$BIC <- bic

#library(xlsx)
#write.xlsx(devExplained, "Deviance Explained All Models.xlsx")

ggplot(devExplained, aes(x = 1:10)) +
        geom_line(aes(y = HetMums)) +
        geom_line(aes(y = RelMums)) +
        geom_line(aes(y = Beach)) +
        geom_line(aes(y = MumPup)) +
        theme_bw()
