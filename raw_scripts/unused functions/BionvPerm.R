## implementing BVStep algorithm from Primer
require(vegan)
source("bio.env.R")
source("bv.step.R")
# already standardized and transformed, transposed abundance matrix
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 

# relatedness matrix, fill up upper tri matrix and diagonal
relatedness <- read.csv(".\\csv_files\\relatedness_41loci.csv",row.names=1)
relInv <- t(relatedness)
relInv[lower.tri(relInv,diag=TRUE)] <- 1
relatedness[upper.tri(relatedness,diag=TRUE)] <- relInv[upper.tri(relInv,diag=TRUE)]

# create empty data frame
allcomps <- data.frame(matrix(nrow=ncol(scent.abundance),ncol=1))
row.names(allcomps) <- names(scent.abundance)
# zero occurences of all compounds
names(allcomps) <- "num"
allcomps[,1] <- 0

# to sample from
samp <- c(1:nrow(scent.abundance))

for (i in c(1:10)) {
        
        indiv <- sample(samp, 41)
        scent <- as.matrix(scent.abundance[indiv, ])
        
        relate <- as.matrix(relatedness[indiv,indiv])
        
        # get bioenv
        res <- bv.step(relate,scent, fix.dist.method="bray", var.dist.method="bray",
                       scale.fix=FALSE, scale.var=FALSE,
                       max.rho=0.95,
                       min.delta.rho=0.001,
                       random.selection=TRUE,
                       prop.selected.var=0.2,
                       num.restarts=10,
                       var.always.include=NULL,
                       var.exclude=NULL,
                       output.best=10) 
        
        # get best compounds
        comps <- (strsplit(res$best.model.vars,","))
        comps <- comps[[1]]
        
        # add compound occurence to results
        allcomps[comps,1] <- allcomps[comps,1] + 1

}

allcomps$names <- row.names(allcomps)
sorted.allcomps <- arrange(allcomps, -num)