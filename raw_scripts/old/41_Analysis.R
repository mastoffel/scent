## all operations from beginning with final matrices

# scent resemblance, relatedness and scent resemblance(PresAbs)
scent <- read.csv("scent41.csv", header=T, row.names=1) 
relate <- read.csv("relate41.csv", header=T, row.names=1) 
scent.pa <- read.csv("scent41_ap.csv", header=T, row.names=1)
scent.select <- read.csv("ResemBackwardDummy.csv", header=T, row.names=1)
#getting a mum pub factor --> computable with dist.object??
#mumpup <- read.csv("mumpubfactor.csv", header=T,row.names=1)
#f <- as.numeric(mumpub[1,])
#mumpup.factor <- factor(f)

scent.dist <- as.dist(scent)
relate.dist <- as.dist(relate)
scent.pa.dist <- as.dist(scent.pa)
scent.select.dist <- as.dist(scent.select)
#graphics
library(ggplot2)
library(lattice)
#xyplot(scent.dist ~ relate.dist | mumpup.factor)

#relevant packages
library(vegan)
library(ecodist)
library(ade4)

## mantel test with ecodist
library(ecodist)
library(vegan)
# for scent based on conc, and scent based on pres abs
# based on non rank
nonranked.mantel  <- mantel(relate.dist ~ scent.dist)     #pval1 = .003 ** !!
nonranked.mantel  <- mantel(relate.dist ~ scent.pa.dist)
nonranked.mantel <- mantel(relate.dist ~ scent.select.dist)
#based on rank --> pearson                                #pval = 0.173
ranked.mantel  <- mantel(relate.dist ~ scent.dist, mrank = TRUE)
ranked.mantel  <- mantel(relate.dist ~ scent.pa.dist, mrank = TRUE)
