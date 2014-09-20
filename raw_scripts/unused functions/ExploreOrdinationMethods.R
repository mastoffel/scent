## script to explore different ordination methods

setwd("C:/Users/Martin/Studium/MSc.Behaviour/Research/Seal Scent/R code/Raw scripts")

# already standardized and transformed, transposed abundance matrix
scent.abundance <- as.data.frame(t(read.csv(".\\csv_files\\scent abundances.csv",row.names=1))) 

# 6 different diversity indices (from primer)
scent.diversity <- read.csv(".\\csv_files\\scent diversity.csv",row.names=1) 

# relatedness matrix
relatedness <- read.csv(".\\csv_files\\relatedness_41loci.csv",row.names=1)

# heterozygosity vector SH
heterozygosity <- read.csv(".\\csv_files\\heterozygosity_41loci.csv", row.names=1) 

# beach and family vector
factors <- read.csv(".\\csv_files\\factors.csv",row.names=1) 

heterozygosity= subset(factors, select=Family)

## exploration heterozygosity
source("DimExplore.R")
resultsFA <- DimExplore(heterozygosity,scent.abundance, method="fa", subset = "mums")
resultsPCA <- DimExplore(heterozygosity,scent.abundance, method="pca", subset = "mums")
resultsPCOA <- DimExplore(heterozygosity,scent.abundance, method="pcoa", subset = "mums")
resultsMDS <- DimExplore(heterozygosity,scent.abundance, method="mds", subset = "mums")

##exploration relatedness
source("DimExploreRel.R")
resultsFA <- DimExploreRel(relatedness,scent.abundance, method="fa", subset = "mums")
resultsPCA <- DimExploreRel(relatedness,scent.abundance, method="pca", subset = "mums")
resultsPCOA <- DimExploreRel(relatedness,scent.abundance, method="pcoa", subset = "mums")
resultsMDS <- DimExploreRel(relatedness,scent.abundance, method="mds", subset = "mums")


#allresults <- c(resultsFA,resultsPCA,resultsPCOA,resultsMDS)
#fac <- as.factor(c(rep(1,10),rep(2,10),rep(3,10),rep(4,10)))
#levels(fac) <- c("FA", "PCA", "PCOA", "MDS")
#hetdfall <- as.data.frame(cbind(allresults,fac))

hetdf <- as.data.frame(cbind(resultsFA,resultsPCA,resultsPCOA,resultsMDS))
cols <- c("FA"="blue","PCA"="red","PCOA"="green","MDS"="black")


theme.paper <-  theme_minimal() +
        theme(strip.text.x = element_text(vjust=1,size = 18),
              #strip.background = element_rect(fill="white",linetype="blank"),
              axis.title.x = element_text(vjust=0.2,size = 16),
              axis.title.y = element_text(vjust=2,size = 16),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              plot.margin = unit(c(1,1,1,1), "cm"), # less margin on down side
              # axis.line = element_line(size=1),
              # panel.border = element_blank() ,
              panel.grid.major = element_blank(),
              # axis.line = element_line(color = 'black'),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour="black",fill=NA,size=1))


ggplot(hetdf, aes(x = 1:10)) +
        geom_line(aes(y = resultsFA, color = "FA"), size = 1.0) +
        geom_line(aes(y = resultsPCA, color = "PCA"), size = 1) +
        geom_line(aes(y = resultsPCOA, color = "PCOA"), size = 1) +
        geom_line(aes(y = resultsMDS, color = "MDS"), size = 1) +
        scale_colour_manual(values=cols,name = "Method") +
        theme.paper +
        xlab("dimensions / factor number") +
        ylab("deviance explained") +
        scale_x_continuous(breaks = 1:10) +
        ggtitle("relatedness") +
        theme(plot.title = element_text(lineheight=0.5 , face="bold", size = 20))
    
        