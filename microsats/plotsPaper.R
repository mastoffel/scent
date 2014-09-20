# seal scent main plots (needs script "Seal Scent Main.R)
library(grid)
# create theme for publication
theme.paper <-  theme_bw() +
                theme(strip.text.x = element_text(vjust=1,size = 18),
                      strip.background = element_rect(fill="white",linetype="blank"),
                      axis.title.x = element_text(vjust=0.1,size = 16),
                      axis.title.y = element_text(vjust=0.1,size = 16),
                      axis.text.x = element_text(size = 14),
                      axis.text.y = element_text(size = 14),
                      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"), # less margin on down side
                      # axis.line = element_line(size=1),
                      # panel.border = element_blank() ,
                      panel.grid.major = element_blank(),
                      # axis.line = element_line(color = 'black'),
                      panel.grid.minor = element_blank(),
                      panel.border = element_rect(colour="black",fill=NA))
        

# Heterozygosity plots Mums
# heterozygosity vs. number of compounds
het.all.plot <- ggplot(het.mums.df, aes(x=het,y=NumComps)) +
        geom_point(colour = "black", size = 2) +
        geom_smooth(method="lm",size = 0.7 ,alpha=0.2, colour="blue") +
        theme.paper +
        theme(axis.text.x = element_blank()) +
        scale_x_continuous(breaks=c(seq(from = 0.8, to = 1.20, by = 0.1)))+
        #geom_text(aes(0.85,80, label="(a) r = 0.34, p = 0.027"),size=4) +
        xlab("") +
        ylab("number of compounds") 
ggsave(file="hetall.svg", width=5, height=4)

# heterozygosity vs. PC2
het.pc2.plot <- ggplot(het.mums.df,aes(x=het,y=PC2)) +
        geom_point(colour = "black", size = 2) +
        geom_smooth(method="lm" ,size = 0.7 ,alpha=0.2,colour="blue") + #fill = "lightblue"        
        theme.paper +    
        scale_x_continuous(breaks=c(seq(from = 0.8, to = 1.20, by = 0.1)))+
        xlab("heterozygosity") +
        ylab("PC2 score") 
ggsave(file="hetpc.svg", width=5, height=4)

#all.het <- multiplot(het.all.plot,het.pc2.plot,cols=1)


# relatedness

rel.pc <- ggplot(relate.mums.df, aes(x = genetic.distance, y = pca.diff)) +
        geom_point(colour="black", size = 2.5,alpha = 0.6 ,show_guide=FALSE) +
        geom_smooth(method="lm",size=1.1,colour="blue", fill="grey",alpha=0.3) +
        #geom_point(size = 3,colour="grey", alpha=0.2, show_guide = F) +          
        theme.paper +
        geom_text(aes(0.72,4.5, label="Mantel R = 0.09**"),size=4) +
        xlab("genetic distance") +
        labs(y=expression("D"[PC6]))

# pc.scores density

pc.dens.plot.1 <- ggplot(pca.scores,aes(x=PC1)) +
        geom_density(alpha=.5, size=0.5, aes(fill = beach)) +
        scale_fill_manual(values = c("green","blue")) +
        guides(fill=guide_legend(title=NULL)) +
        theme.paper + 
        xlab("PC1 Score") +
        ylab("Density")
pc.dens.plot.2 <- ggplot(pca.scores,aes(x=PC2)) +
        geom_density(alpha=.5, size=0.5, aes(fill = beach)) +
        scale_fill_manual(values = c("green","blue")) +
        guides(fill=guide_legend(title=NULL,direction="horizontal")) +
        theme.paper + 
        xlab("PC2 Score") +
        ylab("Density")

multiplot(pc.dens.plot.1,pc.dens.plot.2,cols=2)

## heterozygosity vs. PC

# get PCA screeplot explained variances
scent.pca <- prcomp(scent.abundance)
screeplot(scent.pca,type = "line",main="Scree Plot", xlab="Components", ylab="Eigenvalue")
Eigen <- scent.pca$sd^2
Eigen <- Eigen[1:10]
exp.var <- scent.pca$sdev^2 / sum(scent.pca$sdev^2)
exp.var.10 <- exp.var[1:10]



# create df
pc.plot.df <- data.frame("screeplot"=exp.var.10,
                         "mums"=results.het.mums$rsquared,
                         "pups"=results.het.pups$rsquared,
                         "PC"=as.integer(1:10))
row.names(pc.plot.df) <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

# plot het vs. pc

ggplot(pc.plot.df,aes(x=PC,y="mums heterozygosity")) +
        geom_line(aes(y = screeplot, colour = "screeplot"), size=1.5) +
        geom_line(aes(y = mums, colour = "mums heterozygosity"),size=1.5) +
        geom_line(aes(y = pups, colour = "pups heterozygosity"),size=1.5) +
        theme_bw(base_size=20) +
        theme(strip.text.x = element_text(vjust=1,size = 14),
              legend.position=c(0.80,0.85),
              strip.background = element_rect(fill="white",linetype="blank"),
              axis.title.x = element_text(vjust=0.1,size = 14),
              axis.title.y = element_text(vjust=0.2,size = 14),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour="black",fill=NA),
              legend.title=element_blank()) +
        scale_x_continuous(breaks = c(seq(from = 1, to = 10, by = 1))) +
        xlab("principal components") +
        ylab("explained variance")

## relatedness vs. pc


pc.plot.df.rel <- data.frame("screeplot"=exp.var.10,
                         "mums"=results.relatedness.mums$mantelR,
                         "pups"=results.relatedness.pups$mantelR,
                         "PC"=as.integer(1:10))
row.names(pc.plot.df.rel) <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

## plot het vs. relatedness
ggplot(pc.plot.df.rel,aes(x=PC,y="mums")) +
        ggtitle("PC´s vs. genetic distance") +
        #geom_line(aes(y = screeplot, colour = "screeplot"), size=1.5) +
        geom_line(aes(y = mums, colour = "mums"),size=1.5) +
        geom_line(aes(y = pups, colour = "pups"),size=1.5) +
        geom_abline(intercept=0.068,slope=0, linetype='dashed') +
        theme.paper +
        theme(legend.title=element_blank()) +
        scale_x_continuous(breaks = c(seq(from = 1, to = 10, by = 1))) +
        xlab("principal components") +
        ylab("Mantel R")
        

# very nice plot: Variable loadings vs. PC
load    <- scent.pca$rotation
sorted.loadings <- load[order(load[, 1]), 1] # change both numbers for PC change
myTitle <- "Loadings Plot for PC1" 
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")
