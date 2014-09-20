# seal scent main plots (needs script "Seal Scent Main.R)
 library(ggthemes)
 library(grid)
# create theme for publication
 # width of plotarea ( without axes labels and stuff) : 117.415, height: 67.261
 # ggsave should be (width = 151.585, height = 111.539)
theme.paper <-  theme_minimal() +
                theme(strip.text.x = element_text(vjust=1,size = 18),
                      #strip.background = element_rect(fill="white",linetype="blank"),
                      axis.title.x = element_text(vjust=0,size = 16),
                      axis.title.y = element_text(vjust=0,size = 16),
                      axis.text.x = element_text(size = 14),
                      axis.text.y = element_text(size = 14),
                      plot.margin = unit(c(0,0,1,1), "cm"), # less margin on down side
                      # axis.line = element_line(size=1),
                      # panel.border = element_blank() ,
                      panel.grid.major = element_blank(),
                      # axis.line = element_line(color = 'black'),
                      panel.grid.minor = element_blank(),
                      panel.border = element_rect(colour="black",fill=NA,size=1))
 
theme.paper.scores <- theme_minimal() +
         theme(strip.text.x = element_text(vjust=1,size = 18),
               #strip.background = element_rect(fill="white",linetype="blank"),
               axis.title.x = element_text(vjust=0,size = 24),
               axis.title.y = element_text(vjust=0,size = 24),
               axis.text.x = element_text(size = 20),
               axis.text.y = element_text(size = 20),
               plot.margin = unit(c(0,0,1,1), "cm"), # less margin on down side
               # axis.line = element_line(size=1),
               # panel.border = element_blank() ,
               panel.grid.major = element_blank(),
               # axis.line = element_line(color = 'black'),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(colour="black",fill=NA,size=1))
 
 
# Heterozygosity plots Mums
# heterozygosity vs. number of compounds
het.all.plot <- ggplot(het.mums.df, aes(x=het,y=NumComps)) +
        geom_point(colour = "black", size = 2.5) +
        geom_smooth(method="lm",size = 0.7 ,alpha=0.15, colour="black") +
        theme.paper +
        #theme(axis.text.x = element_blank()) +
        scale_x_continuous(breaks=c(seq(from = 0.8, to = 1.20, by = 0.1))) +
        #geom_text(aes(0.85,80, label="(a) r = 0.34, p = 0.027"),size=4) +
        xlab("") +
        ylab("") 
ggsave(file="hetall.pdf", width=151.585, height=111.539, units=c("mm"),useDingbats=FALSE)

# heterozygosity vs. F1F2
het.F1F2.plot <- ggplot(het.df,aes(x=het,y=F1F2)) +
        geom_point(colour = "black", size = 2.5) +
        geom_smooth(method="lm" ,size = 0.7 ,alpha=0.15,colour="black") + #fill = "lightblue"        
        theme.paper +    
        scale_x_continuous(breaks=c(seq(from = 0.8, to = 1.20, by = 0.1))) +
        xlab("") +
        ylab("") 
ggsave(file="hetF1F2.pdf", width=151.585, height=111.539, units=c("mm"),useDingbats=FALSE)
#all.het <- multiplot(het.all.plot,het.pc2.plot,cols=1)


# relatedness

rel.pc <- ggplot(relate.df, aes(x = F3, y = relatedness)) +
        geom_point(colour="black", size = 2.5,alpha = 0.6 ,show_guide=FALSE) +
        geom_smooth(method="lm",size=1.1,colour="blue", fill="grey",alpha=0.3) +
        #geom_point(size = 3,colour="grey", alpha=0.2, show_guide = F) +          
        theme.paper +
        #geom_text(aes(0.72,4.5, label="Mantel R = 0.09**"),size=4) +
        ylab("relatedness") +
        labs(x=expression("Diff"[F1]))

# fa.scores density
scores$beach <- as.factor(factors$Beach)
fa.dens.plot.1 <- ggplot(scores,aes(x = F1)) +
        geom_density(alpha=0.8, size=0.5, aes(fill = beach),adjust=1.5) +
        scale_fill_manual(values = c("blue","red")) +
        guides(fill=guide_legend(title=NULL)) +
        theme.paper.scores + 
        theme(legend.position="none") +
        scale_x_continuous(breaks = c(seq(from = -1, to = 6, by = 1))) + 
        scale_y_continuous(breaks = c(seq(from = 0, to = 1.4, by = 0.4))) +  
        #scale_y_continuous(breaks = c(seq(from = 0, to = 1, by = 0.2))) +
        xlab("Factor 1") +
        ylab("Density")
 ggsave(file="F1 scores.pdf", width=151.585, height=111.539, units=c("mm"),useDingbats=FALSE)
 
 fa.dens.plot.2 <- ggplot(scores,aes(x=F2)) +
         geom_density(alpha=0.8, size=0.5, aes(fill = beach),adjust=1.5) +
         scale_fill_manual(values = c("blue","red")) +
         guides(fill=guide_legend(title=NULL)) +
         theme.paper.scores + 
         theme(legend.position="none") +
         scale_y_continuous(breaks = c(seq(from = 0, to = 0.8, by = 0.2))) +
         xlab("Factor 2") +
         ylab("Density")
 ggsave(file="F2 scores.pdf", width=151.585, height=111.539, units=c("mm"),useDingbats=FALSE)
 
 fa.dens.plot.3 <- ggplot(scores,aes(x=F3)) +
         geom_density(alpha=0.8, size=0.5, aes(fill = beach),adjust=1.3) +
         scale_fill_manual(values = c("blue","red")) +
         guides(fill=guide_legend(title=NULL)) +
         theme.paper.scores + 
         theme(legend.position="none") +
         scale_y_continuous(breaks = c(seq(from = 0, to = 2, by = 0.4))) +
         scale_x_continuous(breaks = c(seq(from = -5, to = 1, by = 1))) + 
         xlab("Factor 3") +
         ylab("Density")
 ggsave(file="F3 scores.pdf", width=151.585, height=111.539, units=c("mm"),useDingbats=FALSE)
 
 
fa.dens.plot.4 <- ggplot(scores,aes(x=F4)) +
        geom_density(alpha=0.8, size=0.5, aes(fill = beach),adjust=1.5) +
        guides(fill=guide_legend(title=NULL,direction="horizontal",colour = guide_legend(override.aes = list(shape=NA)))) +
        scale_fill_manual(values = c("blue","red"),labels=c("Special Study Beach  ","Freshwater Beach")) +
        theme.paper.scores  + 
        theme(legend.position = c(0.70,0.8), 
              legend.text = element_text(size = 24)) +
        scale_y_continuous(breaks = c(seq(from = 0, to = 1.2, by = 0.2))) +
        scale_x_continuous(breaks = c(seq(from = -1, to = 7, by = 1))) + 
        xlab("Factor 4") +
        ylab("Density")
 ggsave(file="F4 scores.pdf", width=151.585, height=111.539, units=c("mm"),useDingbats=FALSE)
 
 
 
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
scent.pca <- prcomp(scent.abundance) 
load    <- scent.pca$rotation
sorted.loadings <- load[order(load[, 2]), 2] # change both numbers for PC change
#sorted.loadings.1 <- load[order(load[, 1]), 1]
myTitle <- "Loadings Plot for PC1" 
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")
 
 
 ## try nicer scores plot
 sco <- scores[,-6]
 sco <- as.vector(as.matrix(sco))
 sco <- data.frame(Fac=sco)
 sco$beach <- rep(scores$beach,5)
 sco$Num <- c(rep(1,82),rep(2,82),rep(3,82),rep(4,82),rep(5,82))
 sco$Num <- as.factor(sco$Num)
 fa.dens.plot<- ggplot(sco,aes(x = Fac)) +
         geom_density(alpha=0.8, size=0.5, aes(fill = beach),adjust=3) +
         scale_fill_manual(values = c("blue","red")) +
         guides(fill=guide_legend(title=NULL)) +
         theme.paper + 
         theme(legend.position="none") +
         xlab("F1 score") +
         ylab("Density")
 
 
 ## 3d plot heterozygosity
 library(lattice)
 library(grid)
 library(locfit)
 library(xlsx)
 ## alternative kernel weight function: "gauss"/"geom", crucial: smoothing pa alpha
 attach(het.df)
 plot(locfit(het ~ F1 + F2, data=het.df, maxk=200, family="gaussian", kern="tcub", alpha = 0.98),
      type = "persp", theta=5,phi=25) 
 data <- locfit(het ~ F1+F2, data=het.df, maxk=200, family="geom", kern="tcub", alpha = 0.98)
 locfitdata <- expand.grid(F1 = seq(-0.860,5.4,0.08), F2=seq(-0.86,4.5, 0.08)) ## basically the F1 F2 ranges
 locfitdata$het <- predict(data, newdata=locfitdata)
 write.xlsx(locfitdata,"locfitdata.xlsx")
 
 