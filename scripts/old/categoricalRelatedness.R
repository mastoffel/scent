library(plyr)
source("summarySE.R")

genclass <- cut(relate.df$relatedness,8)
relate.df$genclass <- as.factor(genclass)
levels(relate.df$genclass) <- c("1","2","3","4","5",
                                "6","7","8")

group1 <- relate.df[relate.df$genclass==1, ]

newdata <- relate.df[order(relate.df$relatedness), ]
equgroups <- rep(1:10, rep(82,10))
relate.df$groups <- equgroups

t <- summarySE(relate.df, measurevar = "F1", groupvars = "groups")


     

boxplot(relate.df$pca.diff ~ relate.df$genclass)

barplot(t$mean)

ggplot(t, aes(x = 1:10,y = F1)) +
        theme_bw(base_size=20) +
        geom_point( fill = "black", size=5) +
        geom_errorbar(aes(ymin=F1-se, ymax=F1+se), width=.1) +
        scale_x_continuous(breaks=c(seq(from=1, to = 8, by =1))) 
        