library(plyr)
library(dplyr)
source("summarySE.R")

genclass <- cut(relate.df$relatedness,6)
relate.df$genclass <- as.factor(genclass)
levels(relate.df$genclass) <- c("1","2","3","4","5",
                                "6")

t <- summarySE(relate.df, measurevar = "F1", groupvars = "genclass")

ggplot(t, aes(x = 1:6,y = F1)) +
        theme_bw(base_size=20) +
        geom_point( fill = "black", size=5) +
        geom_errorbar(aes(ymin=F1-se, ymax=F1+se), width=.1) +
        scale_x_continuous(breaks=c(seq(from=1, to = 8, by =1))) 



group1 <- relate.df[relate.df$genclass==1, ]
group2 <- relate.df[relate.df$genclass==2, ]
group8 <- relate.df[relate.df$genclass==8, ]
group7 <- relate.df[relate.df$genclass==7, ]

newrelate <- relate.df[order(relate.df$relatedness), ]
equgroups <- rep(1:10, rep(82,10))
newrelate$groups <- equgroups

t <- summarySE(newrelate, measurevar = "F1", groupvars = "groups")

newr <- arrange(newrelate, groups, relatedness)



ggplot(t, aes(x = 1:10,y = F1)) +
        theme_bw(base_size=20) +
        geom_point( fill = "black", size=5) +
        geom_errorbar(aes(ymin=F1-se, ymax=F1+se), width=.1) +
        xlab("relatedness class") +
        ggtitle("equal sized groups") +
        scale_x_continuous(breaks=c(seq(from=1, to = 10, by =1))) 
       
start <- newr$relatedness[1]
stop1 <- vector()
for (i in 1:10) {
        if (i != 10) {
        start <- append(start, newr$relatedness[i * 82 + 1])
        }
        stop1 <- append(stop1, newr$relatedness[i * 82])
}
dfclass <- data.frame("start" = start, "stop" = stop1)

boxplot(relate.df$pca.diff ~ relate.df$genclass)

        