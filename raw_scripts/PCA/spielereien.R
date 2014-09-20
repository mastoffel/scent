library("lattice")
NumComps <- het.mums.df$NumComps
het <- het.mums.df$het
t <- xyplot(NumComps ~ het, panel=function(x,y,...) {
        panel.xyplot(het, NumComps, ...)
        panel.lmline(het,NumComps) 
})
t


cutpoints <- quantile(relate.mums.df$genetic.distance, seq(0,1,length=4),na.rm=TRUE)
relate.mums.df$fac <- cut(relate.mums.df$genetic.distance,cutpoints)

g <- ggplot(relate.mums.df,aes(genetic.distance,pca.diff)) +
        geom_point(alpha=1/3) +
        #facet_wrap( ~fac, nrow=3) +
        geom_smooth(method="lm")

