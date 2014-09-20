# further analyses

# sum of explained variance per element
sum.var <- apply(final.data,1,sum,na.rm=TRUE)
plot(sum.var)
axis(1, at=1:12,labels = names(sum.var))
