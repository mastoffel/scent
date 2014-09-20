# sum of explained variance per pair through simper elements
groups <- simper.list[[2]]

sum.var <- apply(best.simper.one,2,sum,na.rm=TRUE)
plot(sum.var, xlab = "group ID", ylab = "Mum-Pup-Similarity explained through elements from Simper")
lines(sum.var)
text(sum.var,groups, cex=0.6, pos = 2, col = "red")

sum.var <- apply(best.simper.two,2,sum,na.rm=TRUE)
plot(sum.var, xlab = "group ID", ylab = "Mum-Pup-Similarity explained through elements from Simper")
lines(sum.var)
text(sum.var,groups, cex=0.6, pos = 2, col = "red")


sum.var <- apply(best.simper.three,2,sum,na.rm=TRUE)
plot(sum.var, xlab = "group ID", ylab = "Mum-Pup-Similarity explained through elements from Simper")
lines(sum.var)
text(sum.var,groups, cex=0.6, pos = 2, col = "red")

sum.var <- apply(best.simper.four,2,sum,na.rm=TRUE)
plot(sum.var, xlab = "group ID", ylab = "Mum-Pup-Similarity explained through elements from Simper")
lines(sum.var)
text(sum.var,groups, cex=0.6, pos = 2, col = "red")

sum.var <- apply(best.simper.five,2,sum,na.rm=TRUE)
plot(sum.var, xlab = "group ID", ylab = "Mum-Pup-Similarity explained through elements from Simper")
lines(sum.var)
text(sum.var,groups, cex=0.6, pos = 3, col = "red")
#


# plot explained variance per element

sum.var <- apply(best.simper.one,1,sum,na.rm=TRUE)
plot(sum.var, xlab = "element", ylab = "contribution")
lines(sum.var)
element.names <- as.character(round(as.numeric(as.vector(rownames(best.simper.one))),2))
text(sum.var,element.names, cex=1, pos = 2, col = "red")

sum.var <- apply(best.simper.two,1,sum,na.rm=TRUE)
plot(sum.var, xlab = "element", ylab = "contribution")
lines(sum.var)
element.names <- as.character(round(as.numeric(as.vector(rownames(best.simper.two))),2))
text(sum.var,element.names, cex=0.8, pos = 3, col = "red")

sum.var <- apply(best.simper.three,1,sum,na.rm=TRUE)
plot(sum.var, xlab = "element", ylab = "contribution")
lines(sum.var)
element.names <- as.character(round(as.numeric(as.vector(rownames(best.simper.three))),2))
text(sum.var,element.names, cex=1, pos = 2, col = "red")

sum.var <- apply(best.simper.four,1,sum,na.rm=TRUE)
plot(sum.var, xlab = "element", ylab = "contribution")
lines(sum.var)
element.names <- as.character(round(as.numeric(as.vector(rownames(best.simper.four))),2))
text(sum.var,element.names, cex=1, pos = 3, col = "red")

sum.var <- apply(best.simper.five,1,sum,na.rm=TRUE)
plot(sum.var, xlab = "element", ylab = "contribution")
lines(sum.var)
element.names <- as.character(round(as.numeric(as.vector(rownames(best.simper.five))),2))
text(sum.var,element.names, cex=0.7, pos = 3, col = "red")