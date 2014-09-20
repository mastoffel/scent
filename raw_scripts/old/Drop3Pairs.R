# compute new matrices

# drop the three pairs which are not in relatedness data

drops <- c("M1","P1","M22","P22","M32","P32")

relate.drop <- relate[!(rownames(relate) %in%drops),!(colnames(relate) %in%drops)]
scent.drop <- scent[!(rownames(scent) %in%drops),!(colnames(scent) %in%drops)]

write.xls(x = relate.drop, file = "relate41.xls", sh.names="relate41",row.names = T)
write.xls(x = scent.drop, file = "scent41.xls", sh.names="scent41",row.names = T)