# extracting moms and pubs matrices

## for relatedness
# extract Mums from Relate Mat by substracting all rows beginning with P
rel.mum <- relate[- grep("^P", rownames(relate)), - grep("^P", colnames(relate))]
rel.pup <- relate[- grep("^M", rownames(relate)), - grep("^M", colnames(relate))]

## for scent conc
scent.mum <- scent[ - grep("^P", rownames(scent)), - grep("^P", colnames(scent))]
scent.pup <- scent[ - grep("^M", rownames(scent)), - grep("^M", colnames(scent))]

## for scent abs/pres
scent.pa.mum <- scent[ - grep("^P", rownames(scent.pa)), - grep("^P", colnames(scent.pa))]
scent.pa.pup <- scent[ - grep("^M", rownames(scent.pa)), - grep("^M", colnames(scent.pa))]

# for Primer: exchanging NA with empty #done in excel
#relMum[is.na(relMum)] <- 
#relPub[is.na(relPub)] <- 
#scentMum[is.na(scentMum)] <- 
#scentPub[is.na(scentPub)] <- 

library(dataframes2xls)
write.xls(x = rel.mum, file = "relatednessMum.xls", sh.names="Mum",row.names = T)
write.xls(x = rel.pup, file = "relatednessPub.xls", sh.names="Pub",row.names = T)

write.xls(x = scent.mum, file = "scentMum.xls", sh.names="Mum",row.names = T)
write.xls(x = scent.pup, file = "scentPub.xls", sh.names="Pub",row.names = T)

write.xls(x = scent.pa.mum, file = "scentPaMum.xls", row.names=T)
write.xls(x = scent.pa.pup, file = "scentPaPup.xls", row.names=T)


