## getting results from simper analysis, figure out which substances 
## were at least 10 times in the best contributing 90% to similarity between
## mum and pubs and get their indices for primer

SIMPER <- read.csv("SIMPER Results.csv", header=T) 
RT <- read.csv("AllRT.csv",header=F)
rt <- SIMPER$Volatile

RT <- as.vector(RT)
hist(table(rt))

best <- sort(table(rt))
best.over10 <- best[best>=10]
newvec <- as.numeric(names(best.over10))
write.xls(x = newvec, file = "mostcontributed.xls")

allind <- numeric()
for (i in seq_along(newvec)) {
  for (k in 1:nrow(RT)) {
    if (newvec[i]==RT[k,1]){
      allind <- append(allind, k)
    }
  }
}

write.xls(allind,"simperbest.xls") 