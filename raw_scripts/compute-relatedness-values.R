library(Demerelate)
gen <- read.table("C://Users//Martin//Studium//MSc.Behaviour//Research//Seal Scent//R code//raw scripts//41_loci.txt", header = TRUE)


rels2 <- Demerelate(gen, object=TRUE, value="rxy" , NA.rm=FALSE, ref.pop = "NA", file.output=FALSE, pairs = 1000)
rel <- emp.calc(gen, value="rxy", ref.pop=gen)


vals <- rels$Empirical_Relatedness

relate <- relnew$BAS

data(demerelpop)
