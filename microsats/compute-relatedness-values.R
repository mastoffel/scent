library(Demerelate)
gen <- read.table("C://Users//Martin//Studium//MSc.Behaviour//Research//Seal Scent//R code//raw scripts//41_loci.txt", header = TRUE)
gentest <- gen[, 1:8]

rel <- emp.calc(gentest, value="rxy", ref.pop=gentest)


vals <- rels$Empirical_Relatedness

relate <- relnew$BAS

data(demerelpop)
