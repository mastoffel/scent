library(Demerelate)
gen <- read.table("C://Users//Martin//Studium//MSc.Behaviour//Research//Seal Scent//R code//raw scripts//41_loci.txt", header = TRUE)
gen <- gen[, -2]

gen2 <- input.txt("C://Users//Martin//Studium//MSc.Behaviour//Research//Seal Scent//R code//raw scripts//41_loci.txt", mod="pop")
gen2 <- gen2[, -2]
rels <- Demerelate(gen2[, 1:6], object=TRUE, value="rxy" , NA.rm=FALSE, ref.pop = "NA", file.output=TRUE)

rels <- Demerelate(inp)

vals <- rels$Empirical_Relatedness

relate <- relnew$BAS

data(demerelpop)
