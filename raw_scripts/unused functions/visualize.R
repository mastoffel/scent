visualize <- function(number_relate, number_scent, title, colour = "black") {
  ## Both inputs should be lower rectangular resemblance matrices 
  ## with labels in first column and first row
  ## first two args are meaning sheet number within excel file
  col <- colour
  ##libraries
  #library(ggplot2)
 # library(ecodist)
  #library(xlsx)
  ## load as data.frame
  scent <- read.xlsx("ALLRESEMBLANCE2.xlsx",number_scent,row.names=1)
  relate <- read.xlsx("ALLRESEMBLANCE2.xlsx",number_relate,row.names=1) 
  
  ## turn into vector
  scent.vec <- as.vector(as.matrix(scent))
  relate.vec <- as.vector(as.matrix(relate))
  
  ## delete NAs
  badscent <- is.na(scent.vec)
  badrelate <- is.na(relate.vec)
  
  scent.vec <- scent.vec[!badscent]
  relate.vec <- relate.vec[!badrelate]
  
  ## lm
  fit <- lm(relate.vec ~ scent.vec)
  
  ##mantel rank
  ranked.mantel  <- mantel(relate.vec ~ scent.vec, mrank = TRUE)
  mantel <- round(ranked.mantel[c("mantelr","pval1")],4)
  
  text  <- paste("RankMantel:","Rho =",mantel[1],"p =",mantel[2],sep=" ")
  
  ##plo
  title <- qplot(scent.vec,relate.vec) +
          geom_smooth(method = lm,size =1.5)+
          geom_text(x = 60, y = -0.1, aes(label = text)) +
          xlab("scent similarity") +
          ylab("relatedness") +
          ggtitle(title) +
          geom_point(colour = col) 
          theme_set(theme_bw(10))
  return(title)
}





