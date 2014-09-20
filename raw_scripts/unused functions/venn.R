## Venn Diagrams
library("VennDiagram")
library("venneuler")
## from best backwards, 41 pairs, C3C2, errorwindow 0.03

#allPA <- c(1,4,6,7,10,14,19,20,22,24,26,41,43,48,52,82,86,90,93,96,104,108,109,110,118,120,126,130,132,138,143,147,153,158,160,169,173,186,190,195,199,204,207,208,210,211)
#mumPA <- c(2,3,4,13,15,16,19,20,22,24,29,33,39,43,44,46,47,48,49,51,53,63,64,67,69,71,74,80,90,91,93,96,99,102,103,106,109,110,111,113,114,115,116,118,124,130,134,135,136,138,140,141,142,144-148,150,152,153,154,156-162,167-170,173,175-177,180,183,186,190,193,197,198,200,203,204,206-208,211,212,213)
#pupPA <- c(3,4,6,7,10,14,17,19,22,26,40,41,42,69,71,77,81,82,83,90,98,104,109,120,123,129,130,138,143,147,153,160,166,169,173,187,190,202,204,207)

#simperbest <- c(91,143,195,323,379,90,170,214,104,128,116,33,213,347,129,174,176,203,209,232,177,192,207,332,120,40,114,256,173,288,79,178)

#venn <- list(allcon,mumcon,pupcon)

## just necessary for package venneuler
#common.all.mom <- Reduce(intersect, list(allcon,mumcon))
#common.all.pup <- Reduce(intersect, list(allcon,pupcon))
#common.mum.pup <- Reduce(intersect, list(mumcon,pupcon))
#common.all.mum.pup <- Reduce(intersect, list(mumcon,pupcon,allcon))


# from VennDiagram
#venn.diagram(list(All = allcon, Mum = mumcon, Pup = pupcon),fill=c("red","green","blue"), filename="VennConc.tiff")
#venn.diagram(list(All = allPA, Mum = mumPA, Pup = pupPA),fill=c("purple","yellow","black"), filename="VennPA.tiff")
#venn.diagram(list(All = allcon, BetweenGroups = simperbest), fill = c("black","white"), filename="VennBestVsSimper.tiff")

# create simper venns for beach comparison

beach.dis <- c(58,68,74,86,90,106,107,164,181,211)
beach.1 <- c(8,13,42,43,58,60,62,74,83,86,90,96,98,106,107,108,110,123,164,181,189)
beach.2 <- c(1,8,13,43,52,58,59,60,62,68,74,80,89,90,98,100,104,110,164,181,189,206,209,211)

venn.diagram(list(Beach.Dissim = beach.dis, Beach.1 = beach.1, Beach.2 = beach.2),fill=c("red","green","blue"), filename="VennSimper5.tiff")

#computing indices of elements for simper best five - three compounds here: beach 1, beach 2, beach dissimilarity
beach.1.2 <- intersect(beach.1,beach.2) # overlap between beaches
just.between.beaches <- !(intersect(beach.1,beach.2)%in% beach.dis ) #this is a logical
just.beach.1.2 <- beach.1.2[just.between.beaches] # overlap between beaches minus elements accounting for beach dissimilarity
beach.1.2.dis <- intersect(intersect(beach.1,beach.2),beach.dis) # overlap between beach 1, beach 2 and beach dissim
