## Venn Diagrams
library("VennDiagram")
library("venneuler")

comp_ind_m  <- comp_ind_m[1:9]
# create simper venns for beach comparison

venn.diagram(list("beach" = simper_beach_ind, "mum pup" = simper_mp_ind, "best mum" = comp_ind_m),
             fill=c("red","green", "blue"), filename="VennSimper.tiff",
             sub.cex = 30)

#computing indices of elements for simper best five - three compounds here: beach 1, beach 2, beach dissimilarity
beach.1.2 <- intersect(beach.1,beach.2) # overlap between beaches
just.between.beaches <- !(intersect(beach.1,beach.2)%in% beach.dis ) #this is a logical
just.beach.1.2 <- beach.1.2[just.between.beaches] # overlap between beaches minus elements accounting for beach dissimilarity
beach.1.2.dis <- intersect(intersect(beach.1,beach.2),beach.dis) # overlap between beach 1, beach 2 and beach dissim
