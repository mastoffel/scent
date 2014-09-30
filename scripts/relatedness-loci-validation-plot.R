# loci validation relatedness plot

library(ggplot2)
library(grid)
source("SumResults.R")

rel_per_loc <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                  "projects\\sealscent\\data_files\\",
                                  "Rdata\\csv_files\\",
                                  "locirelate.csv", sep = ""),
                                  row.names=1)

# get mean, sd and se per loci-number in table
loc_sum <- SumResults(rel_per_loc)

ggplot(loc_sum, aes(x = locnum, y = cormean)) +
        geom_errorbar(aes(ymin = cormean-se, ymax = cormean+se),
                      width=0.8, alpha=0.7, size = 0.8) +
        geom_point(size = 3, shape = 16) +
        geom_line(size = 0.8) +
        theme_minimal(base_size = 24) +
        theme(
                axis.title.x = element_text(vjust= -2 ,size = 26),
                axis.title.y = element_text(vjust=3,size = 26),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                plot.margin = (unit(c(.5, .5, 2, 2), "cm"))) +
        #geom_hline(yintercept=0.305, linetype = "") +
        ylab(expression(mantelsR [mean+-se])) +
        xlab("number of loci") +
        theme(legend.title=element_blank(),
              legend.position="none") +
        scale_colour_brewer(type="qual",palette="Set2") +
        scale_shape_manual(values=c(1,19,0,15)) 