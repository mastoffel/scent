# calculating relative average concentration and variation of concentration for 
# important compounds
library(grid)
library(ggplot2)
theme.paper.scores <- theme_minimal() +
        theme(strip.text.x = element_text(vjust=1,size = 18),
              #strip.background = element_rect(fill="white",linetype="blank"),
              axis.title.x = element_text(vjust=0,size = 24),
              axis.title.y = element_text(vjust=0,size = 24),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20),
              plot.margin = unit(c(1,1,1,1), "cm"), # less margin on down side
              # axis.line = element_line(size=1),
              # panel.border = element_blank() ,
              panel.grid.major = element_blank(),
              # axis.line = element_line(color = 'black'),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour="black",fill=NA,size=1))


scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "projects\\sealscent\\data_files\\",
                                                  "Rdata\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))

simper_mp_ind <- c(58,60,68,74,86,90,96,107,164,181,189,209)
comp_ind_m <- c(36,52,86,88,96,103,110,203,206,207)
simper_beach_ind <- c(58,62,68,74,86,89,90,98,106,107,110,164,181,189,211)


length(intersect(simper_mp_ind, simper_beach_ind)) # 9
length(intersect(simper_mp_ind, comp_ind_m)) # 2
length(intersect(simper_beach_ind, comp_ind_m))

# extract mean concentration from samples where substances are non-zero
getmean <- function(x) {
        #out <- sd(x[x!=0])
        out <- mean(x)
}

getsd <- function(x) {
        #out <- sd(x[x!=0])
        out <- sd(x)
}

mp_conc <- sapply(scent_abundance[, simper_mp_ind], getmean)
beach_conc <- sapply(scent_abundance[, simper_beach_ind], getmean)
relate_conc <- sapply(scent_abundance[, comp_ind_m],  getmean)


mp <- data.frame("mp" = mp_conc)
beach <- data.frame("beach" = beach_conc)
rel <- data.frame("rel" = relate_conc)

# preparing df for boxplot
alldf <- matrix(rep(NA, 74), ncol = 2)
alldf[1:15] <- beach_conc
alldf[16:27] <- mp_conc
alldf[28:37] <- relate_conc

alldf[, 2] <- c(rep(1, 15), rep(2, 12), rep(3, 10))
df <- as.data.frame(alldf)
names(df) <- c("conc", "substances")
df$substances <- factor(df$substances, labels = c("beach", "mumpup", "relatedness"))

# fill = substances
ggplot(df, aes(x = substances, y = conc)) +
        geom_boxplot(colour = "black", lwd = 1) +
        theme.paper.scores +
        # coord_flip() +
        ylab("sd relative conc") +
        guides(fill = FALSE) +
        theme(plot.margin = unit(c(2,2,2,2), "cm"),
               axis.title.y = element_text(vjust = 1.7),
               axis.title.x = element_text(vjust = 0.02))
