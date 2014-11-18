# shannon

scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "projects\\sealscent\\data_files\\",
                                                  "Rdata\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))

scent_diversity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                  "projects\\sealscent\\data_files\\",
                                  "Rdata\\csv_files\\",
                                  "scent diversity.csv", sep = ""),
                                row.names=1)

heterozygosity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                 "projects\\sealscent\\data_files\\",
                                 "Rdata\\csv_files\\",
                                 "heterozygosity_41loci.csv", sep = ""),
                                row.names=1) 

scent_abund <- read.csv("scent_abundance_untransformed.csv", row.names = 1)
scent_abund <- as.data.frame(t(scent_abund))

#
library(vegan)
?diversity
div2 <- diversity(scent_abund, MARGIN = 1, "simpson")
het_df <- cbind(heterozygosity, div, div2)
names(het_df) <- c("het", "Shannon", "InvSimpson")
het_df <- het_df[1:41, ]


mod_shan <- lm(het_df$Shannon ~ het_df$het) 
mod_simp <- lm(het_df$InvSimpson ~ het_df$het) 


div_mum <- div[1:41]
het_mum <- heterozygosity[1:41,]
mod <- lm(div_mum ~ het_mum)
plot(div_mum ~ het_mum)
abline(mod)

fm_S <- lm(div_mum$S ~ het_mum)
fm_d <- lm(div_mum$d ~ het_mum)
fm_Shannon <- lm(div_mum$Shannon ~ het_mum)
fm_Simpson <- lm(div_mum$Simpson ~ het_mum)

plot(fm_Shannon)

plot(div_mum$Shannon ~ het_mum)
abline(fm_Shannon)
summary(fm_Shannon)
summary(fm_Simpson)

# calculating diversity with R



