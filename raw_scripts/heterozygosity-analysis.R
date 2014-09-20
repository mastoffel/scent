# Load heterozygosity data (82*2 data.frame) and diversity data (82 * 14) data.frame
# construct plots and fit linear models with ggplot2


div <- read.csv(".\\csv_files\\Diversity.csv",row.names=1)
het <- read.csv(".\\csv_files\\heterozygosity2.csv",row.names=1)

library(ggplot2)
source("multiplot.R")
library(extrafont)
library(extrafontdb)
# Needed only on Windows - run once per R session
# Adjust the path to match your installation of Ghostscript
Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.14/bin/gswin64c.exe")

#devide into mums and pups
div_mum <- div[1:41,]
het_mum <- het[1:41,]
div_pup <- div[42:81,]
het_pup <- het[42:81,]

# creating one data frame
all <- data.frame(het$SH,div[,2:7])

# create age factor
age <- c(rep(1,41),rep(2,41))
age  <- factor(age)
levels(age) <- c("Mothers","Pups")
all$age <- age

#linear models for mums het vs. div
fm_S <- lm(div_mum$S ~ het_mum$SH)
fm_d <- lm(div_mum$d ~ het_mum$SH)
fm_Shannon <- lm(div_mum$H..loge. ~ het_mum$SH)
fm_Simpson <- lm(div_mum$X1.Lambda. ~ het_mum$SH)

#linear models for pups het vs. div
fm_S <- lm(div_pup$S ~ het_pup$SH)
fm_d <- lm(div_pup$d ~ het_pup$SH)
fm_Shannon <- lm(div_pup$H..loge. ~ het_pup$SH)
fm_Simpson <- lm(div_pup$X1.Lambda. ~ het_pup$SH)

summary(fm_S)
cor(div_mum$S,het_mum$SH)
plot(div_mum$X1.Lambda.~het_mum$SH)
abline(fm_Simpson)


#plotting everything

S <- qplot(het.SH,S, data=all) +
  geom_point(colour = "black", size = 2) +
  geom_smooth(method="lm",size = 1 ,col="black",alpha=0.15) +
  #geom_point(size = 3,colour="grey", alpha=0.2, show_guide = F) +          
  facet_wrap(~age, ncol=1) +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(vjust=1,size = 18, family="Arial"),
        strip.background = element_rect(fill="white",linetype="blank"),
        axis.title.x = element_text(vjust=0.1,size = 16, family="Arial"),
        axis.title.y = element_text(vjust=0.1,size = 16, family="Arial"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"),
        
       # axis.line = element_line(size=1),
       # panel.border = element_blank() ,
        panel.grid.major = element_blank(),
        axis.line = element_line(color = 'black'),
        panel.grid.minor = element_blank())+     
  geom_text(aes(0.8,80, label="r = 0.34**"),data=all) +
  xlab("Heterozygosity") +
  ylab("Scent profile diversity (number of compounds)") 

#ggsave("test.pdf",S)
#embed_fonts("test.pdf",outfile="test_embed.pdf")

        
d <- qplot(het.SH,d, data=all,color=factor(age)) +
        geom_smooth(method = lm,size=1.5,stat="smooth") +
        xlab("heterozygosity") +
        ylab("diversity (species richness)") +
        theme_set(theme_bw(20))

H <- qplot(het.SH,H..loge., data=all,color=factor(age)) +
        geom_smooth(method = lm,size=1.5,stat="smooth") +
        xlab("heterozygosity") +
        ylab("diversity (Shannon log)") +
        theme_set(theme_bw(20))

Simp <- qplot(het.SH,X1.Lambda., data=all,color=factor(age)) +
        geom_smooth(method = lm,size=1.5,stat="smooth") +
        xlab("heterozygosity") +
        ylab("diversity (Simpson)") +
        theme_set(theme_bw(20)) 

multiplot(S,d,H,Simp,cols=2)
