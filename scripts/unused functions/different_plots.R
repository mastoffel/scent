S1 <- qplot(het.SH,S, data=all) +
  geom_smooth(method="lm",size = 1 ,col="black",alpha=0.1,fill=element_blank())+
  geom_point(aes(size = 1.5), show_guide = F) +          
  facet_wrap(~age) +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 18, family="Arial"),
        strip.background = element_rect(fill="white",linetype="blank"),
        axis.title.x = element_text(size = 16, family="Arial"),
        axis.title.y = element_text(size = 16, family="Arial"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        # axis.line = element_line(size=1),
        panel.border = element_blank() ,
        panel.grid.major = element_blank(),
        axis.line = element_line(color = 'black'))+
  #panel.grid.minor = element_blank())+            
  xlab("Heterozygosity") +
  ylab("Scent profile diversity") 

S2 <- qplot(het.SH,S, data=all) +
  geom_smooth(method="lm",size = 1 ,col="black",alpha=0.1,fill=element_blank())+
  geom_point(aes(size = 1.5), show_guide = F) +          
  facet_wrap(~age) +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 18, family="Arial"),
        strip.background = element_rect(fill="white",linetype="blank"),
        axis.title.x = element_text(size = 16, family="Arial"),
        axis.title.y = element_text(size = 16, family="Arial"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
        # axis.line = element_line(size=1),
        
  #panel.grid.minor = element_blank())+            
  xlab("Heterozygosity") +
  ylab("Scent profile diversity") 


S3 <- qplot(het.SH,S, data=all) +
  geom_smooth(method="lm",size = 1 ,col="black",alpha=0.1,fill=element_blank())+
  geom_point(aes(size = 1.5), show_guide = F) +          
  facet_wrap(~age) +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 18, family="Arial"),
        strip.background = element_rect(fill="white",linetype="blank"),
        axis.title.x = element_text(size = 16, family="Arial"),
        axis.title.y = element_text(size = 16, family="Arial"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  # axis.line = element_line(size=1),
  #panel.grid.minor = element_blank())+            
  xlab("Heterozygosity") +
  ylab("Scent profile diversity") 

S4 <- qplot(het.SH,S, data=all) +
  geom_smooth(method="lm",size = 1 ,col="black",alpha=0.1,fill=element_blank())+
  geom_point(aes(size = 1.5), show_guide = F) +          
  facet_wrap(~age) +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_text(size = 18, family="Arial"),
        strip.background = element_rect(fill="white",linetype="blank"),
        axis.title.x = element_text(size = 16, family="Arial"),
        axis.title.y = element_text(size = 16, family="Arial"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black",size=1.5))+
  # axis.line = element_line(size=1),
  #panel.grid.minor = element_blank())+            
  xlab("Heterozygosity") +
  ylab("Scent profile diversity") 