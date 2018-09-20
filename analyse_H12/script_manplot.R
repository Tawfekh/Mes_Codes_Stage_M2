#############################################################################################################
## Author:Abdou Rahmane WADE                                                             Date:      2018   ##
##                                                                                                         ##
##          ##################################################################################             ##
##          ##    Titre: Script pour générer un scan genomique de statistique H12           ##             ## 
##          ##                           Description:                                       ##             ##
##          ##################################################################################             ##
##                                                                                                         ##
#############################################################################################################

setwd("/home/utlf0032/Bureau/zap/centerH12/")
rm(list = ls())
Hchr1 <- read.table("center_chr1_h12_output_400_50.txt")
Hchr2 <- read.table("center_chr2_h12_output_400_50.txt")
Hchr3 <- read.table("center_chr3_h12_output_400_50.txt")
Hchr4 <- read.table("center_chr4_h12_output_400_50.txt")
Hchr5 <- read.table("center_chr5_h12_output_400_50.txt")
Hchr6 <- read.table("center_chr6_h12_output_400_50.txt")
Hchr7 <- read.table("center_chr7_h12_output_400_50.txt")
CHR <- rep(1, nrow(Hchr1))
chr1 <- data.frame(cbind(CHR,Hchr1[,1], Hchr1[,9]))
names(chr1)<-c("CHR","BP","P")
                     
CHR <- rep(2, nrow(Hchr2))
chr2 <- data.frame(cbind(CHR,Hchr2[,1], Hchr2[,9]))
names(chr2)<-c("CHR","BP","P")

CHR <- rep(3, nrow(Hchr3))
chr3 <- data.frame(cbind(CHR,Hchr3[,1], Hchr3[,9]))
names(chr3)<-c("CHR","BP","P")

CHR <- rep(4, nrow(Hchr4))
chr4 <- data.frame(cbind(CHR,Hchr4[,1], Hchr4[,9]))
names(chr4)<-c("CHR","BP","P")

CHR <- rep(5, nrow(Hchr5))
chr5 <- data.frame(cbind(CHR,Hchr5[,1], Hchr5[,9]))
names(chr5)<-c("CHR","BP","P")

CHR <- rep(6, nrow(Hchr6))
chr6 <- data.frame(cbind(CHR,Hchr6[,1], Hchr6[,9]))
names(chr6)<-c("CHR","BP","P")

CHR <- rep(7, nrow(Hchr7))
chr7 <- data.frame(cbind(CHR,Hchr7[,1], Hchr7[,9]))
names(chr7)<-c("CHR","BP","P")

chr <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7)

don1 <-chr[chr$P>=0.03556,]

library(dplyr)


#don <- don1 %>% 
don <- chr %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
#  left_join(don1, ., by=c("CHR"="CHR")) %>%
  
  left_join(chr, ., by=c("CHR"="CHR")) %>%
  
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)


axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

library(ggplot2)

#pdf(file = "mahattan_H12_ALL_seuil_1pargenome_03556.pdf", width = 14, height = 7)
pdf(file = "mahattan_H12_ALL_seuil_5percent_0211.pdf", width = 14, height = 7)
ggplot(don, aes(x=BPcum, y=P)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.5) +
  geom_hline(yintercept = 0.03556, color="red")+
  scale_color_manual(values = rep(c("grey", "skyblue"), 7 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(breaks = seq(0,0.3,0.05) ) +     # remove space between plot area and x axis
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  labs(x="Chromosomes", y="H12")
dev.off()

quantile(chr$P,.95)

