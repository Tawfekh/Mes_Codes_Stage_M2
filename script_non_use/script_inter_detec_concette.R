#########################################################################################################################
## @author: Abdou Rahmane WADE                                                                                         ##
##                                                                                                                     ##
##                           ######################################################################                    ##
##                           ##       Titre: Script pour detecter les gènes communs entre        ##                    ##
##                           ##         ma detection et celle de Burgarella et al 2018           ##                    ##
##                           ######################################################################                    ##
##                                                                                                                     ## 
#########################################################################################################################

#script Gene Annotation
rm(list = ls())
#dir <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/données_ATGC/H12_all/peaks_all"
dir <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/données_ATGC/H12_center/peaks_center"
setwd(dir)
resrefall <- read.table("resrefall",h=T)

hard_sel_tot_tab <- read.csv2("Regions_under_selection_concetta.csv")

hard_sel <- hard_sel_tot_tab[,c(1:6)]

hard_sel <- hard_sel[hard_sel$Cultivated_group=="cult_C",]

names(hard_sel)[2] <- names(resrefall)[4]

intersect <- merge(resrefall, hard_sel, by=names(resrefall)[4])

intersect_chr1 <- intersect[intersect$chr == "chr1",]

library("RColorBrewer")

plot(21501:21700, rep(2,200 ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "" , xaxt="n")
segments((intersect_chr1$gene_ref_start / 1000), 1, x1=(intersect_chr1$gene_ref_stop / 1000), y1=1, size =10, xlim = c(21501,21700),col='red')
l=0
for(i in 1:length(unique(intersect_chr1$l_pos_hap))){
  segments((unique(intersect_chr1$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr1$r_pos_hap)[i] / 1000), y1=1.002 +l, size =10, xlim = c(21501,21700),col= "green")
  l= l + 0.001
}
axis(1, at = seq(21500,21700, by = 10), las=2)
title(main = "chr1")

##

intersect_chr2 <- intersect[intersect$chr == "chr2",]

library("RColorBrewer")

plot(119001:119050, rep(2,50 ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "" )
segments((intersect_chr2$gene_ref_start / 1000), 1, x1=(intersect_chr2$gene_ref_stop / 1000), y1=1, size =10, xlim = c(119001,119050),col='red')
l=0
for(i in 1:length(unique(intersect_chr2$l_pos_hap))){
  segments((unique(intersect_chr2$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr2$r_pos_hap)[i] / 1000), y1=1.002 +l, 
           size =10, xlim = c(119001,119050),col= "green")
  l= l + 0.001
}
axis(1, at = seq(c(119000,119050), by = 1), las=2)
title(main = "chr2")

##

intersect_chr3 <- intersect[intersect$chr == "chr3",]

library("RColorBrewer")

plot(96190:96750, rep(2,561 ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "" )
segments((intersect_chr3$gene_ref_start / 1000), 1, x1=(intersect_chr3$gene_ref_stop / 1000), y1=1, size =10, xlim = c(96190,96750),col='red')
l=0
for(i in 1:length(unique(intersect_chr3$l_pos_hap))){
  segments((unique(intersect_chr3$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr3$r_pos_hap)[i] / 1000), y1=1.002 +l, 
           size =10, xlim = c(96190,96750),col= "green")
  l= l + 0.001
}
title(main = "chr3")


##

intersect_chr4 <- intersect[intersect$chr == "chr4",]

library("RColorBrewer")
 m <- min(intersect_chr4[,c(3:6)] )
 M <- max(intersect_chr4[,c(3:6)] )
 m <- ((m/1000) - 20) 
 M <- ((M/1000) - 20)
 L <- length(c(m:M))

plot(m:M, rep(2,L ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "" )
segments((intersect_chr4$gene_ref_start / 1000), 1, x1=(intersect_chr4$gene_ref_stop / 1000), y1=1, size =10, xlim = c(m,M),col='red')
l=0
for(i in 1:length(unique(intersect_chr4$l_pos_hap))){
  segments((unique(intersect_chr4$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr4$r_pos_hap)[i] / 1000), y1=1.002 +l, 
           size =10, xlim = c(m,M),col= "green")
  l= l + 0.001
}
axis(1, at = seq(c(m,M), by = 1), las=2)
title(main = "chr4")

plot(m:35790, rep(2,63 ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "" )
segments((intersect_chr4$gene_ref_start / 1000), 1, x1=(intersect_chr4$gene_ref_stop / 1000), y1=1, size =10, xlim = c(m,M),col='red')
l=0
for(i in 1:length(unique(intersect_chr4$l_pos_hap))){
  segments((unique(intersect_chr4$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr4$r_pos_hap)[i] / 1000), y1=1.002 +l, 
           size =10, xlim = c(m,M),col= "green")
  l= l + 0.001
}
axis(1, at = seq(c(m,M), by = 1), las=2)
title(main = "chr4")

plot(89600:(M+20), rep(2,225 ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "" )
segments((intersect_chr4$gene_ref_start / 1000), 1, x1=(intersect_chr4$gene_ref_stop / 1000), y1=1, size =10, xlim = c(m,M),col='red')
l=0
for(i in 1:length(unique(intersect_chr4$l_pos_hap))){
  segments((unique(intersect_chr4$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr4$r_pos_hap)[i] / 1000), y1=1.002 +l, 
           size =10, xlim = c(m,M),col= "green")
  l= l + 0.001
}
axis(1, at = seq(c(m,M), by = 1), las=2)
title(main = "chr4")



##
intersect_chr5 <- intersect[intersect$chr == "chr5",]

library("RColorBrewer")
m <- min(intersect_chr5[,c(3:6)] )
M <- max(intersect_chr5[,c(3:6)] )
m <- ((m/1000) - 20) 
M <- ((M/1000) + 20)
L <- length(c(m:M))

plot(m:M, rep(2,L ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "" )
segments((intersect_chr5$gene_ref_start / 1000), 1, x1=(intersect_chr5$gene_ref_stop / 1000), y1=1, size =10, xlim = c(m,M),col='red')
l=0
for(i in 1:length(unique(intersect_chr5$l_pos_hap))){
  segments((unique(intersect_chr5$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr5$r_pos_hap)[i] / 1000), y1=1.002 +l, 
           size =10, xlim = c(m,M),col= "green")
  l= l + 0.001
}
axis(1, at = seq(c(m,M), by = 1), las=2)
title(main = "chr5")



##
intersect_chr6 <- intersect[intersect$chr == "chr6",]

library("RColorBrewer")
m <- min(intersect_chr6[,c(3:6)] )
M <- max(intersect_chr6[,c(3:6)] )
m <- ((m/1000) - 20) 
M <- ((M/1000) + 20)
L <- length(c(m:M))

plot(m:M, rep(2,L ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "",xaxt="n" )
segments((intersect_chr6$gene_ref_start / 1000), 1, x1=(intersect_chr6$gene_ref_stop / 1000), y1=1, size =10, xlim = c(m,M),col='red')
l=0
for(i in 1:length(unique(intersect_chr6$l_pos_hap))){
  segments((unique(intersect_chr6$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr6$r_pos_hap)[i] / 1000), y1=1.002 +l, 
           size =10, xlim = c(m,M),col= "green")
  l= l + 0.001
}
axis(1, at = seq(m,M, by = 1000), las=2)
title(main = "chr6")


  plot(trunc(m):86770, rep(2,length(trunc(m):86770)), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "",xaxt="n" )
  segments((intersect_chr6$gene_ref_start / 1000), 1, x1=(intersect_chr6$gene_ref_stop / 1000), y1=1, size =10, xlim = c(m,M),col='red')
  l=0
  for(i in 1:length(unique(intersect_chr6$l_pos_hap))){
    segments((unique(intersect_chr6$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr6$r_pos_hap)[i] / 1000), y1=1.002 +l, 
             size =10, xlim = c(m,M),col= "green")
    l= l + 0.001
  }
  axis(1, at = seq(trunc(m),86770, by = 10), las=2)
  title(main = "chr6")
  
  
  plot(221970:222022, rep(2,length(221970:222022)), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "",xaxt="n" )
  segments((intersect_chr6$gene_ref_start / 1000), 1, x1=(intersect_chr6$gene_ref_stop / 1000), y1=1, size =10, xlim = c(m,M),col='red')
  l=0
  for(i in 1:length(unique(intersect_chr6$l_pos_hap))){
    segments((unique(intersect_chr6$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr6$r_pos_hap)[i] / 1000), y1=1.002 +l, 
             size =10, xlim = c(m,M),col= "green")
    l= l + 0.001
  }
  axis(1, at = seq(221970,222022, by = 10), las=2)
  title(main = "chr6")
  
  plot(227350:trunc(M), rep(2,length(227350:trunc(M))), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "",xaxt="n" )
  segments((intersect_chr6$gene_ref_start / 1000), 1, x1=(intersect_chr6$gene_ref_stop / 1000), y1=1, size =10, xlim = c(m,M),col='red')
  l=0
  for(i in 1:length(unique(intersect_chr6$l_pos_hap))){
    segments((unique(intersect_chr6$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr6$r_pos_hap)[i] / 1000), y1=1.002 +l, 
             size =10, xlim = c(m,M),col= "green")
    l= l + 0.001
  }
  axis(1, at = seq(227350,trunc(M), by = 100), las=2)
  title(main = "chr6")
  
  
  ##
  intersect_chr7 <- intersect[intersect$chr == "chr7",]
  
  library("RColorBrewer")
  m <- min(intersect_chr7[,c(3:6)] )
  M <- max(intersect_chr7[,c(3:6)] )
  m <- ((m/1000) - 20) 
  M <- ((M/1000) + 20)
  L <- length(c(m:M))
  
  plot(m:M, rep(2,L ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "",xaxt="n" )
  segments((intersect_chr7$gene_ref_start / 1000), 1, x1=(intersect_chr7$gene_ref_stop / 1000), y1=1, size =10, xlim = c(m,M),col='red')
  l=0
  for(i in 1:length(unique(intersect_chr7$l_pos_hap))){
    segments((unique(intersect_chr7$l_pos_hap)[i] / 1000), 1.002 + l, x1=(unique(intersect_chr7$r_pos_hap)[i] / 1000), y1=1.002 +l, 
             size =10, xlim = c(m,M),col= "green")
    l= l + 0.001
  }
  axis(1, at = seq(m,M, by = 10), las=2)
  title(main = "chr7")
  