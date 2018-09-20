#########################################################################################################################
## @author: Abdou Rahmane WADE                                                                                         ##
##                                                                                                                     ##
##                           ######################################################################                    ##
##                           ##       Titre: Script pour detecter les gènes qui match avec       ##                    ##
##                           ##         les haplotypes détectés pot sous sélection               ##                    ##
##                           ######################################################################                    ##
##                                                                                                                     ## 
#########################################################################################################################
#script Gene Annotation
rm(list = ls())
#dir <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/données_ATGC/H12_all/peaks_all"
dir <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/données_ATGC/H12_center/peaks_center"
setwd(dir)

chr1 <- read.table("chr1_peaks")
chr2 <- read.table("chr2_peaks")
chr3 <- read.table("chr3_peaks")
chr4 <- read.table("chr4_peaks")
chr5 <- read.table("chr5_peaks")
chr6 <- read.table("chr6_peaks")
chr7 <- read.table("chr7_peaks")


names(chr1) <- c("l_pos_hap", "r_pos_hap", "nb_hap", "H12", "H2H1")
names(chr2) <- c("l_pos_hap", "r_pos_hap", "nb_hap", "H12", "H2H1")
names(chr3) <- c("l_pos_hap", "r_pos_hap", "nb_hap", "H12", "H2H1")
names(chr4) <- c("l_pos_hap", "r_pos_hap", "nb_hap", "H12", "H2H1")
names(chr5) <- c("l_pos_hap", "r_pos_hap", "nb_hap", "H12", "H2H1")
names(chr6) <- c("l_pos_hap", "r_pos_hap", "nb_hap", "H12", "H2H1")
names(chr7) <- c("l_pos_hap", "r_pos_hap", "nb_hap", "H12", "H2H1")

chr1 <- chr1[chr1$H12 >= 0.032,]
chr2 <- chr2[chr2$H12 >= 0.032,]
chr3 <- chr3[chr3$H12 >= 0.032,]
chr4 <- chr4[chr4$H12 >= 0.032,]
chr5 <- chr5[chr5$H12 >= 0.032,]
chr6 <- chr6[chr6$H12 >= 0.032,]
chr7 <- chr7[chr7$H12 >= 0.032,]

ref <- read.csv2("annotation.csv")

names(ref)

refchr1 <- ref[ref$chr=="chr1",]
write.table(refchr1, file="refchr1", sep="\t", quote=FALSE, col.names=F, row.names=F)
refchr2 <- ref[ref$chr=="chr2",]
write.table(refchr2, file="refchr2", sep="\t", quote=FALSE, col.names=F, row.names=F)
refchr3 <- ref[ref$chr=="chr3",]
write.table(refchr3, file="refchr3", sep="\t", quote=FALSE, col.names=F, row.names=F)
refchr4 <- ref[ref$chr=="chr4",]
write.table(refchr4, file="refchr4", sep="\t", quote=FALSE, col.names=F, row.names=F)
refchr5 <- ref[ref$chr=="chr5",]
write.table(refchr5, file="refchr5", sep="\t", quote=FALSE, col.names=F, row.names=F)
refchr6 <- ref[ref$chr=="chr6",]
write.table(refchr6, file="refchr6", sep="\t", quote=FALSE, col.names=F, row.names=F)
refchr7 <- ref[ref$chr=="chr7",]
write.table(refchr7, file="refchr7", sep="\t", quote=FALSE, col.names=F, row.names=F)

name <- c("chr", "gene_ref_start", "gene_ref_stop", "gene_ref_name", "l_pos_hap", "r_pos_hap", "nb_hap", "H12", "H2H1")
write.table(t(name), file="resrefchr1", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
write.table(t(name), file="resrefchr2", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
write.table(t(name), file="resrefchr3", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
write.table(t(name), file="resrefchr4", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
write.table(t(name), file="resrefchr5", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
write.table(t(name), file="resrefchr6", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
write.table(t(name), file="resrefchr7", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)

# chr1
refchr1  <- file("refchr1", open = "r")

while (length(oneLine <- readLines(refchr1, n = 1, warn = FALSE)) > 0) {
  
  oneLine<-unlist(strsplit(oneLine,split="\t"))
  
  for(i in 1:nrow(chr1)){
    if(((as.numeric(oneLine[4]) >= (chr1$l_pos_hap[i])) && (as.numeric(oneLine[4]) <= (chr1$r_pos_hap[i]))) || ((as.numeric(oneLine[5]) >= (chr1$l_pos_hap[i])) && (as.numeric(oneLine[5]) <= (chr1$r_pos_hap[i]))))
      write.table(t(c(oneLine[c(1,4,5,10)], chr1[i,])), file="resrefchr1", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  }
  
}


close(refchr1)

#chr2 

refchr2  <- file("refchr2", open = "r")

while (length(oneLine <- readLines(refchr2, n = 1, warn = FALSE)) > 0) {
  
  oneLine<-unlist(strsplit(oneLine,split="\t"))
  
  for(i in 1:nrow(chr2)){
    if(((as.numeric(oneLine[4]) >= (chr2$l_pos_hap[i])) && (as.numeric(oneLine[4]) <= (chr2$r_pos_hap[i]))) || ((as.numeric(oneLine[5]) >= (chr2$l_pos_hap[i])) && (as.numeric(oneLine[5]) <= (chr2$r_pos_hap[i]))))
      write.table(t(c(oneLine[c(1,4,5,10)], chr2[i,])), file="resrefchr2", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  }
  
}

close(refchr2)

#chr3 

refchr3  <- file("refchr3", open = "r")

while (length(oneLine <- readLines(refchr3, n = 1, warn = FALSE)) > 0) {
  
  oneLine<-unlist(strsplit(oneLine,split="\t"))
  
  for(i in 1:nrow(chr3)){
    if(((as.numeric(oneLine[4]) >= (chr3$l_pos_hap[i])) && (as.numeric(oneLine[4]) <= (chr3$r_pos_hap[i]))) || ((as.numeric(oneLine[5]) >= (chr3$l_pos_hap[i])) && (as.numeric(oneLine[5]) <= (chr3$r_pos_hap[i]))))
      write.table(t(c(oneLine[c(1,4,5,10)], chr3[i,])), file="resrefchr3", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  }
  
}

close(refchr3)


#chr4 

refchr4  <- file("refchr4", open = "r")

while (length(oneLine <- readLines(refchr4, n = 1, warn = FALSE)) > 0) {
  
  oneLine<-unlist(strsplit(oneLine,split="\t"))
  
  for(i in 1:nrow(chr4)){
    if(((as.numeric(oneLine[4]) >= (chr4$l_pos_hap[i])) && (as.numeric(oneLine[4]) <= (chr4$r_pos_hap[i]))) || ((as.numeric(oneLine[5]) >= (chr4$l_pos_hap[i])) && (as.numeric(oneLine[5]) <= (chr4$r_pos_hap[i]))))
      write.table(t(c(oneLine[c(1,4,5,10)], chr4[i,])), file="resrefchr4", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  }
  
}

close(refchr4)


#chr5 

refchr5  <- file("refchr5", open = "r")

while (length(oneLine <- readLines(refchr5, n = 1, warn = FALSE)) > 0) {
  
  oneLine<-unlist(strsplit(oneLine,split="\t"))
  
  for(i in 1:nrow(chr5)){
    if(((as.numeric(oneLine[4]) >= (chr5$l_pos_hap[i])) && (as.numeric(oneLine[4]) <= (chr5$r_pos_hap[i]))) || ((as.numeric(oneLine[5]) >= (chr5$l_pos_hap[i])) && (as.numeric(oneLine[5]) <= (chr5$r_pos_hap[i]))))
      write.table(t(c(oneLine[c(1,4,5,10)], chr5[i,])), file="resrefchr5", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  }
  
}

close(refchr5)


#chr6 

refchr6  <- file("refchr6", open = "r")

while (length(oneLine <- readLines(refchr6, n = 1, warn = FALSE)) > 0) {
  
  oneLine<-unlist(strsplit(oneLine,split="\t"))
  
  for(i in 1:nrow(chr6)){
    if(((as.numeric(oneLine[4]) >= (chr6$l_pos_hap[i])) && (as.numeric(oneLine[4]) <= (chr6$r_pos_hap[i]))) || ((as.numeric(oneLine[5]) >= (chr6$l_pos_hap[i])) && (as.numeric(oneLine[5]) <= (chr6$r_pos_hap[i]))))
      write.table(t(c(oneLine[c(1,4,5,10)], chr6[i,])), file="resrefchr6", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  }
  
}

close(refchr6)


#chr7 

refchr7  <- file("refchr7", open = "r")

while (length(oneLine <- readLines(refchr7, n = 1, warn = FALSE)) > 0) {
  
  oneLine<-unlist(strsplit(oneLine,split="\t"))
  
  for(i in 1:nrow(chr7)){
    if(((as.numeric(oneLine[4]) >= (chr7$l_pos_hap[i])) && (as.numeric(oneLine[4]) <= (chr7$r_pos_hap[i]))) || ((as.numeric(oneLine[5]) >= (chr7$l_pos_hap[i])) && (as.numeric(oneLine[5]) <= (chr7$r_pos_hap[i]))))
      write.table(t(c(oneLine[c(1,4,5,10)], chr7[i,])), file="resrefchr7", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  }
  
}

close(refchr7)


resrefchr1 <- read.table("resrefchr1", h=1)
resrefchr2 <- read.table("resrefchr2", h=1)
resrefchr3 <- read.table("resrefchr3", h=1)
resrefchr4 <- read.table("resrefchr4", h=1)
resrefchr5 <- read.table("resrefchr5", h=1)
resrefchr6 <- read.table("resrefchr6", h=1)
resrefchr7 <- read.table("resrefchr7", h=1)

resrefall <- rbind(resrefchr1,
                   resrefchr2,
                   resrefchr3,
                   resrefchr4,
                   resrefchr5,
                   resrefchr6,
                   resrefchr7)
names(resrefall) <- name
write.table(resrefall, file="resrefall", sep="\t", quote=FALSE, col.names=T, row.names=F)



resrefall <- read.table("resrefall",h=T)

tab <- data.frame(matrix(nrow = 8,ncol = 7))
tab[,1] <- rep("center",8)

tab[,2] <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","Allchr")

tab[,3] <- c(length(unique(resrefchr1$l_pos_hap)),
             length(unique(resrefchr2$l_pos_hap)),
             length(unique(resrefchr3$l_pos_hap)),
             length(unique(resrefchr4$l_pos_hap)),
             length(unique(resrefchr5$l_pos_hap)),
             length(unique(resrefchr6$l_pos_hap)),
             length(unique(resrefchr7$l_pos_hap)),
             length(unique(resrefall$l_pos_hap))
)

tab[,4] <- c(length(unique(chr1$l_pos_hap)),
             length(unique(chr2$l_pos_hap)),
             length(unique(chr3$l_pos_hap)),
             length(unique(chr4$l_pos_hap)),
             length(unique(chr5$l_pos_hap)),
             length(unique(chr6$l_pos_hap)),
             length(unique(chr7$l_pos_hap)),0)
tab[8,4] <- sum(tab[1:7,4])

chr1_sf <- read.table("chr1_peaks")
chr2_sf <- read.table("chr2_peaks")
chr3_sf <- read.table("chr3_peaks")
chr4_sf <- read.table("chr4_peaks")
chr5_sf <- read.table("chr5_peaks")
chr6_sf <- read.table("chr6_peaks")
chr7_sf <- read.table("chr7_peaks")

tab[,5] <- c(length(unique(chr1_sf$V1)),
             length(unique(chr2_sf$V1)),
             length(unique(chr3_sf$V1)),
             length(unique(chr4_sf$V1)),
             length(unique(chr5_sf$V1)),
             length(unique(chr6_sf$V1)),
             length(unique(chr7_sf$V1)),0)
tab[8,5] <- sum(tab[1:7,5])

tab[,6] <- c(length(unique(resrefchr1$gene_ref_name)),
             length(unique(resrefchr2$gene_ref_name)),
             length(unique(resrefchr3$gene_ref_name)),
             length(unique(resrefchr4$gene_ref_name)),
             length(unique(resrefchr5$gene_ref_name)),
             length(unique(resrefchr6$gene_ref_name)),
             length(unique(resrefchr7$gene_ref_name)),
             length(unique(resrefall$gene_ref_name))
)

refchr1 <- ref[ref$chr=="chr1",]
refchr2 <- ref[ref$chr=="chr2",]
refchr3 <- ref[ref$chr=="chr3",]
refchr4 <- ref[ref$chr=="chr4",]
refchr5 <- ref[ref$chr=="chr5",]
refchr6 <- ref[ref$chr=="chr6",]
refchr7 <- ref[ref$chr=="chr7",]

tab[,7] <- c(length(unique(refchr1$name)),
             length(unique(refchr2$name)),
             length(unique(refchr3$name)),
             length(unique(refchr4$name)),
             length(unique(refchr5$name)),
             length(unique(refchr6$name)),
             length(unique(refchr7$name)),
             length(unique(ref$name))
)
names(tab) <- c("ech","chr","nb Wind mapped", "nb Wind sup 0.014", "nb Wind tot", "nb genes mapped", "nb gene ref")
write.table(tab, file="tabresult", sep="\t", quote=FALSE, col.names=T, row.names=F)


tab <- read.table("tabresult", sep="\t", h=T)
colonn <- c("mapped", "nomapped")
lign <- c("wind", "gene")
mat_chr1 <- matrix(nrow = 2,ncol = 2)
mat_chr2 <- matrix(nrow = 2,ncol = 2)
mat_chr3 <- matrix(nrow = 2,ncol = 2)
mat_chr4 <- matrix(nrow = 2,ncol = 2)
mat_chr5 <- matrix(nrow = 2,ncol = 2)
mat_chr6 <- matrix(nrow = 2,ncol = 2)
mat_chr7 <- matrix(nrow = 2,ncol = 2)
mat_tot <- matrix(nrow = 2,ncol = 2)
colnames(mat_chr1) <- colonn
colnames(mat_chr2) <- colonn
colnames(mat_chr3) <- colonn
colnames(mat_chr4) <- colonn
colnames(mat_chr5) <- colonn
colnames(mat_chr6) <- colonn
colnames(mat_chr7) <- colonn
colnames(mat_tot) <- colonn
row.names(mat_chr1) <- lign
row.names(mat_chr2) <- lign
row.names(mat_chr3) <- lign
row.names(mat_chr4) <- lign
row.names(mat_chr5) <- lign
row.names(mat_chr6) <- lign
row.names(mat_chr7) <- lign
row.names(mat_tot) <- lign

#chr1
mat_chr1[1,1]<- tab$nb.Wind.sup.0.014[1]
mat_chr1[1,2]<- tab$nb.Wind.tot[1] - tab$nb.Wind.sup.0.014[1]
mat_chr1[2,1]<- tab$nb.genes.mapped[1]
mat_chr1[2,2]<- tab$nb.gene.ref[1] - tab$nb.genes.mapped[1]
mat_chr1

library(DescTools)

GTest(mat_chr1, correct = "none")

#chr2
mat_chr2[1,1]<- tab$nb.Wind.sup.0.014[2]
mat_chr2[1,2]<- tab$nb.Wind.tot[2] - tab$nb.Wind.sup.0.014[2]
mat_chr2[2,1]<- tab$nb.genes.mapped[2]
mat_chr2[2,2]<- tab$nb.gene.ref[2] - tab$nb.genes.mapped[2]

library(DescTools)

GTest(mat_chr2, correct = "none")

#chr3
mat_chr3[1,1]<- tab$nb.Wind.sup.0.014[3]
mat_chr3[1,2]<- tab$nb.Wind.tot[3] - tab$nb.Wind.sup.0.014[3]
mat_chr3[2,1]<- tab$nb.genes.mapped[3]
mat_chr3[2,2]<- tab$nb.gene.ref[3] - tab$nb.genes.mapped[3]
mat_chr3

library(DescTools)

GTest(mat_chr3, correct = "none")

#chr4
mat_chr4[1,1]<- tab$nb.Wind.sup.0.014[4]
mat_chr4[1,2]<- tab$nb.Wind.tot[4] - tab$nb.Wind.sup.0.014[4]
mat_chr4[2,1]<- tab$nb.genes.mapped[4]
mat_chr4[2,2]<- tab$nb.gene.ref[4] - tab$nb.genes.mapped[4]

library(DescTools)

GTest(mat_chr4, correct = "none")

#chr5
mat_chr5[1,1]<- tab$nb.Wind.sup.0.014[5]
mat_chr5[1,2]<- tab$nb.Wind.tot[5] - tab$nb.Wind.sup.0.014[5]
mat_chr5[2,1]<- tab$nb.genes.mapped[5]
mat_chr5[2,2]<- tab$nb.gene.ref[5] - tab$nb.genes.mapped[5]
mat_chr5

library(DescTools)

GTest(mat_chr5, correct = "none")

#chr6
mat_chr6[1,1]<- tab$nb.Wind.sup.0.014[6]
mat_chr6[1,2]<- tab$nb.Wind.tot[6] - tab$nb.Wind.sup.0.014[6]
mat_chr6[2,1]<- tab$nb.genes.mapped[6]
mat_chr6[2,2]<- tab$nb.gene.ref[6] - tab$nb.genes.mapped[6]
mat_chr6

library(DescTools)

GTest(mat_chr6, correct = "none")

#chr7
mat_chr7[1,1]<- tab$nb.Wind.sup.0.014[7]
mat_chr7[1,2]<- tab$nb.Wind.tot[7] - tab$nb.Wind.sup.0.014[7]
mat_chr7[2,1]<- tab$nb.genes.mapped[7]
mat_chr7[2,2]<- tab$nb.gene.ref[7] - tab$nb.genes.mapped[7]
mat_chr7

library(DescTools)

GTest(mat_chr7, correct = "none")

#tot
mat_tot[1,1]<- tab$nb.Wind.sup.0.014[8]
mat_tot[1,2]<- tab$nb.Wind.tot[8] - tab$nb.Wind.sup.0.014[8]
mat_tot[2,1]<- tab$nb.genes.mapped[8]
mat_tot[2,2]<- tab$nb.gene.ref[8] - tab$nb.genes.mapped[8]
mat_tot

library(DescTools)

GTest(mat_tot, correct = "none")




#####################################################################################""""
pdf("fig1.pdf",height = 8, width = 60)
par(mfrow=c(1,1))
plot(1:280000, rep(2,280000 ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "" , xaxt="n")
segments((resrefchr1$gene_ref_start / 1000), 1, x1=(resrefchr1$gene_ref_stop / 1000), y1=1, size =10, xlim = c(0,280000),col='red')
segments((resrefchr1$l_pos_hap / 1000), 1.002, x1=(resrefchr1$r_pos_hap / 1000), y1=1.002, size =10, xlim = c(0,280000),col='green')
title(main = "chr1")

plot(1:3000, rep(2,3000 ), ylim = c(0.9,1.1), yaxt="n" , xlab = "Position (kb)", ylab = "" , xaxt="n")
segments((resrefchr1$gene_ref_start / 1000), 1, x1=(resrefchr1$gene_ref_stop / 1000), y1=1, size =10, xlim = c(1,3000),col='red')
segments((resrefchr1$l_pos_hap / 1000), 1.002, x1=(resrefchr1$r_pos_hap / 1000), y1=1.002, size =10, xlim = c(1,3000),col='green')
title(main = "chr1")