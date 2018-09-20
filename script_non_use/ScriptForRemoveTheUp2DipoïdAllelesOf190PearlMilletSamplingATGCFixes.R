#' @author Abdou Rahmane WADE
#' @title Script for remove the 'up to dipoîd' alleles of 190 pearl millet sampling for 7 chr
#' @input .geno_File
#' @output .geno_File

####Chr1

rm(list=ls())

setwd("/home/utlf0032/Documents/mon_doc_zap/mil_data")

inputFile <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr1_sample_190_of_allSNP_chr1.txt"

outputFile<- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr1_sample_190_of_DIPSNP_chr1.txt"

enleve <- read.table("chr1_tri-four_alleles.txt", sep = " ", h =1)

position <- as.character(enleve$position)

chr1  <- file(inputFile, open = "r")

oneLine <- readLines(chr1, n = 1, warn = FALSE)

Sample.ID <- unlist(strsplit(oneLine, split = " "))

writeLines(Sample.ID, con = outputFile, sep = " ")

write.table("", file = outputFile, sep = "/n", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

while (length(oneLine <- unlist(strsplit(readLines(chr1, n = 1, warn = FALSE),split = " "))) > 0) 
  {
  
  
  if (min(oneLine[2]!=position)>0)
    {
      text<-t(c(oneLine[1],as.numeric(oneLine[2]),oneLine[3:length(oneLine)]))
      write.table(text, file = outputFile, sep = " ", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
    }
  }

close(chr1) 


####Chr2

rm(list=ls)

setwd("/home/utlf0032/Documents/mon_doc_zap/mil_data")

inputFile <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr2_sample_190_of_allSNP_chr2.txt"

outputFile<- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr2_sample_190_of_DIPSNP_chr2.txt"

enleve <- read.table("chr2_tri-four_alleles.txt", sep = " ", h =1)

position <- as.character(enleve$position)

chr2  <- file(inputFile, open = "r")

oneLine <- readLines(chr2, n = 1, warn = FALSE)

Sample.ID <- unlist(strsplit(oneLine, split = " "))

writeLines(Sample.ID, con = outputFile, sep = " ")

write.table("", file = outputFile, sep = "/n", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

while (length(oneLine <- unlist(strsplit(readLines(chr2, n = 1, warn = FALSE),split = " "))) > 0) 
{
  
  
  if (min(oneLine[2]!=position)>0)
  {
    text<-t(c(oneLine[1],as.numeric(oneLine[2]),oneLine[3:length(oneLine)]))
    write.table(text, file = outputFile, sep = " ", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
  }
}

close(chr2) 


####chr3

rm(list=ls)

setwd("/home/utlf0032/Documents/mon_doc_zap/mil_data")

inputFile <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr3_sample_190_of_allSNP_chr3.txt"

outputFile<- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr3_sample_190_of_DIPSNP_chr3.txt"

enleve <- read.table("chr3_tri-four_alleles.txt", sep = " ", h =1)

position <- as.character(enleve$position)

chr3  <- file(inputFile, open = "r")

oneLine <- readLines(chr3, n = 1, warn = FALSE)

Sample.ID <- unlist(strsplit(oneLine, split = " "))

writeLines(Sample.ID, con = outputFile, sep = " ")

write.table("", file = outputFile, sep = "/n", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

while (length(oneLine <- unlist(strsplit(readLines(chr3, n = 1, warn = FALSE),split = " "))) > 0) 
{
  
  
  if (min(oneLine[2]!=position)>0)
  {
    text<-t(c(oneLine[1],as.numeric(oneLine[2]),oneLine[3:length(oneLine)]))
    write.table(text, file = outputFile, sep = " ", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
  }
}

close(chr3) 



####chr4

rm(list=ls)

setwd("/home/utlf0032/Documents/mon_doc_zap/mil_data")

inputFile <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr4_sample_190_of_allSNP_chr4.txt"

outputFile<- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr4_sample_190_of_DIPSNP_chr4.txt"

enleve <- read.table("chr4_tri-four_alleles.txt", sep = " ", h =1)

position <- as.character(enleve$position)

chr4  <- file(inputFile, open = "r")

oneLine <- readLines(chr4, n = 1, warn = FALSE)

Sample.ID <- unlist(strsplit(oneLine, split = " "))

writeLines(Sample.ID, con = outputFile, sep = " ")

write.table("", file = outputFile, sep = "/n", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

while (length(oneLine <- unlist(strsplit(readLines(chr4, n = 1, warn = FALSE),split = " "))) > 0) 
{
  
  
  if (min(oneLine[2]!=position)>0)
  {
    text<-t(c(oneLine[1],as.numeric(oneLine[2]),oneLine[3:length(oneLine)]))
    write.table(text, file = outputFile, sep = " ", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
  }
}

close(chr4) 


####chr5

rm(list=ls)

setwd("/home/utlf0032/Documents/mon_doc_zap/mil_data")

inputFile <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr5_sample_190_of_allSNP_chr5.txt"

outputFile<- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr5_sample_190_of_DIPSNP_chr5.txt"

enleve <- read.table("chr5_tri-four_alleles.txt", sep = " ", h =1)

position <- as.character(enleve$position)

chr5  <- file(inputFile, open = "r")

oneLine <- readLines(chr5, n = 1, warn = FALSE)

Sample.ID <- unlist(strsplit(oneLine, split = " "))

writeLines(Sample.ID, con = outputFile, sep = " ")

write.table("", file = outputFile, sep = "/n", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

while (length(oneLine <- unlist(strsplit(readLines(chr5, n = 1, warn = FALSE),split = " "))) > 0) 
{
  
  
  if (min(oneLine[2]!=position)>0)
  {
    text<-t(c(oneLine[1],as.numeric(oneLine[2]),oneLine[3:length(oneLine)]))
    write.table(text, file = outputFile, sep = " ", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
  }
}

close(chr5) 



####chr6

rm(list=ls)

setwd("/home/utlf0032/Documents/mon_doc_zap/mil_data")

inputFile <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr6_sample_190_of_allSNP_chr6.txt"

outputFile<- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr6_sample_190_of_DIPSNP_chr6.txt"

enleve <- read.table("chr6_tri-four_alleles.txt", sep = " ", h =1)

position <- as.character(enleve$position)

chr6  <- file(inputFile, open = "r")

oneLine <- readLines(chr6, n = 1, warn = FALSE)

Sample.ID <- unlist(strsplit(oneLine, split = " "))

writeLines(Sample.ID, con = outputFile, sep = " ")

write.table("", file = outputFile, sep = "/n", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

while (length(oneLine <- unlist(strsplit(readLines(chr6, n = 1, warn = FALSE),split = " "))) > 0) 
{
  
  
  if (min(oneLine[2]!=position)>0)
  {
    text<-t(c(oneLine[1],as.numeric(oneLine[2]),oneLine[3:length(oneLine)]))
    write.table(text, file = outputFile, sep = " ", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
  }
}

close(chr6) 


####chr7

rm(list=ls)

setwd("/home/utlf0032/Documents/mon_doc_zap/mil_data")

inputFile <- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr7_sample_190_of_allSNP_chr7.txt"

outputFile<- "/home/utlf0032/Documents/mon_doc_zap/mil_data/chr7_sample_190_of_DIPSNP_chr7.txt"

enleve <- read.table("chr7_tri-four_alleles.txt", sep = " ", h =1)

position <- as.character(enleve$position)

chr7  <- file(inputFile, open = "r")

oneLine <- readLines(chr7, n = 1, warn = FALSE)

Sample.ID <- unlist(strsplit(oneLine, split = " "))

writeLines(Sample.ID, con = outputFile, sep = " ")

write.table("", file = outputFile, sep = "/n", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

while (length(oneLine <- unlist(strsplit(readLines(chr7, n = 1, warn = FALSE),split = " "))) > 0) 
{
  
  
  if (min(oneLine[2]!=position)>0)
  {
    text<-t(c(oneLine[1],as.numeric(oneLine[2]),oneLine[3:length(oneLine)]))
    write.table(text, file = outputFile, sep = " ", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
  }
}

close(chr7) 
