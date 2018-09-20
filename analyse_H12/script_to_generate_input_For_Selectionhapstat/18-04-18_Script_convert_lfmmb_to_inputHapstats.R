
#########################################################################################################################
## @author: Abdou Rahmane WADE                                                                    Date:18-04-2018      ##
##                                                                                                                     ##
##                           ######################################################################                    ##
##                           ## Titre: Script for generate SelectionHapstats (Garud et al 2015)  ##                    ##
##                           ##                       Description:                               ##                    ##
##                           ######################################################################                    ##
##                                                                                                                     ##
##             - Haplodiploidise ei Fixe un alléles au hasard pour les locus diploïdes                                 ##
##             - Code les données manquantes par N                                                                     ##
##             - Separateur = ,                                                                                        ##
##             - Exemple format fichier d'entrée:                                                                      ##
##                                                                                                                     ##
##                CHR1 Position ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8 ind9 ind10 ind11 ind12 ind13                   ##
##                chr1 35914 G G G G G G G G G G G G G                                                                 ##
##                chr1 35959 T K G T G G T G T G K T G                                                                 ##
##                chr1 36170 T - T T T T T T T T T T -                                                                 ##
##                chr1 36202 A - A A A A A A A - A A A                                                                 ##
##                chr1 36206 A - A A A A A A A - A A A                                                                 ##
##                chr1 36207 C - C C C C C C C - C C C                                                                 ##
##                chr1 36226 M - C A C C A C A - C C C                                                                 ##
##                chr1 36259 T - T T - T T T R - T T T                                                                 ##
##                chr1 36284 C - - - C C - C C - C C C                                                                 ##
##                                                                                                                     ##
##             - Exemple format fichier de sortie:                                                                     ##
##                                                                                                                     ##
##                35914,G,G,G,G,G,G,G,G,G,G,G,G,G                                                                      ##
##                35959,T,T,G,T,G,G,T,G,T,G,G,T,G                                                                      ##
##                36170,T,N,T,T,T,T,T,T,T,T,T,T,N                                                                      ##
##                36202,A,N,A,A,A,A,A,A,A,N,A,A,A                                                                      ##
##                36206,A,N,A,A,A,A,A,A,A,N,A,A,A                                                                      ##
##                36207,C,N,C,C,C,C,C,C,C,N,C,C,C                                                                      ##
##                36226,A,N,C,A,C,C,A,C,A,N,C,C,C                                                                      ##
##                36259,T,N,T,T,N,T,T,T,G,N,T,T,T                                                                      ##
##                36284,C,N,N,N,C,C,N,C,C,N,C,C,C                                                                      ##
##                                                                                                                     ##
##                                                                                                                     ##
##                                                                                                                     ##
##            InputFile: each files in 190_pearl_millet_Haploid_WithSomeDiploidsLoci_data_lfmmb.tar.gz                 ##
##            OutputFile: each files in 20-04-18_190_pearl_millet_Haploid_data_HapstatsInputFiles.tar.gz.gz            ##
##                                                                                                                     ##
##            @ Script modifié à partir de convert_ch_to_lfmmb.R (auteur: ??)                                          ##
##                                                                                                                     ## ##                                                                                                                     ##
#########################################################################################################################



args <- commandArgs(TRUE)

# name inputFile
inputFile <- args[1]

# name outputFile
L <- length(unlist(strsplit(inputFile,split = "")))
outputFile<- paste(substring(inputFile, 1 , L-4),"_h12.txt",sep="")

# open connection
f1  <- file(inputFile, open = "r")

# read first line
readLines(f1, n = 1, warn = FALSE)

i<-1
while (length(oneLine <- readLines(f1, n = 1, warn = FALSE)) > 0) {
  oneLine<- unlist(strsplit(oneLine, split="\t"))
  oneLine<-unlist(strsplit(oneLine,split=" "))

  # Replace missing data code from "-" to "N"
  oneLine[oneLine == "-"]<- "N"
  
  # Replace ambiguous bases
  oneLine[oneLine == "R"] = sample(c("A","G"), sum(oneLine[-c(1,2)] == "R"),replace=T)
  oneLine[oneLine == "Y"] = sample(c("C","T"), sum(oneLine[-c(1,2)] == "Y"),replace=T)
  oneLine[oneLine == "W"] = sample(c("A","T"), sum(oneLine[-c(1,2)] == "W"),replace=T)
  oneLine[oneLine == "S"] = sample(c("C","G"), sum(oneLine[-c(1,2)] == "S"),replace=T)
  oneLine[oneLine == "M"] = sample(c("A","C"), sum(oneLine[-c(1,2)] == "M"),replace=T)
  oneLine[oneLine == "K"] = sample(c("T","G"), sum(oneLine[-c(1,2)] == "K"),replace=T)
  oneLine[oneLine == "U"] = "T"
  
  #### Printing a counter for conversion ####
  print(paste("converted line",i))
  i <- i+1
  
  write.table(t(oneLine[-1]), file = outputFile, sep = ",", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
  
}

close(f1) 
