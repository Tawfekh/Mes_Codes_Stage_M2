
#########################################################################################################################
##  Author:Abdou Rahmane WADE                                                                        Date:     2018   ##
##                                                                                                                     ##
##                      ############################################################################                   ##
##                      ## Titre: Script for convert FastEPRR InputFile(vcf) from genotype format ##                   ##
##                      ##                           Description:                                 ##                   ##
##                      ############################################################################                   ##
##                                                                                                                     ##
##             From genotype to vcf format                                                                             ##
##             - haploid data                                                                                          ##
##             - per chromosome                                                                                        ##
##                                                                                                                     ##
##                                                                                                                     ##
##             - Exemple format fichier de sortie:                                                                     ##
##                                                                                                                     ##
##                CHR1 Position ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8 ind9 ind10 ind11 ind12                         ##
##                chr1 35914 G G G G G G G G G G G G                                                                   ##
##                chr1 35959 K G T G G T G T G K T G                                                                   ##
##                chr1 36170 - T T T T T T T T T T -                                                                   ##
##                chr1 36202 - A A A A A A A - A A A                                                                   ##
##                chr1 36206 T A A A A A A A - A A A                                                                   ##
##                chr1 36207 - C C C C C C C - C C C                                                                   ##
##                chr1 36226 W C A C C A C A - C C C                                                                   ##
##                chr1 36259 - T T - T T T R - T T T                                                                   ##
##                chr1 36284 R - - C C - C C - C C C                                                                   ##
##                                                                                                                     ##
##                                                                                                                     ##
##             - Exemple format fichier de sortie: vcf file input file for FastEPRR                                    ##
##                                                                                                                     ##
##                ##fileformat=VCFv4.1                                                                                 ##
##                #CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT  ind1 ind2 ind3 ind4 ind5 ind6                    ##
##                 ind7 ind8 ind9 ind10 ind11 ind12                                                                    ##
##                chr1  35959 . T A . . . GT  0|0 1|1 0|0 1|1 1|1 0|0 1|1 0|0 1|1 1|1 0|0 1|1                          ##
##                chr1  36170 . T A . . . GT  ?|? 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 0|0 ?|?                          ##
##                chr1  36202 . A T . . . GT  ?|? 0|0 0|0 0|0 0|0 0|0 0|0 0|0 ?|? 0|0 0|0 0|0                          ##
##                chr1  36206 . T A . . . GT  0|0 1|1 1|1 1|1 1|1 1|1 1|1 1|1 ?|? 1|1 1|1 1|1                          ##
##                chr1  36207 . C T . . . GT  ?|? 0|0 0|0 0|0 0|0 0|0 0|0 0|0 ?|? 0|0 0|0 0|0                          ##
##                chr1  36226 . T A . . . GT  0|0 1|1 1|1 1|1 1|1 1|1 1|1 1|1 ?|? 1|1 1|1 1|1                          ##
##                chr1  36259 . T A . . . GT  ?|? 0|0 0|0 ?|? 0|0 0|0 0|0 1|1 ?|? 0|0 0|0 0|0                          ##
##                chr1  36284 . A T . . . GT  0|0 ?|? ?|? 1|1 1|1 ?|? 1|1 1|1 ?|? 1|1 1|1 1|1                          ##
##                                                                                                                     ##
##                                                                                                                     ##
##                                                                                                                     ##
##            InputFile: each files in 190_pearl_millet_Haploid_WithSomeDiploidsLoci_data_lfmmb.tar.gz                 ##
##            OutputFile: each files in 20-04-18_190_pearl_millet_Haploid_vcfFile                                      ##
##                                                                                                                     ##
##        @ Script modifié à partir de vcf_from_genotype_format_Concetta (auteur: Concetta Burgarella 2016)            ##
##                                                                                                                     ##
#########################################################################################################################



args <- commandArgs(TRUE)

# Prepare vcf file (haploid data, ie sample size=N) chromosome by chromosome


input <- args[1]

### Write output file's name 

L <- length(unlist(strsplit(input,split = "")))

output<- paste(substring(input, 1 , L-4),".vcf",sep="")

print(input)
print(output)

### Write chromosome's number
chr <- substring(input, L-7 , L-4)

# write the firSt line of input file with pop names

first_line <- "##fileformat=VCFv4.1"
write.table(first_line, file=output, sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)


# open connection
genotypes  <- file(input, open = "r")

oneLine <- readLines(genotypes, n = 1, warn = FALSE)

Description.SampleID <- oneLine<-unlist(strsplit(oneLine, split=" "))

second_line <- t(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",as.character(Description.SampleID[3:length(Description.SampleID)])))
write.table(second_line, file=output, sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)



# read genetic data (column 3 is the reference genotype)




l <- 0


while (length(oneLine <- readLines(genotypes, n = 1, warn = FALSE)) > 0) {
  
  oneLine<-unlist(strsplit(oneLine,split=" "))
  
  l <- l+1
  ### Write chromosome and input line's numbers
  print(paste (chr, l, sep=" "))
  
  # random sampling of one allele in heterocygous loci
  
  oneLine[][oneLine[]== "R"] <- sample(c("G","A"), length(oneLine[][oneLine[]== "R"]), replace=T)
  oneLine[][oneLine[]== "Y"] <- sample(c("T","C"), length(oneLine[][oneLine[]== "Y"]), replace=T)
  oneLine[][oneLine[]== "K"] <- sample(c("G","T"), length(oneLine[][oneLine[]== "K"]), replace=T)
  oneLine[][oneLine[]== "M"] <- sample(c("A","C"), length(oneLine[][oneLine[]== "M"]), replace=T)
  oneLine[][oneLine[]== "S"] <- sample(c("G","C"), length(oneLine[][oneLine[]== "S"]), replace=T)
  oneLine[][oneLine[]== "W"] <- sample(c("A","T"), length(oneLine[][oneLine[]== "W"]), replace=T)
  
  
  # define the alternative allele (fictif)
  not_miss <- which(oneLine != "-")[-c(1:2)]
  
  if(oneLine[not_miss][1] == "A" || oneLine[not_miss][1] == "C" || oneLine[not_miss][1] == "G" ){
    alternative <- "T"
  }else{
    alternative <- "A" 
  }
  
  # define the reference allele
  ref <- oneLine[not_miss][1]
  
  
  # inds equal to reference
  oneLine[not_miss][oneLine[not_miss] == oneLine[not_miss][1]] <- "0|0"
  
  # inds different from reference
  same_ref <- which(oneLine == "0|0")
  not_ref <- setdiff(not_miss, same_ref)
  oneLine[not_ref] <- "1|1"

  #coding missing values as "?" instead of "."
  oneLine[][oneLine[]== "-"] <- "?|?"
  
  
  locus_line <- t(c(oneLine[1:2],".",ref,alternative,".",".",".","GT",oneLine[3:length(oneLine)]))
  
  write.table(locus_line, file=output, sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  
  
} 

write.table(paste (chr, l, sep=" "), file=paste ("nb _line",chr, sep="_"), sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)

close(genotypes)




