#########################################################################################################################
##  Author:Abdou Rahmane WADE                                                                        Date:      2018   ##
##                                                                                                                     ##
##                           ######################################################################                    ##
##                           ##               Titre: Script for generate vcf file                ##                    ##
##                           ##           for ms simulations  with 25% missing dat               ##                    ##
##                           ######################################################################                    ##
#########################################################################################################################
args <- commandArgs(TRUE)

inputFile <- args[1]

nsam <- args[2]

nrep <- args[3]

#outputFile<- paste(args[1],".txt",sep="")

dir.create('./190cult_vcf')

dir.create('./25wild_vcf')

if (is.numeric(nsam)){
  cat(paste("nombre individus =", nsam, sep=" "))
} else {
    nsam <- as.numeric(nsam)
    cat(paste("nombre individus =", nsam, sep=" "))
  }
  


if (is.numeric(nrep)){
  cat(paste("nombre individus =", nrep, sep=" "))
} else{
  nrep <- as.numeric(nrep)
    cat(paste("nombre individus =", nrep, sep=" "))
  }
  



File  <- file(inputFile, open = "r")

oneLine1 <- readLines(File, n = 1, warn = FALSE)
oneLine2 <- readLines(File, n = 1, warn = FALSE)
readLines(File, n = 1, warn = FALSE)

write.table(oneLine1, file = paste("log",inputFile,sep="_"), sep = "\n", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
write.table(oneLine2, file = paste("log",inputFile,sep="_"), sep = "\n", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

i=1


while (i <= nrep) {
  
  first_line <- "##fileformat=VCFv4.1"
  
  write.table(first_line, file=paste("./25wild_vcf/outvcf_rep",i,".vcf",sep=""), sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)

  write.table(first_line, file=paste("./190cult_vcf/outvcf_rep",i,".vcf",sep=""), sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  
  second_line_wild <- t(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",(paste("Ind",rep(1:25),sep="_"))))
  
  write.table(second_line_wild, file=paste("./25wild_vcf/outvcf_rep",i,".vcf",sep=""), sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)

  second_line_cult <- t(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",(paste("Ind",rep(1:(nsam-25)),sep="_"))))
  
  write.table(second_line_cult, file=paste("./190cult_vcf/outvcf_rep",i,".vcf",sep=""), sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  
  readLines(File, n = 1, warn = FALSE)
  
  segsite <- readLines(File, n = 1, warn = FALSE)
  segsite <- unlist(strsplit(segsite,split=" "))
  segsite <- as.numeric(segsite[2])
  write.table(t(c("segsites",i,segsite)), file="segsites.txt", sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
  
  ref <- "A"
  
  alt <- "T"
  
  position <- readLines(File, n = 1, warn = FALSE)
  
  position <- unlist(strsplit(position,split=" "))
  
  position <- as.numeric(position[-1])
  
  position <- (position * 100000)
  
  mat <- matrix(ncol = nsam, nrow = segsite)
  l=1
  if(i<nrep){
    while ((oneLine <- readLines(File, n = 1, warn = FALSE)) != "") {
      oneLine <- unlist(strsplit(oneLine,split=""))
      for(j in 1:segsite){
        oneLine[j] <- paste(oneLine[j],"|",oneLine[j],sep = "")
      }
      mat[,l] <- (oneLine)
      l=l+1
    }
  } else {
    while (length(oneLine <-readLines(File, n = 1, warn = FALSE)) > 0) {
      oneLine <- unlist(strsplit(oneLine,split=""))
      for(j in 1:segsite){
        oneLine[j] <- paste(oneLine[j],"|",oneLine[j],sep = "")
      }
      mat[,l] <- (oneLine)
      l=l+1
    }
  }
	missx <- sample(c(1:segsite), (0.25*(nsam*segsite)), replace=TRUE)
	missy <- sample(c(1:nsam), (0.25*(nsam*segsite)), replace=TRUE)
	
	for(w in c(1: (0.25*(nsam*segsite)))){
	mat[missx[w],missy[w]] <- ".|."	
	}

for(k in 1:segsite){
      locus_linewild <- t(c("chr1",position[k],".",ref,alt,".",".",".","GT", mat[k,c(1:25)]))

	write.table(locus_linewild, file=paste("./25wild_vcf/outvcf_rep",i,".vcf",sep=""), sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)

      locus_line <- t(c("chr1",position[k],".",ref,alt,".",".",".","GT", mat[k,c(26:nsam)]))

	write.table(locus_line, file=paste("./190cult_vcf/outvcf_rep",i,".vcf",sep=""), sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
    }

  print(paste ("in process ", ((i/nrep)*100),"%",sep=" "))

  i=i+1
}
close(File) 
