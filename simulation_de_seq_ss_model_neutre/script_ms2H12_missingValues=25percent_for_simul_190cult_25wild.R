#########################################################################################################################
##  Author:Abdou Rahmane WADE                                                                        Date:      2018   ##
##                                                                                                                     ##
##                           ######################################################################                    ##
##                           ## Titre: Script for generate SelectionHapstats (Garud et al 2015)  ##                    ##
##                           ##           for ms simulations  with 25% missing dat               ##                    ##
##                           ######################################################################                    ##
#########################################################################################################################
args <- commandArgs(TRUE)

inputFile <- args[1]

nsam <- args[2]

nrep <- args[3]

#outputFile<- paste(args[1],".txt",sep="")

dir.create('./190cult')

dir.create('./25wild')

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
      oneLine[oneLine == "0"] <- ref
      oneLine[oneLine == "1"] <- alt
      mat[,l] <- (oneLine)
      l=l+1

    }
  } else {
    while (length(oneLine <-readLines(File, n = 1, warn = FALSE)) > 0) {
      oneLine <- unlist(strsplit(oneLine,split=""))*
      oneLine[oneLine == "0"] <- ref
      oneLine[oneLine == "1"] <- alt
      mat[,l] <- (oneLine)
      l=l+1
    }
  }
	missx <- sample(c(1:segsite), (0.25*(nsam*segsite)), replace=TRUE)
	missy <- sample(c(1:nsam), (0.25*(nsam*segsite)), replace=TRUE)
	
	for(w in c(1: (0.25*(nsam*segsite)))){
	mat[missx[w],missy[w]] <- "N"	
	}
      
for(k in 1:segsite){
      locus_line <- t(c(position[k],mat[k,c(1:25)]))
      write.table(locus_line, file=paste("./25wild/outH12_rep",i,".txt",sep=""), sep=",", quote=FALSE, col.names=F, row.names=F, append=TRUE)

      locus_line <- t(c(position[k],mat[k,c(26:nsam)]))
      write.table(locus_line, file=paste("./190cult/outH12_rep",i,".txt",sep=""), sep=",", quote=FALSE, col.names=F, row.names=F, append=TRUE)
      }

  print(paste ("in process ", ((i/nrep)*100),"%",sep=" "))

  i=i+1
}
close(File) 
