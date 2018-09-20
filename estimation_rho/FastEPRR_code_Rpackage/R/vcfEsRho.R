
vcfEsRho <- function(paras) {
  cat("Parsing file", paras$fp, "\n\n")
  erS <- paras$es
  erE <- paras$ee
  if(is.numeric(erS) && is.numeric(erE)) {
    readVCF(paras, erS, erE)
  } else if(is.numeric(erS) && !is.numeric(erE)) {
    readVCF(paras, erS, Inf)
  } else if(!is.numeric(erS) && is.numeric(erE)) {
    readVCF(paras, -Inf, erE)
  } else {
    readVCF(paras, -Inf, Inf)
  }    
}

readVCF <- function(paras, rStart, rEnd) {
  
  gapDF <- NULL
  if(!is.null(paras$gapfp)) {
    tryCatch({
      gapDF <- read.table(paras$gapfp, header=TRUE, colClasses="numeric", comment.char="")
      if(ncol(gapDF) > 2) {
        stop("Invalid gap file format !")
      }
    }, error=function(err) {
      stop("'gapFilePath' parameter is not String!")
    }, warning=function(warn) {
      stop("'gapFilePath' parameter is not String!")
    },finally={
      
    })  
  }
  
  if(is.null(paras$vcfis)) {
    readAllVCF(paras, rStart, rEnd, gapDF)
  } else {
    readSpecifiedVCF(paras, rStart, rEnd, gapDF)
  }
  rm(gapDF)
}

readSpecifiedVCF <- function(paras, rStart, rEnd, gapDF) {
  qual <- paras$vcfqt
  winL <- paras$wl
  if(is.numeric(paras$sl)) {
    stepL <- paras$sl
  } else {
    stepL <- winL
  }
  vcfCon <- file(paras$fp, open="r")
  CHROMIndex <- 0
  POSIndex <- 0
  REFIndex <- 0
  ALTIndex <- 0
  QUALIndex <- 0
  FORMATIndex <- 0
  indDF <- NULL
  totalIndiviNum <- 0
  allInds <- 0
  positions <- c()
  positionsCopy <- c()
  miniThres <- 0
  leftP <- rStart
  rightP <- leftP + winL - 1
  toEnd <- FALSE
  currPos <- 0
  currChr <- 0
  comindex <- NULL
  while(length(lines<-readLines(vcfCon, n=3000, warn = FALSE)) > 0) {
    linesLength <- length(lines)
    for(i in 1:linesLength) {
      tempLine <- lines[i]
      if(!grepl("^##.*", tempLine)) { 
        if(grepl("^#CHR.*", tempLine)) { # header
          headerInfo <- unlist(strsplit(tempLine, " |\t"))
          CHROMIndex <- grep("CHROM", headerInfo)
          POSIndex <- grep("POS", headerInfo)
          REFIndex <- grep("REF", headerInfo)
          ALTIndex <- grep("ALT", headerInfo)
          QUALIndex <- grep("QUAL", headerInfo)
          FORMATIndex <- grep("FORMAT", headerInfo)
          if((length(POSIndex)==0) || (length(FORMATIndex)==0) || (length(REFIndex)==0) || (length(ALTIndex)==0) || (length(QUALIndex)==0)) {
            stop("This VCF file has a wrong format, please check !\n")
          }
          indis <- paras$vcfis
          
          comIndi <- intersect(indis, headerInfo)
          if(length(comIndi) != length(indis)) {            
            stop(cat("There are some individuals not in VCF file !\n"))
          }
          
          specfor <- paras$vcfsf
          totalIndiviNum <- length(indis)
          for(ac in 1:totalIndiviNum) {
            tempFor <- specfor[ac]
            if(tempFor == 0) {
              allInds <- allInds + 2
            } else {
              allInds <- allInds + 1
            }
          }
          
          if(allInds < 6) {
            stop("The sample size required must be greater than or equal to six !")
          }
          
          if((winL != -1) && (winL <= floor(log(allInds)))) {
            stop(cat("The minimum value of winLength for sample size ", nsam, " is ", (floor(log(nsam))+1) ,"!", sep=""))
          }
          
          indDF <- data.frame(row.names=1:allInds)
          comindex <- which(headerInfo %in% comIndi)   
          miniThres <- floor(log(allInds))
          
        } else if(nchar(tempLine) > 0) {
          
          snpInfo <- unlist(strsplit(tempLine, " |\t"))
          #currChr <- as.integer(snpInfo[CHROMIndex])
		      currChr <- snpInfo[CHROMIndex]
          currPos <- as.numeric(snpInfo[POSIndex])
          refAllele <- snpInfo[REFIndex]
          altAllele <- snpInfo[ALTIndex]
          quality <- snpInfo[QUALIndex]

          
          if(leftP == -Inf) {
            leftP <- currPos
            rightP <- leftP + winL - 1
          }
          
          if((rEnd!=Inf) && (currPos>rEnd)) {
            interS <- intersect(leftP:rEnd, positionsCopy)
            interIndex <- which(positionsCopy %in% interS)
            if((length(positionsCopy[interIndex]) > miniThres) && (!insecGap(leftP, rEnd, gapDF))) {
              winRhoVCF(currChr, allInds, as.matrix(indDF[ ,interIndex]), leftP, rEnd)
            }
            toEnd <- TRUE
            break
          }
          
          if(currPos >= rStart) {
            passInfo <- FALSE
            if(quality == ".") {
              passInfo <- TRUE
            } else {
              if(as.numeric(quality)>qual) {
                passInfo <- TRUE
              } 
            }
            
          
            if(!(altAllele == ".") && (!isInsDel(refAllele, altAllele)) && passInfo) {
              individuals <- snpInfo[comindex]
              currCol <- c()
              hasBreak <- FALSE
              
              for(i in 1:totalIndiviNum) {
                currIndi <- individuals[i]
                currFor <- specfor[which(indis == headerInfo[comindex[i]])]
                if(currFor == 0) {
                  if(substr(currIndi, 2, 2) == "/") {
                    cat("This VCF file has unphased genotype data !\n")
                    hasBreak <- TRUE
                  }
                  firstAllele <- substr(currIndi, 1, 1)
                  checkAllel(firstAllele)
                  currCol <- c(currCol, firstAllele)
                  secondAllels <- substr(currIndi, 3, 3)
                  checkAllel(secondAllels)
                  currCol <- c(currCol, secondAllels)
                } else if(currFor == -1) {
                  if(substr(currIndi, 2, 2) == "/") {
                    cat("This VCF file has unphased genotype data !\n")
                    hasBreak <- TRUE
                  }
                  firstAllele <- substr(currIndi, 1, 1)
                  checkAllel(firstAllele)
                  currCol <- c(currCol, firstAllele)
                } else if(currFor == 1) {
                  if(substr(currIndi, 2, 2) == "/") {
                    cat("This VCF file has unphased genotype data !\n")
                    hasBreak <- TRUE
                  }
                  secondAllels <- substr(currIndi, 3, 3)
                  checkAllel(secondAllels)
                  currCol <- c(currCol, secondAllels)
                }
              }

              
              if(hasBreak == FALSE) {
                positionsCopy <- c(positionsCopy, currPos)
                indDF[ ,ncol(indDF)+1] <- currCol
              }
            }
          }
          
          if(currPos > rightP) {
            
            while((length(positionsCopy)>0) && (rightP <= positionsCopy[length(positionsCopy)])) {
              
              interS <- intersect(leftP:rightP, positionsCopy)
              interIndex <- which(positionsCopy %in% interS)
              
              if((length(positionsCopy[interIndex]) > miniThres) && (!insecGap(leftP, rightP, gapDF))) {  #
                winRhoVCF(currChr, allInds, as.matrix(indDF[ ,interIndex]), leftP, rightP)
              }
              
              revIndex <- (positionsCopy >= leftP)
              positionsCopy <- positionsCopy[revIndex]
              indDF <- as.data.frame(indDF[ ,revIndex])
              leftP <- leftP + stepL
              rightP <- rightP + stepL
              
            }
          }
        }
      }
    }
    if(toEnd) {
      break
    }
  }
  
  if((rEnd == Inf) || (rEnd == currPos)) {
    interS <- intersect(leftP:positionsCopy[length(positionsCopy)], positionsCopy)
    interIndex <- which(positionsCopy %in% interS)
    if((length(positionsCopy[interIndex]) > miniThres) && (!insecGap(leftP, positionsCopy[length(positionsCopy)], gapDF))) {
      winRhoVCF(currChr, allInds, as.matrix(indDF[ ,interIndex]), leftP, positionsCopy[length(positionsCopy)])
    }
  }
  
  close(vcfCon)
}

readAllVCF <- function(paras, rStart, rEnd, gapDF) {
  allfor <- paras$vcfaf
  qual <- paras$vcfqt
  winL <- paras$wl
  if(is.numeric(paras$sl)) {
    stepL <- paras$sl
  } else {
    stepL <- winL
  }
  vcfCon <- file(paras$fp, open="r")
  CHROMIndex <- 0
  POSIndex <- 0
  REFIndex <- 0
  ALTIndex <- 0
  QUALIndex <- 0
  FORMATIndex <- 0
  indDF <- NULL
  totalIndiviNum <- 0
  allInds <- 0
  positions <- c()
  positionsCopy <- c()
  miniThres <- 0
  leftP <- rStart
  rightP <- leftP + winL - 1
  toEnd <- FALSE
  currPos <- 0
  currChr <- 0
  while(length(lines<-readLines(vcfCon, n=3000, warn = FALSE)) > 0) {
    linesLength <- length(lines)
    for(i in 1:linesLength) {
      tempLine <- lines[i]
      if(!grepl("^##.*", tempLine)) { 
        if(grepl("^#CHR.*", tempLine)) { # header
          headerInfo <- unlist(strsplit(tempLine, " |\t"))
          CHROMIndex <- grep("CHROM", headerInfo)
          POSIndex <- grep("POS", headerInfo)
          REFIndex <- grep("REF", headerInfo)
          ALTIndex <- grep("ALT", headerInfo)
          QUALIndex <- grep("QUAL", headerInfo)
          FORMATIndex <- grep("FORMAT", headerInfo)
          if((length(POSIndex)==0) || (length(FORMATIndex)==0) || (length(REFIndex)==0) || (length(ALTIndex)==0) || (length(QUALIndex)==0)) {
            stop("This VCF file has a wrong format, please check !\n")
          }
          totalIndiviNum <- (length(headerInfo)-FORMATIndex) # init
          if(allfor == 0) {
            allInds <- totalIndiviNum*2
          } else {
            allInds <- totalIndiviNum
          }
          if(allInds < 6) {
            stop("The sample size required must be greater than or equal to six !")
          }
          if((winL != -1) && (winL <= floor(log(allInds)))) {
            stop(cat("The minimum value of winLength for sample size ", nsam, " is ", (floor(log(nsam))+1) ,"!", sep=""))
          }
          
          indDF <- data.frame(row.names=1:allInds)
          miniThres <- floor(log(allInds))
        } else if(nchar(tempLine) > 0) {
          snpInfo <- unlist(strsplit(tempLine, " |\t"))
          #currChr <- as.integer(snpInfo[CHROMIndex])
		      currChr <- snpInfo[CHROMIndex]
          currPos <- as.numeric(snpInfo[POSIndex])
          refAllele <- snpInfo[REFIndex]
          altAllele <- snpInfo[ALTIndex]
          quality <- snpInfo[QUALIndex]
          
          if(leftP == -Inf) {
            leftP <- currPos
            rightP <- leftP + winL - 1
          }
          
          if((rEnd!=Inf) && (currPos>rEnd)) {
            interS <- intersect(leftP:rEnd, positionsCopy)
            interIndex <- which(positionsCopy %in% interS)
            if((length(positionsCopy[interIndex]) > miniThres) && (!insecGap(leftP, rEnd, gapDF))) {
              winRhoVCF(currChr, allInds, indDF[, interIndex], leftP, rEnd)
            }
            toEnd <- TRUE
            break
          }
          
          if(currPos >= rStart) {
            passInfo <- FALSE
            if(quality == ".") {
              passInfo <- TRUE
            } else {
              if(as.numeric(quality)>qual) {
                passInfo <- TRUE
              } 
            }
            if(!(altAllele == ".") && (!isInsDel(refAllele, altAllele)) && passInfo) {
              individuals <- snpInfo[(FORMATIndex+1):length(snpInfo)]
              currCol <- c()
              hasBreak <- FALSE
              for(i in 1:totalIndiviNum) {
                currIndi <- individuals[i]
                if(allfor == 0) {
                  if(substr(currIndi, 2, 2) == "/") {
                    cat("This VCF file has unphased genotype data !\n")
                    hasBreak <- TRUE
                  }
                  
                  firstAllele <- substr(currIndi, 1, 1)
                  checkAllel(firstAllele)
                  currCol <- c(currCol, firstAllele)
                  
                  secondAllels <- substr(currIndi, 3, 3)
                  checkAllel(secondAllels)
                  currCol <- c(currCol, secondAllels)     
                } else if(allfor == -1) {
                  if(substr(currIndi, 2, 2) == "/") {
                    cat("This VCF file has unphased genotype data !\n")
                    hasBreak <- TRUE
                  }
                  firstAllele <- substr(currIndi, 1, 1)
                  checkAllel(firstAllele)
                  currCol <- c(currCol, firstAllele) 
                } else if(allfor == 1) {
                  if(substr(currIndi, 2, 2) == "/") {
                    cat("This VCF file has unphased genotype data !\n")
                    hasBreak <- TRUE
                  }
                  secondAllels <- substr(currIndi, 3, 3)
                  checkAllel(secondAllels)
                  currCol <- c(currCol, secondAllels)    
                }
              }
              if(hasBreak == FALSE) {
                positionsCopy <- c(positionsCopy, currPos)
                indDF[ ,ncol(indDF)+1] <- currCol
              }
            }
          }
          
          if(currPos > rightP) {
            while((length(positionsCopy)>0) && (rightP <= positionsCopy[length(positionsCopy)])) {
              interS <- intersect(leftP:rightP, positionsCopy)
              interIndex <- which(positionsCopy %in% interS)
              if((length(positionsCopy[interIndex]) > miniThres) && (!insecGap(leftP, rightP, gapDF))) {  #
                winRhoVCF(currChr, allInds, indDF[, interIndex], leftP, rightP)
              }
              revIndex <- (positionsCopy >= leftP)
              positionsCopy <- positionsCopy[revIndex]
			        indDF <- as.data.frame(indDF[, revIndex])
              leftP <- leftP + stepL
              rightP <- rightP + stepL
            }
          }
        }
      }
    }
    if(toEnd) {
      break
    }
  }
  
  if((rEnd == Inf) || (rEnd == currPos)) {
    interS <- intersect(leftP:positionsCopy[length(positionsCopy)], positionsCopy)
    interIndex <- which(positionsCopy %in% interS)
    if((length(positionsCopy[interIndex]) > miniThres) && (!insecGap(leftP, positionsCopy[length(positionsCopy)], gapDF))) {
      winRhoVCF(currChr, allInds, indDF[, interIndex], leftP, positionsCopy[length(positionsCopy)])
    }
  }

  close(vcfCon)
}

checkAllel <- function(firorSec) {
  if((firorSec == "0") || (firorSec == "1") || (firorSec == "2") || (firorSec == "?")) {
    
  } else {
    stop("There are some characters falling outside the range of '0', '1', '2' and '?'! \n")
  }
}

isInsDel <- function(refAllele, altAllele) {
  if(nchar(refAllele) > 1) {
    return(TRUE)
  }
  if(nchar(altAllele) > 1) {
    return(TRUE)
  }
  alts <- unlist(strsplit(altAllele, ","))
  for(i in 1:length(alts)) {
    if(nchar(alts[i]) > 1) {
      return(TRUE)
    }
  }
  return(FALSE)
}


insecGap <- function(leftP, rightP, gapDF) {
  hasInsec <- FALSE
  if(!is.null(gapDF)) {
    gapRow <- nrow(gapDF)
    for(i in 1:gapRow) {
      left <- gapDF[i, 1]
      right <- gapDF[i, 2]
      
      if((left<=rightP) && (leftP<=right)) {
        hasInsec <- TRUE
        break
      }
      
      if(left > rightP) {
        break
      }
    }
  }
  return(hasInsec)
}