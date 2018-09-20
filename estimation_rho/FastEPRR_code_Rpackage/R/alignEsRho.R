
alignEsRho <- function(paras) {
  cat("Parsing file", paras$fp, "\n\n")
  srcData <- tryCatch({
    read.alignment(file=paras$fp, format=paras$fm, forceToLower = FALSE)
  }, error=function(err) {
    stop("Parse alignment file error, please check your file !")
  }, warning=function(warn) {
    stop("Parse alignment file error, please check your file !")
  },finally={
    
  })  
  sequences <- as.matrix(srcData$seq)
  nsam <- as.integer(srcData$nb)
  if(nsam < 6) {
    stop("The sample size required must be greater than or equal to six !")
  }
  

  if(is.integer(paras$wl) && (paras$wl <= floor(log(nsam)))) {
    stop(cat("The minimum value of winLength for sample size ", nsam, " is ", (floor(log(nsam))+1), sep=""))
  }
  
  seqMxSrc <- matrix(unlist(strsplit(sequences,"")), nrow = nsam, byrow=TRUE)  
  
  rm(sequences)
  
  seqLength <- ncol(seqMxSrc)
  
  if(seqLength <= 1) {
    stop("There are no SNP in your alignment file!")
  }
  
  erS <- paras$es
  erE <- paras$ee
  
  if(is.integer(erS) && (erS >= seqLength)) {
    stop(cat("The estimated rho Start position is larger than or equal to sequence length: ", seqLength, "!", sep=""))
  }
  
  if(is.integer(erE) && (erE >= seqLength)) {
    warning(cat("The estimated rho End position is larger than or equal to sequence length: ", seqLength, ". FastEPRR use sequence length as erEnd.", sep=""))
  }
  
  if(is.integer(erS) && is.integer(erE)) {
    pStart <- erS
    pEnd <- min(seqLength, erE)
  } else if(is.integer(erS) && !is.integer(erE)) {
    pStart <- erS
    pEnd <- seqLength
  } else if(!is.integer(erS) && is.integer(erE)) {
    pStart <- 1
    pEnd <- min(seqLength, erE)
  } else {
    pStart <- 1
    pEnd <- seqLength
  }
  
  seqMxPro <- seqMxSrc[ ,pStart:pEnd]
  rm(seqMxSrc)
  
  winL <- paras$wl
  stepL <- paras$sl
  if(!is.integer(winL) && !is.integer(stepL)) {
    alnWinRho(nsam, seqMxPro, pStart, pEnd)
  } else if(is.integer(winL) && !is.integer(stepL)) {
    oneSeqLen <- ncol(seqMxPro)
    if(winL >= oneSeqLen) {
      warning("Because the winLength is larger than or equal to alignments need to be estimated, there is no sliding-window!")
      alnWinRho(nsam, seqMxPro, pStart, pEnd)
    } else {
      leftP <- 1
      rightP <- winL
      pEnd <- pStart + winL - 1
      while(rightP <= oneSeqLen) {
        alnWinRho(nsam, seqMxPro[ ,leftP:rightP], pStart, pEnd)
        leftP <- leftP + winL
        rightP <- rightP + winL
        pStart <- pStart + winL
        pEnd <- pEnd + winL          
      }
      if((leftP+1) < oneSeqLen) { # last segment
        pEnd <- pStart + oneSeqLen - leftP
        alnWinRho(nsam, seqMxPro[ ,leftP:oneSeqLen], pStart, pEnd)
      }
    }
  } else if(is.integer(winL) && is.integer(stepL)) {
    oneSeqLen <- ncol(seqMxPro)
    if(winL >= oneSeqLen) {
      warning("Because the winLength is larger than or equal to alignments need to be estimated, there is no sliding-window!")
      alnWinRho(nsam, seqMxPro, pStart, pEnd)
    } else {
      leftP <- 1
      rightP <- winL
      pEnd <- pStart + winL - 1
      while(rightP <= oneSeqLen) {
        alnWinRho(nsam, seqMxPro[ ,leftP:rightP], pStart, pEnd)
        leftP <- leftP + stepL
        rightP <- rightP + stepL
        pStart <- pStart + stepL
        pEnd <- pEnd + stepL   
      }
      if((leftP+1) < oneSeqLen) { # last segment
        pEnd <- pStart + oneSeqLen - leftP 
        alnWinRho(nsam, seqMxPro[ ,leftP:oneSeqLen], pStart, pEnd)
      }
    }
  }
}

