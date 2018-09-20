
FastEPRR_ALN <- function(alnFilePath=NULL, format=-1, erStart=NULL, erEnd=NULL, winLength=NULL, stepLength=NULL, demoParameter=NULL, outputFilePath=NULL) {

  
  
 # dyn.load("FastEPRR.dll")
  ################# check alnFilePath ###############
  if(is.null(alnFilePath) || is.na(alnFilePath) || !is.character(alnFilePath)) {
    stop("Please input valid character value 'alnFilePath' (The absolute path of input file) (ends in [.fas|.aln|.phy])!")
  }
  if(!file.exists(alnFilePath) || !grepl(".*[/|\\].*(.fas|.aln|.phy)$", alnFilePath, ignore.case=TRUE)) {
    stop("Please input valid character value 'alnFilePath' (The full path of input file) (ends in [.fas|.aln|.phy])!")
  }
  ################# check alnFilePath ###############
  
  
  
  
  ################# check format ###############
  if(is.null(format) || is.na(format) || !grepl("^[1-3]$",format) || !is.numeric(format)) {
    stop("Please input int value 'format' (format number) of alignment file  (1:fasta; 2:clustal; 3:phylip)!")  
  }
  if((format == 1) && !grepl(".*[/|\\].*(.fas)$", alnFilePath)) {
    stop("Please input int value 'format' (format number) of alignment file  (1:fasta; 2:clustal; 3:phylip)!")  
  }
  if((format == 2) && !grepl(".*[/|\\].*(.aln)$", alnFilePath)) {
    stop("Please input int value 'format' (format number) of alignment file  (1:fasta; 2:clustal; 3:phylip)!")  
  }
  if((format == 3) && !grepl(".*[/|\\].*(.phy)$", alnFilePath)) {
    stop("Please input int value 'format' (format number) of alignment file  (1:fasta; 2:clustal; 3:phylip)!")  
  }
  ################# check format ###############
  
  
  
  
  ################# check erStart ###############
  options(scipen=200)
  if(!is.null(erStart) && !is.character(erStart)) {
    stop("Please input valid character value 'erStart' (The start position of estimating rho)!")
  }
  if(!is.null(erStart) && is.character(erStart)) {
    if(!grepl("^[0-9,]+$", erStart)) {
      stop("Please input valid character value 'erStart' (The start position of estimating rho)!")
    } 
    erStart <- as.integer(gsub(",", "", erStart))
    if(erStart < 1) {
      stop("The minimum value of erStart is 1!")
    }
  }
  ################# check erStart ###############
  
  
  
  
  ################# check erEnd ###############
  if(!is.null(erEnd) && (!is.character(erEnd))) {
    stop("Please input valid character value 'erEnd' (The end position of estimating rho)!")
  }
  if(!is.null(erEnd) && is.character(erEnd)) {
    if(!grepl("^[0-9,]+$", erEnd)) {
      stop("Please input valid character value 'erEnd' (The end position of estimating rho)!")
    } 
    erEnd <- as.integer(gsub(",", "", erEnd))
    if(erEnd < 1) {
      stop("Please input valid character value 'erEnd' (The end position of estimating rho)!")
    }
  }   
  if(is.integer(erStart) && (is.integer(erEnd)) && (erEnd <= erStart)) {
    stop("The 'erEnd' must be greater than 'erStart'!")
  }
  ################# check erEnd ###############
  

  
  
  ################# check winLength ###############
  if(is.null(winLength) && !is.null(stepLength)) {
    stop("Because you have specified 'stepLength', you must specify 'winLength' first!")
  }
  if(!is.null(winLength) && (!is.character(winLength))) {
    stop("Please input valid character value 'winLength' (The length of sliding-window)!")
  }
  if(!is.null(winLength) && is.character(winLength)) {
    if(!grepl("^[0-9,]+$", winLength)) {
      stop("Please input valid character value 'winLength' (The length of sliding-window)!")
    } 
    winLength <- as.integer(gsub(",", "", winLength))
    if(winLength == 0) {
      stop("Please input valid character value 'winLength' (The length of sliding-window)!")
    }
  }    
  ################# check winLength ###############

  
  
  
  ################# check stepLength ###############
  if(!is.null(stepLength) && (!is.character(stepLength))) {
    stop("Please input valid character value 'stepLength' (The length of sliding step) (kb))!")
  }
  if(!is.null(stepLength) && is.character(stepLength)) {
    if(!grepl("^[0-9,]+$", stepLength)) {
      stop("Please input valid character value 'stepLength' (The length of sliding step) (kb)!")
    } 
    stepLength <- as.integer(gsub(",", "", stepLength))
    if(stepLength == 0) {
      stop("Please input valid character value 'stepLength' (The length of sliding step) (kb)!")
    }
    if(stepLength > winLength) {
      warning("Because the 'stepLength' is larger than 'winLength', there is no overlapping sliding-windows!")
    }
  }   
  ################# check stepLength ###############
  

  ################# check demoParameter ###############
  if(!is.null(demoParameter) && (!is.character(demoParameter))) {
    stop("Please input valid character value 'demoParameter' (The population demographic parameters)!")
  }
  if(!is.null(demoParameter)) {
    demoInfo <- unlist(strsplit(demoParameter, " "))
    if((length(which(demoInfo == "-f")) > 0) || (length(which(demoInfo == "-seeds")) > 0) || (length(which(demoInfo == "-t")) > 0) ||
         (length(which(demoInfo == "-s")) > 0) || (length(which(demoInfo == "-T")) > 0) || (length(which(demoInfo == "-L")) > 0) ||
         (length(which(demoInfo == "-p")) > 0) || (length(which(demoInfo == "-r")) > 0) || (length(which(demoInfo == "-c")) > 0)) {
      stop("'-f -seeds -t -s -T -L -p -r -c' parameters are not allowed!")
    }
    if((length(which(demoInfo == "f")) > 0) || (length(which(demoInfo == "seeds")) > 0) || (length(which(demoInfo == "t")) > 0) ||
         (length(which(demoInfo == "s")) > 0) || (length(which(demoInfo == "T")) > 0) || (length(which(demoInfo == "L")) > 0) ||
         (length(which(demoInfo == "p")) > 0) || (length(which(demoInfo == "r")) > 0) || (length(which(demoInfo == "c")) > 0)) {
      stop("'-f -seeds -t -s -T -L -p -r -c' parameters are not allowed!")
    }
    assign("FastEPRRdemoParas", demoParameter, envir = .GlobalEnv)
  } else {
    assign("FastEPRRdemoParas", NULL, envir = .GlobalEnv)
  }
  ################# check demoParameter ###############
  
  
  
  ################# outputFilePath demoParameter ###############
  if(is.null(outputFilePath) || is.na(outputFilePath) || !is.character(outputFilePath)) {
    stop("Please input valid character value 'outputFilePath' (The absolute path of output file)!")
  }
  if(!file.exists(dirname(outputFilePath))) {
    stop("The 'outputFilePath' is not valid!")
  }
  if(file.exists(outputFilePath)) {
    file.remove(outputFilePath)
  }
  ################# outputFilePath demoParameter ###############
  

  installPkes("seqinr")    
  
  paras <- list(fp=alnFilePath, fm=format, es=erStart, ee=erEnd, wl=winLength, sl=stepLength)
  
  assign("winSNPThre", NULL, envir = .GlobalEnv)
  
  assign("FastEPRRrepN", 100, envir = .GlobalEnv)
  
  assign("FastEPRRDataInfos", data.frame(nsam=integer(), wdou=integer(), wxton=integer(), wH=integer(), wHetero=numeric(), wavsk2=numeric(), wavr2=numeric(),
                                         misInfo=character(), startPos=integer(), endPos=integer(), rho=numeric(), rhoL=numeric(), rhoR=numeric(), flag=integer(), stringsAsFactors=FALSE), envir = .GlobalEnv)
  alignEsRho(paras)
  startBoosting(outputFilePath)
}

startBoosting <- function(outputFilePath) {
  installPkes("mboost")
  
  while(nrow(FastEPRRDataInfos) > 0) {
    
    cat("window size:", nrow(FastEPRRDataInfos), "\n")
    currWin <- FastEPRRDataInfos[1, ]
    
    cat("Position ", currWin$startPos, "-", currWin$endPos, ":", "\n", sep="", file=outputFilePath, append=TRUE)
    
    
    currDou <- currWin$wdou
    currXto <- currWin$wxton
    currH <- currWin$wH
    
    douchkMis <- TRUE
    
    if(is.na(currWin$misInfo)) { 
      currNsam <- currWin$nsam
      missDaIDs <- integer(0) 
    } else { # have missing data
      missAndsam <- as.integer(unlist(strsplit(currWin$misInfo, ",")))
      currNsam <- missAndsam[1]
      missDaIDs <- missAndsam[-1]
      if(length(missDaIDs) == 0) {
        douchkMis <- FALSE
      }
    }
    
    if(is.na(currWin$flag)) {
      if((length(missDaIDs)==0) & douchkMis) {
        sameConfigALN(currNsam, currDou, currXto)
        currWin <- FastEPRRDataInfos[1, ]
        if((rhoOutputFor(as.numeric(currWin$rhoL)) == "0.00") & (rhoOutputFor(as.numeric(currWin$rhoR)) == "0.00")) {
          cat("Rho:", rhoOutputFor(as.numeric(currWin$rho)), "\n", sep="", file=outputFilePath, append=TRUE)
        } else {
          cat("Rho:", rhoOutputFor(as.numeric(currWin$rho)), " CIL:", rhoOutputFor(as.numeric(currWin$rhoL)), " CIR:", rhoOutputFor(as.numeric(currWin$rhoR)), "\n", sep="", file=outputFilePath, append=TRUE)
        }
      } else { # have missing data
        mfsAndSS <- list(samSize=currNsam, mfsInfo=c(currDou, currXto), fourSSInfo=c(currWin$wH, currWin$wHetero, currWin$wavsk2, currWin$wavr2), MDataIDs=missDaIDs)
        erhos <- missEstiRho(mfsAndSS)
        if((rhoOutputFor(as.numeric(erhos[2])) == "0.00") & (rhoOutputFor(as.numeric(erhos[3])) == "0.00")) {
          cat("Rho:", rhoOutputFor(as.numeric(erhos[1])), "\n", sep="", file=outputFilePath, append=TRUE)
        } else {
          cat("Rho:", rhoOutputFor(as.numeric(erhos[1])), " CIL:", rhoOutputFor(as.numeric(erhos[2])), " CIR:", rhoOutputFor(as.numeric(erhos[3])), "\n", sep="", file=outputFilePath, append=TRUE)
        }
      }
      assign("FastEPRRDataInfos", FastEPRRDataInfos[-1, ], envir = .GlobalEnv)
    } else {
      if((rhoOutputFor(as.numeric(currWin$rhoL)) == "0.00") & (rhoOutputFor(as.numeric(currWin$rhoR)) == "0.00")) {
        cat("Rho:", rhoOutputFor(as.numeric(currWin$rho)), "\n", sep="", file=outputFilePath, append=TRUE)
      } else {
        cat("Rho:", rhoOutputFor(as.numeric(currWin$rho)), " CIL:", rhoOutputFor(as.numeric(currWin$rhoL)), " CIR:", rhoOutputFor(as.numeric(currWin$rhoR)), "\n", sep="", file=outputFilePath, append=TRUE)
      }
      assign("FastEPRRDataInfos", FastEPRRDataInfos[-1, ], envir = .GlobalEnv)
    }
  }
}



sameConfigALN <- function(nsam, doubleton, xton) {
  missDaIDs <- integer(0)
  wZero <- getZeroLowHigh(nsam, doubleton, xton, missDaIDs)
  zeroLow <- wZero[1]
  zeroHigh <- wZero[2]
  
  bmlAndH <- firstModel(nsam, doubleton, xton, missDaIDs)
  
  thrH <- bmlAndH$thrH
  firBml <- bmlAndH$firMod
  secBml <- NULL
  
  commons <- which((FastEPRRDataInfos$wdou == doubleton) & (FastEPRRDataInfos$wxton == xton) & is.na(FastEPRRDataInfos$misInfo))
  for(i in commons) {
    currWin <- FastEPRRDataInfos[i, ]
    currH <- currWin$wH
    if(currH <= zeroLow) {
      FastEPRRDataInfos[i, 11] <<- 0.00
      FastEPRRDataInfos[i, 12] <<- 0.00
      FastEPRRDataInfos[i, 13] <<- 0.00
      FastEPRRDataInfos[i, 14] <<- 1
    } else if(currH <= thrH) {
      mfsAndSS <- list(samSize=currWin$nsam, mfsInfo=c(doubleton, xton), fourSSInfo=c(currWin$wH, currWin$wHetero, currWin$wavsk2, currWin$wavr2), MDataIDs=integer(0))
      returnRho <- alphaCorrect(mfsAndSS, firBml)
      
      CIul <- getCI(mfsAndSS, firBml, returnRho, 170.0)
      upCI <- max(CIul)
      if((currH >= zeroLow) && (currH <= zeroHigh)) {
        downCI <- 0.0
      } else {
        downCI <- min(CIul)
      }
      if(downCI < 0.0) {
        downCI <- 0
      }
      if((rhoOutputFor(as.numeric(upCI)) == "0.00") & (rhoOutputFor(as.numeric(downCI)) == "0.00")) {
        FastEPRRDataInfos[i, 11] <<- returnRho
        FastEPRRDataInfos[i, 12] <<- 0.00
        FastEPRRDataInfos[i, 13] <<- 0.00
        FastEPRRDataInfos[i, 14] <<- 1
      } else {
        FastEPRRDataInfos[i, 11] <<- returnRho
        FastEPRRDataInfos[i, 12] <<- downCI
        FastEPRRDataInfos[i, 13] <<- upCI
        FastEPRRDataInfos[i, 14] <<- 1
      }
      rm(mfsAndSS)
    } else {
      mfsAndSS <- list(samSize=currWin$nsam, mfsInfo=c(doubleton, xton), fourSSInfo=c(currWin$wH, currWin$wHetero, currWin$wavsk2, currWin$wavr2), MDataIDs=integer(0))
      
      if(is.null(secBml)) {
        secBml <- secModel(nsam, doubleton, xton, missDaIDs, bmlAndH$firTD)
      }
      
      returnRho <- alphaCorrect(mfsAndSS, secBml)
      
      CIul <- getCI(mfsAndSS, secBml, returnRho, 350.0)
      upCI <- max(CIul)
      if((currH >= zeroLow) && (currH <= zeroHigh)) {
        downCI <- 0.0
      } else {
        downCI <- min(CIul)
      }
      if(downCI < 0.0) {
        downCI <- 0
      }
      
      if((rhoOutputFor(as.numeric(upCI)) == "0.00") & (rhoOutputFor(as.numeric(downCI)) == "0.00")) {
        FastEPRRDataInfos[i, 11] <<- returnRho
        FastEPRRDataInfos[i, 12] <<- 0.00
        FastEPRRDataInfos[i, 13] <<- 0.00
        FastEPRRDataInfos[i, 14] <<- 1
      } else {
        FastEPRRDataInfos[i, 11] <<- returnRho
        FastEPRRDataInfos[i, 12] <<- downCI
        FastEPRRDataInfos[i, 13] <<- upCI
        FastEPRRDataInfos[i, 14] <<- 1
      }
      rm(mfsAndSS)
    }
  }
}

alnWinRho <- function(nsam, seqs, leftP, rightP) {
  
  #cat("pos:", leftP, "-" ,rightP, "\n")
  mfsAndH <- getWinMFSAndH(nsam, seqs, ncol(seqs), "ALN")
      
  if(!is.null(mfsAndH)) {
    fourSSObs <- mfsAndH$fourSSInfo   
    mfs <- mfsAndH$mfsInfo
    tempRow <- c(nsam, mfs[1], mfs[2], fourSSObs[1], fourSSObs[2], fourSSObs[3], fourSSObs[4])
    if(mfsAndH$MDataOrnot) { # have missing data
      tempRow <- c(tempRow, paste(c(mfsAndH$samSize, mfsAndH$MDataIDs), collapse=","), leftP, rightP, NA, NA, NA, NA)
    } else {
      tempRow <- c(tempRow, NA, leftP, rightP, NA, NA, NA, NA)
    }
  	dInfosRow <- nrow(FastEPRRDataInfos) + 1
    FastEPRRDataInfos[dInfosRow, ] <<- tempRow	
  }
}

