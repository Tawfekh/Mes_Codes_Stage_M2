FastEPRR_VCF_step1 <- function(vcfFilePath=NULL, erStart=NULL, erEnd=NULL, winLength=NULL, stepLength=NULL, winDXThreshold=30, idvlConsidered="all", idvlChrFormat="[1:1]", qualThreshold=20, gapFilePath=NULL, srcOutputFilePath=NULL) {
  
  #dyn.load("FastEPRR.dll")
  
  ################# check vcfFilePath ###############
  if(is.null(vcfFilePath) || is.na(vcfFilePath) || !is.character(vcfFilePath)) {
    stop("Please input valid character value 'vcfFilePath' (The absolute path of input file) (ends in [.gz|.vcf])!")
  }
  if(!file.exists(vcfFilePath) || !grepl(".*[/|\\].*(.gz|.vcf)$", vcfFilePath, ignore.case=TRUE)) {
    stop("Please input valid character value 'vcfFilePath' (The full path of input file) (ends in [.gz|.vcf])!")
  }
  ################# check vcfFilePath ###############
  
  
  
  ################# check erStart ###############
  options(scipen=200)
  if(!is.null(erStart) && !is.character(erStart)) {
    stop("Please input valid character value 'erStart' (The start position of estimating rho) (kb)!")
  }
  if(!is.null(erStart) && is.character(erStart)) {
    if(!grepl("^[0-9,.]+$", erStart)) {
      stop("Please input valid character value 'erStart' (The start position of estimating rho) (kb)!")
    } 
    erStart <- as.numeric(gsub(",", "", erStart))
    if(erStart == 0) {
      stop("Please input valid character value 'erStart' (The start position of estimating rho) (kb)!")
    }
    erStart <- erStart*1000
  }
  ################# check erStart ###############
  
  
  
  
  ################# check erEnd ###############
  if(!is.null(erEnd) && (!is.character(erEnd))) {
    stop("Please input valid character value 'erEnd' (The end position of estimating rho) (kb)!")
  }
  if(!is.null(erEnd) && is.character(erEnd)) {
    if(!grepl("^[0-9,.]+$", erEnd)) {
      stop("Please input valid character value 'erEnd' (The end position of estimating rho) (kb)!")
    } 
    erEnd <- as.numeric(gsub(",", "", erEnd))
    if(erEnd == 0) {
      stop("Please input valid character value 'erEnd' (The end position of estimating rho) (kb)!")
    }
    erEnd <- erEnd*1000
  }   
  if(is.numeric(erStart) && is.numeric(erEnd) && (erEnd <= erStart)) {
    stop("The 'erEnd' must be greater than 'erStart'!")
  }
  ################# check erEnd ###############
  
  
  
  
  ################# check winLength ###############
  if(is.null(winLength) || is.na(winLength)) {
    stop("Please input valid character value 'winLength' (The length of sliding-window) (kb)!")
  }
  if(is.null(winLength) && !is.null(stepLength)) {
    stop("Because you have specified 'stepLength', you must specify 'winLength' first!")
  }
  if(!is.null(winLength) && (!is.character(winLength))) {
    stop("Please input valid character value 'winLength' (The length of sliding-window) (kb)!")
  }
  if(!is.null(winLength) && is.character(winLength)) {
    if(!grepl("^[0-9,.]+$", winLength)) {
      stop("Please input valid character value 'winLength' (The length of sliding-window) (kb)!")
    } 
    winLength <- as.numeric(gsub(",", "", winLength))
    if(winLength == 0) {
      stop("Please input valid character value 'winLength' (The length of sliding-window) (kb)!")
    }
    winLength <- winLength*1000
  } 
  ################# check winLength ###############
  
  
  ################# check stepLength ###############
  if(!is.null(stepLength) && (!is.character(stepLength))) {
    stop("Please input valid character value 'stepLength' (The length of sliding step) (kb))!")
  }
  if(!is.null(stepLength) && is.character(stepLength)) {
    if(!grepl("^[0-9,.]+$", stepLength)) {
      stop("Please input valid character value 'stepLength' (The length of sliding step) (kb)!")
    } 
    stepLength <- as.numeric(gsub(",", "", stepLength))
    if(stepLength == 0) {
      stop("Please input valid character value 'stepLength' (The length of sliding step) (kb)!")
    }
    stepLength <- stepLength*1000
    if(stepLength > winLength) {
      warning("Because the 'stepLength' is larger than 'winLength', there is no overlapping sliding-windows!")
    }
  }   
  ################# check stepLength ###############
  
  
  
  
  ################# check winDXThreshold ###############
  if(!is.null(winDXThreshold) && !is.numeric(winDXThreshold)) {
    stop("Please input valid integer value 'winDXThreshold' (The minimum number of doubleton and xton of each window)!")
  }
  if(!is.null(winDXThreshold) && is.numeric(winDXThreshold)) {
    if(!grepl("^[0-9]+$", winDXThreshold)) {
      stop("Please input valid integer value 'winDXThreshold' (The minimum number of doubleton and xton of each window)!")
    } 
    winDXThreshold <- as.integer(winDXThreshold)
    if(winDXThreshold < 0) {
      stop("Please input valid integer value 'winDXThreshold' (The minimum number of doubleton and xton of each window)!")
    }
    if(winDXThreshold >= 0) {
      assign("winSNPThre", winDXThreshold, envir = .GlobalEnv)
    } else {
      assign("winSNPThre", NULL, envir = .GlobalEnv)
    }
  }
  ################# check winDXThreshold ###############
  
  
  
  
  ################# check idvlConsidered & idvlChrFormat###############
  specifiedIndiNames <- NULL
  allFormat <- NULL
  specifiedFormat <- NULL
  if(!is.null(idvlConsidered) && !is.character(idvlConsidered)) {
    stop("Please input valid character value 'idvlConsidered' (The considered individuals)!")
  }
  if(!is.null(idvlChrFormat) && !is.character(idvlChrFormat)) {
    stop("Please input valid character value 'idvlChrFormat' (The chromosome format of individual)!")
  }
  
  if(!is.null(idvlConsidered) && is.character(idvlConsidered)) {
    commaCount <- length(grep(":", idvlChrFormat))
    leftBCount <- length(grep("\\[", idvlChrFormat))
    rightBCount <- length(grep("\\]", idvlChrFormat))
    if((commaCount != leftBCount) || (leftBCount != rightBCount) || (commaCount != rightBCount)) {
      stop("Please input valid character value 'idvlChrFormat' (The chromosome format of individual)!")
    }
    if((commaCount != 1) || (leftBCount != 1) || (rightBCount != 1)) {
      stop("Please input valid character value 'idvlChrFormat' (The chromosome format of individual)!")
    } 
    firChr <- substr(idvlChrFormat, 2, 2)
    secChr <- substr(idvlChrFormat, 4, 4)
    if((firChr != "0") && (firChr != "1")) {
      stop("Please input valid character value 'idvlChrFormat' (The chromosome format of individual)!")
    }
    if((secChr != "0") && (secChr != "1")) {
      stop("Please input valid character value 'idvlChrFormat' (The chromosome format of individual)!")
    }
    
    if((firChr == "0") && (secChr == "0")) {
      stop("Please input valid character value 'idvlChrFormat' (The chromosome format of individual)!")
    } else if((firChr == "0") && (secChr == "1")) {
      allFormat <- 1
    } else if((firChr == "1") && (secChr == "0")) {
      allFormat <- -1
    } else if((firChr == "1") && (secChr == "1")) {
      allFormat <- 0
    }
    
    if(toupper(idvlConsidered) == toupper("all")) {
      
    } else {
      if(!grepl(";", idvlConsidered)) {
        stop("Please input valid character value 'idvlConsidered' (The considered individuals)!")
      }
      
      idvals <- unlist(strsplit(idvlConsidered, ";"))
      
      for(i in idvals) { 
        if(is.null(i) || is.na(i) || nchar(i) == 0) {
          stop("Please input valid character value 'idvlConsidered' (The considered individuals)!")
        }
      }
      
      specifiedIndiNames <- vector("character", length=length(idvals))
      
      if(grepl(":", idvlConsidered) && grepl("\\[", idvlConsidered) && grepl("\\]", idvlConsidered)) {
        hasUnifiedFormat <- FALSE
        commaCount <- length(grep(":", idvals))
        leftBCount <- length(grep("\\[", idvals))
        rightBCount <- length(grep("\\]", idvals))
        if((commaCount != leftBCount) || (leftBCount != rightBCount) || (commaCount != rightBCount)) {
          stop("Please input valid 'idvlConsidered' (The considered individuals)!")
        } 
        
        specifiedChrFor <- vector("integer", length=length(idvals))
        for(k in 1:length(idvals)) {
          currIndi <- idvals[k]
          if(grepl("\\[0:1\\]|\\[1:1\\]|\\[1:0\\]", currIndi)) {
            if(!grepl("\\[0:1\\]$|\\[1:1\\]$|\\[1:0\\]$", currIndi)) {
              stop(cat(currIndi, " is not valid!", sep=""))
            }
            commaCount <- length(grep(":", currIndi))
            leftBCount <- length(grep("\\[", currIndi))
            rightBCount <- length(grep("\\]", currIndi))
            if((commaCount != leftBCount) || (leftBCount != rightBCount) || (commaCount != rightBCount)) {
              stop(cat(currIndi, " is not valid!", sep=""))
            } 
            if((commaCount != 1) || (leftBCount != 1) || (rightBCount != 1)) {
              stop(cat(currIndi, " is not valid!", sep=""))
            } 
            forIndex <- regexpr("\\[0:1\\]|\\[1:1\\]|\\[1:0\\]", currIndi)[1]
            currFormat <- substr(currIndi, forIndex, forIndex+5)
            currIndiName <- substr(currIndi, 1, forIndex-1)
            if(nchar(currIndiName) <= 0) {
              stop(cat(currIndi, " is not valid!", sep=""))
            }
            specifiedIndiNames[k] <- trim(currIndiName)
            firChr <- substr(currFormat, 2, 2)
            secChr <- substr(currFormat, 4, 4)
            if((firChr != "0") && (firChr != "1")) {
              stop(cat(currIndi, " is not valid!", sep=""))
            }
            if((secChr != "0") && (secChr != "1")) {
              stop(cat(currIndi, " is not valid!", sep=""))
            }
            
            if((firChr == "0") && (secChr == "0")) {
              stop(cat(currIndi, " is not valid!", sep=""))
            } else if((firChr == "0") && (secChr == "1")) {
              specifiedChrFor[k] <- 1
            } else if((firChr == "1") && (secChr == "0")) {
              specifiedChrFor[k] <- -1
            } else if((firChr == "1") && (secChr == "1")) {
              specifiedChrFor[k] <- 0
            }
            
          } else {
            specifiedChrFor[k] <- NA
            specifiedIndiNames[k] <- trim(currIndi)
          }
        }
      } else { # unified format
        hasUnifiedFormat <- TRUE
        for(k in 1:length(idvals)) {
          specifiedIndiNames[k] <- trim(idvals[k])
        }
      }
      
      if(hasUnifiedFormat) {
        specifiedFormat <- rep(allFormat, length(idvals))
      } else {
        specifiedChrFor[which(is.na(specifiedChrFor))] <- allFormat
        specifiedFormat <- specifiedChrFor
      }
    }
  }
  ################# check idvlConsidered & idvlChrFormat###############
  
  
  ################# check qualThreshold ###############
  if(is.null(qualThreshold) || is.na(qualThreshold) || !is.numeric(qualThreshold) || !grepl("^[0-9]+$", qualThreshold)) {
    stop("Please input valid int value 'qualThreshold' (VCF quality threshold)!")
  } 
  ################# check qualThreshold ###############
  
  
  ################# check gapFilePath ###############
  if(!is.null(gapFilePath) && (!is.character(gapFilePath))) {
    stop("Please input valid absolute path of 'gapFilePath' (gap file path) or let it be NULL!")
  }
  if(!is.null(gapFilePath) && !file.exists(gapFilePath)) {
    stop("The 'gapFilePath' is not valid!")
  }
  ################# check gapFilePath ###############
  
  
  
  ################# check srcOutputFilePath ###############
  if(is.null(srcOutputFilePath) || is.na(srcOutputFilePath) || !is.character(srcOutputFilePath)) {
    stop("Please input valid character value 'srcOutputFilePath' (The absolute path of output file)!")
  }
  if(!file.exists(dirname(srcOutputFilePath))) {
    stop("The 'srcOutputFilePath' is not valid!")
  }
  if(file.exists(srcOutputFilePath)) {
    file.remove(srcOutputFilePath)
  }
  ################# check srcOutputFilePath ###############
  
  
  paras <- list(fp=vcfFilePath, es=erStart, ee=erEnd, wl=winLength, sl=stepLength, vcfis=specifiedIndiNames, vcfaf=allFormat, vcfsf=specifiedFormat, vcfqt=qualThreshold, gapfp=gapFilePath, sofp=srcOutputFilePath)
  
  headInfo <- c("chr nsam wdou wxton wH wHetero wavsk2 wavr2 misInfo startPos endPos")
  cat(headInfo, "\n", file=srcOutputFilePath, append=TRUE)
  
  
  assign("FastEPRRSavePaths", srcOutputFilePath, envir = .GlobalEnv)
  
  
  vcfEsRho(paras)  
}


winRhoVCF <- function(currChr, nsam, seqs, leftP, rightP) {
  
  mfsAndH <- getWinMFSAndH(nsam, seqs, ncol(seqs), "VCF")
  
  if(!is.null(mfsAndH)) {
    
    fourSSObs <- mfsAndH$fourSSInfo   
    mfs <- mfsAndH$mfsInfo
    
    if(mfsAndH$MDataOrnot) { # have missing data
      tempRow <- c(currChr, nsam, mfs[1], mfs[2], fourSSObs[1], fourSSObs[2], fourSSObs[3], fourSSObs[4], paste(c(mfsAndH$samSize, mfsAndH$MDataIDs), collapse=","), leftP, rightP)
    } else {
      tempRow <- c(currChr, nsam, mfs[1], mfs[2], fourSSObs[1], fourSSObs[2], fourSSObs[3], fourSSObs[4], NA, leftP, rightP)
    }
    
    cat(tempRow, "\n", file=FastEPRRSavePaths, append=TRUE)
    rm(mfsAndH)
  }
  
}

