FastEPRR_VCF_step2 <- function(srcFolderPath=NULL, jobNumber=1, currJob=1, demoParameter=NULL, DXOutputFolderPath=NULL) {
  

  #dyn.load("FastEPRR.dll")
  
  assign("FastEPRRrepN", as.integer(100), envir = .GlobalEnv)
  
  ################# check srcFilePath ###############
  if(is.null(srcFolderPath) || is.na(srcFolderPath) || !is.character(srcFolderPath)) {
    stop("Please input valid character value 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1.")
  }
  if(!file.exists(srcFolderPath)) {
    stop("Please input valid character value 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1.")
  }
  ################# check srcFilePath ###############
  
  
  
  ################# check jobNumber  ###############
  tryCatch({
    if(is.null(jobNumber ) || is.na(jobNumber ) || !is.numeric(jobNumber )) {
      stop("Please input intger value 'jobNumber ' (The number of jobs)!")  
    }  
  }, error=function(err) {
    stop("Please input intger value 'jobNumber ' (The number of jobs)!")  
  }, warning=function(warn) {
    stop("Please input intger value 'jobNumber ' (The number of jobs)!")    
  },finally={
    
  }) 
  ################# check jobNumber  ###############
  
  
  
  
  ################# check currJob ###############
  tryCatch({
    if(is.null(currJob) || is.na(currJob) || !is.numeric(currJob)) {
      stop("Please input intger value 'currJob' (The current job number)!")  
    }  
  }, error=function(err) {
    stop("Please input intger value 'currJob' (The current job number)!")    
  }, warning=function(warn) {
    stop("Please input intger value 'currJob' (The current job number)!")     
  },finally={
    
  }) 
  if(currJob > jobNumber ) {
    stop("currJob must be less than or equal to jobNumber !")  
  }
  ################# check currJob ###############
  
  
  
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
  
  
  
  ################# check DXOutputFolderPath ###############
  if(is.null(DXOutputFolderPath) || is.na(DXOutputFolderPath) || !is.character(DXOutputFolderPath)) {
    stop("Please input valid character value 'DXOutputFolderPath' (The folder path stores files which named by DX_doubleton_xton)!")
  }
  if(!file.exists(DXOutputFolderPath)) {
    stop("The 'DXOutputFolderPath' is not valid!")
  }
  ################# check DXOutputFolderPath ###############
  
  
  allChrs <- list.files(path=srcFolderPath, full.names=TRUE)
  if(length(allChrs) == 0) {
    stop("Please input valid character value 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1.")
  }
  
  allSrcData <- NULL
  
  srcList <- list()

  tryCatch({
    for(ele in allChrs) {
      allSrcData <- rbind(allSrcData, read.table(file=ele, header=T, stringsAsFactors=F))
      srcList[[length(srcList)+1]] <- list(currFile=ele, currData=read.table(file=ele, header=TRUE, stringsAsFactors=F))
    } 
    
  }, error=function(err) {
    stop(paste("File", ele, "in", srcFolderPath,"is not valid!"))    
  }, warning=function(warn) {
    stop(paste("File", ele, "in", srcFolderPath, "is not valid!"))   
  },finally={
    
  }) 
  
  
  wuniMFS <- allSrcData[, c(2,3,4)]
  colnames(wuniMFS) <- c("nsam", "douton", "xton")
  uniqueMFS <- unique(wuniMFS)
  uniqRow <- nrow(uniqueMFS)
  
  rm(wuniMFS)
  rm(allSrcData)
  
  if(uniqRow < jobNumber) {
    stop(paste("The maximum value of jobNumber  is ", uniqRow, sep=""))
  }
  
  cat("Total job: ", jobNumber , "\n", "Current job: ", currJob, "\n", sep="")
  
  eachFileDXs <- floor(uniqRow/jobNumber)
  
  currDXs <- NULL
  left <- 1
  right <- eachFileDXs
  i <- 1
  while(i <= jobNumber) {
    if(i == jobNumber) {
      currDXs <- as.data.frame(uniqueMFS[left:uniqRow, ])
    } else if(i == currJob) {
      currDXs <- as.data.frame(uniqueMFS[left:right, ])
      break
    }
    left <- right + 1
    right <- right + eachFileDXs
    i <- i + 1 
  }
  
  installPkes("mboost")
  
  while(nrow(currDXs) > 0) {
    cat("DX size of currJob: ", nrow(currDXs), "\n", sep="")
    
    currMFS <- currDXs[1, ]
    currNsam <- as.integer(currMFS[1])
    currDou <- as.integer(currMFS[2])
    currXto <- as.integer(currMFS[3])    
    
    saveFile <- paste(DXOutputFolderPath, "/DX", "_", currDou, "_", currXto, sep="")
    
    if(file.exists(saveFile)) {
      file.remove(saveFile)
    }
    
    headInfo <- c("chrName startPos endPos rho rhoL rhoR")
    cat(headInfo, "\n", file=saveFile, append=TRUE)
    
    cat("Curr dton: ", currDou, " xton: ", currXto, "\n", sep="")
    
    if(length(srcList)>0) {
      sameConfig(currNsam, currDou, currXto, srcList, saveFile)
    }    
    currDXs <- currDXs[-1, ]
  }
}



sameConfig <- function(nsam, doubleton, xton, srcList, saveFile) {
  
  missDaIDs <- integer(0)
  wZero <- getZeroLowHigh(nsam, doubleton, xton, missDaIDs)
  zeroLow <- wZero[1]
  zeroHigh <- wZero[2]
  
  bmlAndH <- firstModel(nsam, doubleton, xton, missDaIDs)
  
  thrH <- bmlAndH$thrH
  firBml <- bmlAndH$firMod
  secBml <- NULL
  
  for(j in 1:length(srcList)) {
    tempData <- srcList[[j]]$currData
    commons <- which((tempData$wdou == doubleton) & (tempData$wxton == xton) & is.na(tempData$misInfo))
    for(i in commons) {
      currWin <- tempData[i, ]
      currH <- currWin$wH
      if(currH <= zeroLow) {
        cat(currWin$chr, currWin$startPos, currWin$endPos, rhoOutputFor(0.00), rhoOutputFor(0.00), rhoOutputFor(0.00), "\n", file=saveFile, append=TRUE)   
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
          cat(currWin$chr, currWin$startPos, currWin$endPos, rhoOutputFor(as.numeric(returnRho)), rhoOutputFor(0.00), rhoOutputFor(0.00), "\n", file=saveFile, append=TRUE)
        } else {
          cat(currWin$chr, currWin$startPos, currWin$endPos, rhoOutputFor(as.numeric(returnRho)), rhoOutputFor(as.numeric(downCI)), rhoOutputFor(as.numeric(upCI)), "\n", file=saveFile, append=TRUE)
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
          cat(currWin$chr, currWin$startPos, currWin$endPos, rhoOutputFor(as.numeric(returnRho)), rhoOutputFor(0.00), rhoOutputFor(0.00), "\n", file=saveFile, append=TRUE)
        } else {
          cat(currWin$chr, currWin$startPos, currWin$endPos, rhoOutputFor(as.numeric(returnRho)), rhoOutputFor(as.numeric(downCI)), rhoOutputFor(as.numeric(upCI)), "\n", file=saveFile, append=TRUE)
        }
        rm(mfsAndSS)
      }
    }
    
    # for missing data
    commons <- which((tempData$wdou == doubleton) & (tempData$wxton == xton) & !is.na(tempData$misInfo))
    for(i in commons) {
      currWin <- tempData[i, ]
      missAndsam <- as.integer(unlist(strsplit(currWin$misInfo, ",")))
      currNsam <- missAndsam[1]
      missDaIDs <- missAndsam[-1]
      mfsAndSS <- list(samSize=currNsam, mfsInfo=c(doubleton, xton), fourSSInfo=c(currWin$wH, currWin$wHetero, currWin$wavsk2, currWin$wavr2), MDataIDs=missDaIDs)
      erhos <- missEstiRho(mfsAndSS)
      if((rhoOutputFor(as.numeric(erhos[2])) == "0.00") & (rhoOutputFor(as.numeric(erhos[3])) == "0.00")) {
        cat(currWin$chr, currWin$startPos, currWin$endPos, rhoOutputFor(as.numeric(erhos[1])), rhoOutputFor(0.00), rhoOutputFor(0.00), "\n", file=saveFile, append=TRUE)
      } else {
        cat(currWin$chr, currWin$startPos, currWin$endPos, rhoOutputFor(as.numeric(erhos[1])), rhoOutputFor(as.numeric(erhos[2])), rhoOutputFor(as.numeric(erhos[3])), "\n", file=saveFile, append=TRUE)
      }
    }
  }
}
