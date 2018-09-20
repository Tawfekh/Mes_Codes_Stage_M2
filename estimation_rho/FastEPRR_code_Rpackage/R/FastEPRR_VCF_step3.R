FastEPRR_VCF_step3 <- function(srcFolderPath=NULL, DXFolderPath=NULL, finalOutputFolderPath=NULL) {
  
  
  ################# check srcFilePath ###############
  if(is.null(srcFolderPath) || is.na(srcFolderPath) || !is.character(srcFolderPath)) {
    stop("Please input valid string value 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1.")
  }
  if(!file.exists(srcFolderPath)) {
    stop("Please input valid string value 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1.")
  }
  ################# check srcFilePath ###############
  
  
  
  ################# check DXFolderPath ###############
  if(is.null(DXFolderPath) || is.na(DXFolderPath) || !is.character(DXFolderPath)) {
    stop("Please input valid string value 'DXFolderPath' (The folder path stores files which named by nsam_doulbeton_xton)!")
  }
  if(!file.exists(DXFolderPath)) {
    stop("The 'DXFolderPath' is not valid!")
  }
  ################# check DXFolderPath ###############
  
  
  
  ################# check finalOutputFolderPath ###############
  if(is.null(finalOutputFolderPath) || is.na(finalOutputFolderPath) || !is.character(finalOutputFolderPath)) {
    stop("Please input valid string value 'finalOutputFolderPath' (The folder path which stores the output files)!")
  }
  if(!file.exists(finalOutputFolderPath)) {
    stop("The 'finalOutputFolderPath' is not exist !")
  }
  ################# check finalOutputFolderPath ###############
  
  
  allChrs <- list.files(path=srcFolderPath, full.names=TRUE)
  if(length(allChrs) == 0) {
    stop("Please input valid string value 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1.")
  }
  
  alldouxData <- list.files(path=DXFolderPath, full.names=TRUE)
  if(length(alldouxData) == 0) {
    stop("Please input valid string value 'DXFolderPath' (The folder path stores files which named by nsam_doulbeton_xton)!")
  }
  
  myData <- NULL
  rho <- NA
  rhoL <- NA
  rhoR <- NA
  
  tryCatch({
    for(ele in allChrs) {
      temp <- read.table(file=ele, header=TRUE, stringsAsFactors=F) 
      temp <- temp[,c(1, 10, 11)]
      temp <- cbind(temp, rho, rhoL, rhoR)
      myData <- rbind(myData, temp)
    }
  }, error=function(err) {
    stop(paste("File", ele, "in", srcFolderPath, "is not valid!"))    
  }, warning=function(warn) {
    stop(paste("File", ele, "in", srcFolderPath, "is not valid!"))     
  },finally={
    
  }) 
  
  
  cat("Total number of DX files:", length(alldouxData), "\n")
  
  
  tryCatch({
    ddd <- 1
    for(ele in alldouxData) {
      cat("Current number of DX files:", ddd, "\n")
      temp <- read.table(file=ele, header=TRUE, stringsAsFactors=F)
      while(nrow(temp) > 0) {
        currRho <- temp[1, ]
        cnn <- currRho$chrName
        spp <- currRho$startPos
        epp <- currRho$endPos
        
        commons <- which((myData$chr == cnn) & (myData$startPos == spp) & (myData$endPos == epp))
        
        #if(length(commons) == 0) {
        #  stop("here is a bug!")  
       # }
        if(length(commons) > 0) {
          myData[commons, 4] <- currRho$rho
          myData[commons, 5] <- currRho$rhoL
          myData[commons, 6] <- currRho$rhoR
        }
        
        temp <- temp[-1, ]
      }
      ddd <- ddd +1
    }
  }, error=function(err) {
    stop(paste("File", ele, "in", DXFolderPath, "is not valid!"))    
  }, warning=function(warn) {
    stop(paste("File", ele, "in", DXFolderPath, "is not valid!"))     
  },finally={
    
  }) 
  
  
  while(nrow(myData) > 0) {
    currWin <- myData[1, ]
    if(!is.na(currWin$rho)) {
      currChr <- currWin$chr
      saveFile <- paste(finalOutputFolderPath, "/chr","_", currChr, sep="")
      cat("Position(kb) ", posOutputFor(currWin$startPos), "-", posOutputFor(currWin$endPos), ":", "\n", sep="", file=saveFile, append=TRUE)
      if((rhoOutputFor(currWin$rhoL) == "0.00") & (rhoOutputFor(currWin$rhoR) == "0.00")) {
        cat("Rho:", rhoOutputFor(currWin$rho), "\n", sep="", file=saveFile, append=TRUE)
      } else {
        cat("Rho:", rhoOutputFor(currWin$rho), " CIL:", rhoOutputFor(currWin$rhoL), " CIR:", rhoOutputFor(currWin$rhoR), "\n", sep="", file=saveFile, append=TRUE)
      }
    }
    myData <- myData[-1, ]
  }
  
  
}
