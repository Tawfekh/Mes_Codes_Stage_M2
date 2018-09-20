

missEstiRho <- function(mfsAndSS) {
  
  currNsam <- mfsAndSS$samSize
  mfs <- mfsAndSS$mfsInfo
  currDou <- mfs[1]
  currXto <- mfs[2]
  missDaIDs <- mfsAndSS$MDataIDs
  currH <- mfsAndSS$fourSSInfo[1]
  
  wZero <- getZeroLowHigh(currNsam, currDou, currXto, missDaIDs)
  zeroLow <- wZero[1]
  zeroHigh <- wZero[2]
  
  bmlAndH <- firstModel(currNsam, currDou, currXto, missDaIDs)
  
  thrH <- bmlAndH$thrH
  firBml <- bmlAndH$firMod
  secBml <- NULL
  
  if(currH <= zeroLow) {
    return(c(0.00, 0.00, 0.00))
  } else if(currH <= thrH) {
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
    return(c(returnRho, downCI, upCI))
  } else {
    if(is.null(secBml)) {
      secBml <- secModel(currNsam, currDou, currXto, missDaIDs, bmlAndH$firTD)
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
    return(c(returnRho, downCI, upCI))
  }
}

firstModel <- function(nsam, currDou, currXto, missDaIDs) {
  
  trainRho <- c(0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 70.0, 110.0, 170.0)
  
  nr <- 0L
  ssAndRho <- matrix(nrow=FastEPRRrepN*length(trainRho), ncol=5) 
  
  for(rho in trainRho) {
    cat(rho)
    seed1 <- sample(1:100000, 1)
    
    if(is.null(FastEPRRdemoParas)) {
      argv <- c("sim", as.character(nsam), as.character(FastEPRRrepN), "-Y", "0", currDou, currXto, "-U", "0", "-H", "1", as.character(10000), "-r", as.character(rho), "10000", "-seeds", as.character(seed1))
    } else {
      argv <- c("sim", as.character(nsam), as.character(FastEPRRrepN), "-Y", "0", currDou, currXto, "-U", "0", "-H", "1", as.character(10000), "-r", as.character(rho), "10000")
      argv <- c(argv, unlist(strsplit(FastEPRRdemoParas, " ")))
      argv <- c(argv, "-seeds", as.character(seed1))
    }
    
    RCallC <- .C("mySimulator", as.character(argv), as.integer(length(argv)), as.integer(vector("integer",length=1)), as.integer(missDaIDs), as.integer(length(missDaIDs)), fourSS=as.numeric(vector("numeric", length=FastEPRRrepN*4)), as.integer(5), PACKAGE="FastEPRR") #
    fourSSs <- RCallC$fourSS
    left <- 1 + nr*FastEPRRrepN
    nr <- nr + 1
    right <- nr*FastEPRRrepN
    fourSSMatrix <- matrix(fourSSs, ncol=4, byrow=T)
    ssAndRho[left:right, 1:4] <- fourSSMatrix
    ssAndRho[left:right, 5] <- rho
  }
  cat("\n")
  firMaxHs <- sort(fourSSMatrix[ ,1])
  firMaxH <- firMaxHs[length(firMaxHs)*0.95]
  firstBml <- boostingModel(ssAndRho)
  
  rm(fourSSMatrix)
  
  return(list(firMod=firstBml, thrH=firMaxH, firTD=ssAndRho))
}

secModel <- function(nsam, currDou, currXto, missDaIDs, firTD) {
  
  nr <- as.integer(nrow(firTD)/FastEPRRrepN)
  trainRhoSec <- c(180.0, 190.0, 200.0, 220.0, 250.0, 300.0, 350.0)
  ssAndRhoSec <- matrix(nrow=(nrow(firTD)+FastEPRRrepN*length(trainRhoSec)), ncol=5) 
  ssAndRhoSec[1:(nr*FastEPRRrepN),] <- firTD
  
  for(rho in trainRhoSec) {
    cat(rho)
    seed1 <- sample(1:100000, 1)
    if(is.null(FastEPRRdemoParas)) {
      argv <- c("sim", as.character(nsam), as.character(FastEPRRrepN), "-Y", "0", currDou, currXto, "-U", "0", "-H", "1", as.character(10000), "-r", as.character(rho), "10000", "-seeds", as.character(seed1))
    } else {
      argv <- c("sim", as.character(nsam), as.character(FastEPRRrepN), "-Y", "0", currDou, currXto, "-U", "0", "-H", "1", as.character(10000), "-r", as.character(rho), "10000")
      argv <- c(argv, unlist(strsplit(FastEPRRdemoParas, " ")))
      argv <- c(argv, "-seeds", as.character(seed1))
    }
    RCallC <- .C("mySimulator", as.character(argv), as.integer(length(argv)), as.integer(vector("integer",length=1)), as.integer(missDaIDs), as.integer(length(missDaIDs)), fourSS=as.numeric(vector("numeric", length=FastEPRRrepN*4)), as.integer(5), PACKAGE="FastEPRR") #
    fourSSs <- RCallC$fourSS
    left <- 1 + nr*FastEPRRrepN
    nr <- nr + 1
    right <- nr*FastEPRRrepN
    fourSSMatrix <- matrix(fourSSs, ncol=4, byrow=T)
    ssAndRhoSec[left:right, 1:4] <- fourSSMatrix
    ssAndRhoSec[left:right, 5] <- rho
  }
  cat("\n")
  secBml <- boostingModel(ssAndRhoSec)
  
  rm(ssAndRhoSec)
  rm(fourSSMatrix)
  
  return(secBml)
}

getZeroLowHigh <- function(nsam, douton, xton, missDaIDs) {
  seed1 <- sample(1:100000, 1)
  
  if(is.null(FastEPRRdemoParas)) {
    argv <- c("sim", as.character(nsam), as.character(1000), "-Y", "0", as.character(douton), as.character(xton), "-U", "0", "-H", "1", as.character(10000), "-r", as.character(0.0), "10000", "-seeds", as.character(seed1))
  } else {
    argv <- c("sim", as.character(nsam), as.character(1000), "-Y", "0", as.character(douton), as.character(xton), "-U", "0", "-H", "1", as.character(10000), "-r", as.character(0.0), "10000")
    argv <- c(argv, unlist(strsplit(FastEPRRdemoParas, " ")))
    argv <- c(argv, "-seeds", as.character(seed1))
  }
  
  RCallC <- .C("mySimulator", as.character(argv), as.integer(length(argv)), yiWanhaps = as.integer(vector("integer",length=1000)), as.integer(missDaIDs), as.integer(length(missDaIDs)), vector("numeric",length=1), as.integer(1), PACKAGE="FastEPRR") #
  
  zeroHs <- sort(RCallC$yiWanhaps)
  zeroLow <- zeroHs[50]
  zeroHigh <- zeroHs[950]
  rm(zeroHs)
  rm(RCallC)
  return(c(zeroLow, zeroHigh))
}

getCI <- function(mfsAndSS, mygam, pvRho, maxRho) {
  nsam <- mfsAndSS$samSize
  missDaIDs <- mfsAndSS$MDataIDs  
  
  mfs <- mfsAndSS$mfsInfo
  
  if(pvRho < 0.0) {
    pvRho <- 0.0
    return(c(0.0, 0.0))
  }
  
  if(pvRho > maxRho) {
    return(c(0.0, 0.0))
  }
  
  pvRho <- round(pvRho, 2)
  
  seed1 <- sample(1:100000, 1)
  
  if(is.null(FastEPRRdemoParas)) {
    argv <- c("sim", as.character(nsam), as.character(FastEPRRrepN), "-Y", "0", as.character(mfs[1]), as.character(mfs[2]), "-U", "0", "-H", "1", as.character(10000), "-r", as.character(pvRho), "10000", "-seeds", as.character(seed1))
  } else {
    argv <- c("sim", as.character(nsam), as.character(FastEPRRrepN), "-Y", "0", as.character(mfs[1]), as.character(mfs[2]), "-U", "0", "-H", "1", as.character(10000), "-r", as.character(pvRho), "10000")
    argv <- c(argv, unlist(strsplit(FastEPRRdemoParas, " ")))
    argv <- c(argv, "-seeds", as.character(seed1))
  }
  
  RCallC <- .C("mySimulator", as.character(argv), as.integer(length(argv)), as.integer(vector("integer",length=1)), as.integer(missDaIDs), as.integer(length(missDaIDs)), fourSS=as.numeric(vector("numeric", length=FastEPRRrepN*4)), as.integer(5), PACKAGE="FastEPRR") #
  fourSSs <- RCallC$fourSS
  fourSSMatrix <- matrix(fourSSs, ncol=4, byrow=T)
  fourSSDF <- data.frame(fourSSMatrix)      
  colnames(fourSSDF) <- c("H", "Hetero", "average_sk2", "average_r2")
  ordfourSSDF<- fourSSDF[with(fourSSDF, order(H)), ]
  rm(RCallC)
  rm(fourSSs)
  rm(fourSSMatrix)
  rm(fourSSDF)
  fourLower <- ordfourSSDF[3, ]
  ciLower <- predict(mygam, data.frame(H=fourLower[1], Hetero=fourLower[2], average_sk2=fourLower[3], average_r2=fourLower[4]))
  fourUpper <- ordfourSSDF[97, ]
  ciUpper <- predict(mygam, data.frame(H=fourUpper[1], Hetero=fourUpper[2], average_sk2=fourUpper[3], average_r2=fourUpper[4]))
  rm(ordfourSSDF)
  
  if((ciUpper >= maxRho) || (ciLower >= maxRho)) {
    ciUpper <- 0
    ciLower <- 0
  } 
  
  if((ciUpper <= pvRho) || (ciLower >= pvRho)) {
    ciUpper <- 0
    ciLower <- 0
  }   
  
  return(c(ciUpper, ciLower))
}



alphaCorrect <- function(mfsAndSS, mygam) {
  nsam <- mfsAndSS$samSize
  missDaIDs <- mfsAndSS$MDataIDs  
  
  mfs <- mfsAndSS$mfsInfo
  preVal <- 0.0
  alpha <- 0.0
  
  fourSSObs <- mfsAndSS$fourSSInfo
  preVal <- predict(mygam, data.frame(H=fourSSObs[1], Hetero=fourSSObs[2], average_sk2=fourSSObs[3], average_r2=fourSSObs[4]))
  
  if(preVal < 0.0) {
    preVal <- 0.0
    alpha <- 1.0
  } else {
    nearRho <- round(preVal, 2)
    seed1 <- sample(1:100000, 1)
    if(is.null(FastEPRRdemoParas)) {
      argv <- c("sim", as.character(nsam), as.character(FastEPRRrepN), "-Y", "0", as.character(mfs[1]), as.character(mfs[2]), "-U", "0", "-H", "1", as.character(10000), "-r", as.character(nearRho), "10000", "-seeds", as.character(seed1))
    } else {
      argv <- c("sim", as.character(nsam), as.character(FastEPRRrepN), "-Y", "0", as.character(mfs[1]), as.character(mfs[2]), "-U", "0", "-H", "1", as.character(10000), "-r", as.character(nearRho), "10000")
      argv <- c(argv, unlist(strsplit(FastEPRRdemoParas, " ")))
      argv <- c(argv, "-seeds", as.character(seed1))
    }
    RCallC <- .C("mySimulator", as.character(argv), as.integer(length(argv)), as.integer(vector("integer",length=1)), as.integer(missDaIDs), as.integer(length(missDaIDs)), fourSS=as.numeric(vector("numeric", length=FastEPRRrepN*4)), as.integer(5), PACKAGE="FastEPRR") #
    fourSSs <- RCallC$fourSS
    fourSSMatrix <- matrix(fourSSs, ncol=4, byrow=T)
    fourSSDF <- data.frame(fourSSMatrix)      
    colnames(fourSSDF) <- c("H", "Hetero", "average_sk2", "average_r2")
    nearRhoPvs <- predict(mygam, fourSSDF)
    meanPvs <- mean(nearRhoPvs)
    alpha <- nearRho/meanPvs
    rm(fourSSs)
    rm(fourSSMatrix)
    rm(fourSSDF)
    rm(nearRhoPvs)
  }
  
  rtnRho <- preVal*alpha
  
  if(rtnRho < 0.0) {
    rtnRho <- 0.0
  }
  
  return(rtnRho)
}

boostingModel <- function(SSAndRho) {
  mygam <- NULL
  
  myTrain <- data.frame(SSAndRho)
  colnames(myTrain) <- c("H", "Hetero", "average_sk2", "average_r2", "rho")
  myctrl <- boost_control(mstop=200)
  mygam <- gamboost(rho ~ ., data=myTrain, control=myctrl)     
  
  cv10f <- cv(model.weights(mygam), type = "kfold")
  cvm <- cvrisk(mygam, folds = cv10f, papply = lapply)
  mygam <- mygam[mstop(cvm)]
  
  rm(myTrain)
  rm(myctrl)
  rm(cv10f)
  rm(cvm)
  
  return(mygam)
}

convertToZeroOne <- function(mxSrc) {
  seqLength <- ncol(mxSrc)
  
  for(j in 1:seqLength) {
    tempSeq <- mxSrc[ , j] # one SNP
    tempTable <- sort(table(tempSeq))
    nuclCount <- as.vector(tempTable)
    nuclName <- names(tempTable)
    
    if('?' %in% nuclName) {
      quesMIndex <- which(nuclName == '?')
      nuclName <- nuclName[-quesMIndex]
      nuclCount <- nuclCount[-quesMIndex]
    }
    
    numOfnucl <- length(nuclName)
    if(numOfnucl == 2) {
      tempRow <- which(tempSeq == nuclName[1])
      mxSrc[tempRow, j] <- '1' 
      
      tempRow <- which(tempSeq == nuclName[2])
      mxSrc[tempRow, j] <- '0'
    } else if(numOfnucl == 3) { # multiple hit
      tempRow <- which(tempSeq == nuclName[1])
      mxSrc[tempRow, j] <- '1'
      
      tempRow <- which(tempSeq == nuclName[2])
      mxSrc[tempRow, j] <- '1'
      
      tempRow <- which(tempSeq == nuclName[3])
      mxSrc[tempRow, j] <- '0'
    } else if(numOfnucl == 4) { # multiple hit
      tempRow <- which(tempSeq == nuclName[1])
      mxSrc[tempRow, j] <- '1'
      
      tempRow <- which(tempSeq == nuclName[2])
      mxSrc[tempRow, j] <- '1'
      
      tempRow <- which(tempSeq == nuclName[3])
      mxSrc[tempRow, j] <- '1'
      
      tempRow <- which(tempSeq == nuclName[4])
      mxSrc[tempRow, j] <- '0'
    }
  }
  return(mxSrc)
}

trim <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
} 

installPkes <- function(packageName) {
  if(!(packageName %in% rownames(installed.packages()))) {
    cat("Installing R package '", packageName ,"'...", "\n\n", sep="")
    if(packageName == "mboost") {
      install.packages("mboost", repos="http://cran.r-project.org")      
    } else if(packageName == "seqinr") {
      install.packages("seqinr", repos="http://cran.r-project.org")
    }
  } 
  
  if(!require(packageName, character.only=TRUE)) {
    # for off-line user, they cannot use !
    stop(cat("Sorry,", "'", packageName , "' package is not available!", sep=""))
  }
}

posOutputFor <- function(x) { 
  x <- as.numeric(x)
  format(round(x/1000.0, 3), nsmall=3) 
}

rhoOutputFor <- function(x) { 
  x <- as.numeric(x)
  format(round(x, 2), nsmall=2) 
}





getWinMFSAndH <- function(nsam, seqMatrix, seqLength, flag) {
  mfs <- vector("numeric", length=2) 
  douOth <- vector("logical", length=seqLength)
  for(j in 1:seqLength) {
    tempSeq <- toupper(seqMatrix[ , j]) # one SNP
    tempTable <- sort(table(tempSeq))
    nuclCount <- as.vector(tempTable)
    nuclName <- names(tempTable)
    if(!('-' %in% nuclName) && !(hasOtherNucl(nuclName))) { # not gap
      
      if('?' %in% nuclName) {
        quesMIndex <- which(nuclName == '?')
        nuclName <- nuclName[-quesMIndex]
        nuclCount <- nuclCount[-quesMIndex]
      }
      
      numOfnucl <- length(nuclName)
      if(numOfnucl == 2) {
        minNum <- min(nuclCount)
        if(minNum == 2) { # doubleton
          mfs[1] <- mfs[1] + 1L
          douOth[j] = TRUE
        } else if(minNum >= 3){ # xton
          mfs[2] <- mfs[2] + 1L
          douOth[j] = TRUE
        }        
      } else if(numOfnucl == 3) { # multiple hit
        minNum <- nuclCount[1]
        if(minNum == 1) { # singleton mark '?'
          tempRow <- which(tempSeq == nuclName[1])
          seqMatrix[tempRow, j] <- '?'
        } else if(minNum == 2) { # doubleton
          mfs[1] <- mfs[1] + 1L
          douOth[j] = TRUE
        } else if(minNum >= 3){ # xton
          mfs[2] <- mfs[2] + 1L
          douOth[j] = TRUE
        }
        
        minNum <- nuclCount[2]
        if(minNum == 1) { # singleton mark '?'
          tempRow <- which(tempSeq == nuclName[2])
          seqMatrix[tempRow, j] <- '?'
        } else if(minNum == 2) { # doubleton
          mfs[1] <- mfs[1] + 1L
          douOth[j] = TRUE
        } else if(minNum >= 3){ # xton
          mfs[2] <- mfs[2] + 1L
          douOth[j] = TRUE
        }       
      } else if(numOfnucl == 4) { # multiple hit
        minNum <- nuclCount[1]
        if(minNum == 1) { # singleton mark '?'
          tempRow <- which(tempSeq == nuclName[1])
          seqMatrix[tempRow, j] <- '?'
        } else if(minNum == 2) { # doubleton
          mfs[1] <- mfs[1] + 1L
          douOth[j] = TRUE
        } else if(minNum >= 3){ # xton
          mfs[2] <- mfs[2] + 1L
          douOth[j] = TRUE
        }
        
        minNum <- nuclCount[2]
        if(minNum == 1) { # singleton mark '?'
          tempRow <- which(tempSeq == nuclName[2])
          seqMatrix[tempRow, j] <- '?'
        } else if(minNum == 2) { # doubleton
          mfs[1] <- mfs[1] + 1L
          douOth[j] = TRUE
        } else if(minNum >= 3){ # xton
          mfs[2] <- mfs[2] + 1L
          douOth[j] = TRUE
        }
        
        minNum <- nuclCount[3]
        if(minNum == 1) { # singleton mark '?'
          tempRow <- which(tempSeq == nuclName[3])
          seqMatrix[tempRow, j] <- '?'
        } else if(minNum == 2) { # doubleton
          mfs[1] <- mfs[1] + 1L
          douOth[j] = TRUE
        } else if(minNum >= 3){ # xton
          mfs[2] <- mfs[2] + 1L
          douOth[j] = TRUE
        }        
      }
    }
  }  
  
  if(is.null(winSNPThre)) {
    if((mfs[1]+mfs[2]) <= floor(log(nsam))) {
      return(NULL)
    }
  } else {
    if((mfs[1]+mfs[2]) < winSNPThre) {
      return(NULL)
    }
  }
  
  
  mxOutSingle <- seqMatrix[ ,douOth]
  
  rm(douOth)
  rm(seqMatrix)
  
  haveMissDa <- FALSE  
  
  missDaIDs <- which(mxOutSingle == '?')
  
  fiveSS <- NULL
  if(length(missDaIDs) > 0) {
    # decide whether all character in one sample are '?' 
    qTure <- (rowSums(mxOutSingle == '?') == ncol(mxOutSingle))
    qTureIndex <- which(qTure == TRUE)
    if(length(qTureIndex) > 0) {
      mxOutSingle <- mxOutSingle[-qTureIndex,]
      missDaIDs <- which(mxOutSingle == '?')
      nsam <- nsam - length(qTureIndex)
    }    
    missDaIDs <- missDaIDs-1
    haveMissDa <- TRUE      
  } 
  
  obserInfo <- cComuSS(nsam, mxOutSingle, length(missDaIDs), flag)
  
  
  rm(mxOutSingle)
  
  return(list(samSize=nsam, mfsInfo=mfs, HInfo=obserInfo$HInfo, fourSSInfo=obserInfo$fourSSInfo, MDataOrnot=haveMissDa, MDataIDs=missDaIDs))
}

cComuSS <- function(nsam, mxOutSingle, missLen, flag) {
  if(flag == "ALN") {
    mxOutSingle <- convertToZeroOne(mxOutSingle)
  } else {
    twoIds <- which(mxOutSingle == '2')
    if(length(twoIds) > 0) {
      mxOutSingle[twoIds] <- '1'
    }
  }  
  srcSeq <- apply(mxOutSingle, 1, paste, collapse="")
  RCallC <- .C("comSS", srcSeq, as.integer(nsam), as.integer(ncol(mxOutSingle)), as.integer(missLen), fourSS=as.numeric(vector("numeric", length=4)), PACKAGE="FastEPRR")
  fourSSs <- RCallC$fourSS
  return(list(samSize=nsam, HInfo=fourSSs[1], fourSSInfo=c(fourSSs[1], fourSSs[2], fourSSs[3], fourSSs[4])))
}

hasOtherNucl <- function(currNuclName) {
  fourNucl <- c("A", "T", "C", "G", "0", "1", "2", "?")    
  donotHave <- FALSE
  for(i in 1:length(currNuclName)) {
    if(!(currNuclName[i] %in% fourNucl)) {
      donotHave <- TRUE
      break
    }
  }
  return(donotHave)
}
