FastEPRR_VAR <- function(varFilePath=NULL, winLength=NULL, stepLength=NULL, gapFilePath=NULL, varOutputFilePath=NULL) {
  
  ################# check varFilePath ###############  
  if(is.null(varFilePath) || is.na(varFilePath) || !file.exists(varFilePath)) {
    stop("Please input 'varFilePath' (The absolute path of input file)!")
  }
  
  tryCatch({
    if(!is.null(varFilePath)) {
      varFilePath <- as.character(varFilePath) 
    }    
  }, error=function(err) {
    stop("'varFilePath' parameter is not String!")
  }, warning=function(warn) {
    stop("'varFilePath' parameter is not String!")
  },finally={
    
  })  
  ################# check varFilePath ###############  
  
  
  ################# check winLength ###############
  if(is.null(winLength) || is.na(winLength)) {
    stop("Please input valid character value 'winLength' (The length of sliding-window) (kb)!")
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
  if(is.null(stepLength) || is.na(stepLength)) {
    stop("Please input valid character value 'stepLength' (The length of sliding step) (kb))!")
  }
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
    tol <- 1e-5
    if(abs((winLength/stepLength)-2) > tol) {
      stop("'winLength' is not the double of 'stepLength'!")
    }
  }   
  ################# check stepLength ###############
  
  
  ################# check gapFilePath ###############
  if(!is.null(gapFilePath) && (!is.character(gapFilePath))) {
    stop("Please input valid absolute path of 'gapFilePath' (gap file path) or let it be NULL!")
  }
  if(!is.null(gapFilePath) && !file.exists(gapFilePath)) {
    stop("The 'gapFilePath' is not valid!")
  }
  ################# check gapFilePath ###############
  
  
  ################# check varOutputFilePath ###############
  if(is.null(varOutputFilePath) || is.na(varOutputFilePath) || !is.character(varOutputFilePath)) {
    stop("Please input valid character value 'varOutputFilePath' (The absolute path of output file)!")
  }
  if(!file.exists(dirname(varOutputFilePath))) {
    stop("The 'varOutputFilePath' is not valid!")
  }
  if(file.exists(varOutputFilePath)) {
    file.remove(varOutputFilePath)
  }
  ################# check varOutputFilePath ###############
  
  
  options(scipen=200)	
  
  rInfor <- readFile(varFilePath, winLength, stepLength, gapFilePath)
  
  posAndRhoList <- rInfor$parho
  
  len <- length(posAndRhoList)
  sc <- 1
  while(sc <= len) {
    rhosList <- obtainRho(posAndRhoList[[sc]], stepLength)
    saveFineRho(rhosList, varOutputFilePath)
    sc <- sc + 1
  }
}

saveFineRho <- function(rhosList, saveFile) {
  len <- length(rhosList)
  if(len > 1) {
    for(i in 1:len) {
      tmp <- rhosList[[i]]
      cat("Position(kb) ", posOutputFor(tmp$lp), "-", posOutputFor(tmp$rp), ":", "\n", sep="", file=saveFile, append=TRUE)
      tmpRhos <- tmp$r
      cat("Rho:", rhoOutputFor(mean(tmpRhos)), "\n", sep="", file=saveFile, append=TRUE)
    }
  }
}

obtainRho <- function(pRList, stepLen) {
  len <- length(pRList)
  cc <- 1
  nullcc <- 1
  delcc <- 1
  innBreak <- FALSE # delete last NULLs
  
  while(cc <= len) {
    if(is.null(pRList[[cc]]$r)) {
      nullcc <- cc+1
      
      if(nullcc > len) {# the last window is NULL
        delcc <- cc - 1
        while(delcc < length(pRList)) {
          pRList[[length(pRList)]] <- NULL
        }
        innBreak <- TRUE
      } else {
        while(is.null(pRList[[nullcc]]$r)) {
          nullcc <- nullcc + 1
          
          if(nullcc > len) {
            delcc <- cc - 1
            while(delcc < length(pRList)) {
              pRList[[length(pRList)]] <- NULL
            }
            innBreak <- TRUE
            break
          }
        }
      }
      
      if(innBreak) {
        break
      }
      
      avgRho <- (pRList[[cc-1]]$r + pRList[[nullcc]]$r)/2.0
      for(i in cc:(nullcc-1)) {
        pRList[[i]]$r <- avgRho
      }
      cc <- nullcc
    } else {
      cc <- cc + 1      
    }
  }
  
  
  comIndex <- 1
  len <- length(pRList)
  
  x1 <- 0
  x2 <- 0
  x3 <- 0
  x4 <- 0
  
  c2 <- 0
  c3 <- 0
  c4 <- 0
  
  L <- 0
  U <- 0
  A <- 0
  B <- 0
  C <- 0
  rhosList <- list()
  
  while((comIndex+2) <= len) {
    rho1 <- pRList[[comIndex]]$r
    rho2 <- pRList[[comIndex+1]]$r
    rho3 <- pRList[[comIndex+2]]$r
    

    if(rho2<(rho1+rho3)) {
      c2 <- rho1
      c3 <- rho2-rho1
      c4 <- rho3-rho2+rho1
      L <- max(0, -c3)
      U <- min(c2, c4)
      A <- c3-c2-c4
      B <- c2*c4 - c3*c4 - c3*c2
      C <- c2*c3*c4
      resL <- main(4, 3*A, 2*B, C)
      if(length(resL) == 2) { # only two solution
        r1 <- resL$x1
        r2 <- resL$x2
        
        if(((r1>L) && (r1<U)) && ((r2>L) && (r2<U))) { # r1 r2
          xVal <- c(L, U, r1, r2)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4), fx(r1, c2, c3, c4), fx(r2, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]
        } else if((r1>L) && (r1<U)) { # r1
          xVal <- c(L, U, r1)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4), fx(r1, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]    
        } else if((r2>L) && (r2<U)) { # r2
          xVal <- c(L, U, r2)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4), fx(r2, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]  
        } else {
          xVal <- c(L, U)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)] 
        }
        
        if(length(x1) > 1) {
          x1 <- sample(x1, 1)
        }
        
        x2 <- c2 - x1
        x3 <- c3 + x1
        x4 <- c4 - x1
        
      } else if(length(resL) == 3) { # have three solution
        r1 <- resL$x1
        r2 <- resL$x2
        r3 <- resL$x3
        
        if(((r1>L) && (r1<U)) && ((r2>L) && (r2<U)) && ((r3>L) && (r3<U))) { # r1 r2 r3
          xVal <- c(L, U, r1, r2, r3)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4), fx(r1, c2, c3, c4), fx(r2, c2, c3, c4), fx(r3, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]  
        } else if(((r1>L) && (r1<U)) && ((r2>L) && (r2<U))) { # r1 r2
          xVal <- c(L, U, r1, r2)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4), fx(r1, c2, c3, c4), fx(r2, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]  
        } else if(((r2>L) && (r2<U)) && ((r3>L) && (r3<U))) { # r2 r3
          xVal <- c(L, U, r2, r3)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4), fx(r2, c2, c3, c4), fx(r3, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]  
        } else if(((r1>L) && (r1<U)) && ((r3>L) && (r3<U))) { # r1 r3
          xVal <- c(L, U, r1, r3)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4), fx(r1, c2, c3, c4), fx(r3, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]  
        } else if((r1>L) && (r1<U)) { # r1
          xVal <- c(L, U, r1)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4), fx(r1, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]  
        } else if((r2>L) && (r2<U)) { # r2
          xVal <- c(L, U, r2)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4), fx(r2, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]  
        } else if((r3>L) && (r3<U)) { # r3
          xVal <- c(L, U, r3)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4), fx(r3, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]  
        } else {
          xVal <- c(L, U)
          yVal <- c(fx(L, c2, c3, c4), fx(U, c2, c3, c4))
          x1 <- xVal[yVal == max(yVal)]          
        }
        
        if(length(x1) > 1) {
          x1 <- sample(x1, 1)
        }
        
        x2 <- c2 - x1
        x3 <- c3 + x1
        x4 <- c4 - x1
      } 
    } else {
      x1 <- 0
      x4 <- 0
      tmp <- (rho2-rho1-rho3)/3
      x2 <- tmp + rho1
      x3 <- tmp + rho3
    }
    
    
    # store x1
    if(x1>=0) {
      x1lp <- pRList[[comIndex]]$lp
      rhosList <- storeRho(rhosList, x1lp, x1lp+stepLen-1, x1)    
    }
    
    # store x2
    if(x2>=0) {
      x2lp <- pRList[[comIndex+1]]$lp
      rhosList <- storeRho(rhosList, x2lp, x2lp+stepLen-1, x2)    
    }
    
    # store x3
    if(x3>=0) {
      x3lp <- pRList[[comIndex+2]]$lp
      rhosList <- storeRho(rhosList, x3lp, x3lp+stepLen-1, x3)  
    }
    
    # store x4
    if((x3>=0) && (x4>=0)) {
      x4lp <- x3lp+stepLen
      rhosList <- storeRho(rhosList, x4lp, x4lp+stepLen-1, x4)
    }
    
    comIndex <- comIndex + 1
  }

  return(rhosList)
}

storeRho <- function(rhosList, leftP, rightP, rho) {
  len <- length(rhosList)
  tol <- 1e-5
  if(len == 0) {
    rhosList[[1]] <- list(lp=leftP, rp=rightP, r=rho)
  } else {
    findIndex <- -1
    for(i in 1:len) {
      tmp <- rhosList[[i]]
      if((abs(tmp$lp-leftP)<=tol) && (abs(tmp$rp-rightP)<=tol)) {
        findIndex <- i
        break
      }
    }
    
    if(findIndex != -1) { # has find
      tmp <- rhosList[[findIndex]]
      tmpRhos <- tmp$r
      tmpRhos <- c(tmpRhos, rho)
      rhosList[[findIndex]] <- list(lp=leftP, rp=rightP, r=tmpRhos)
    } else { # not find
      rhosList[[len+1]] <- list(lp=leftP, rp=rightP, r=rho)
    }
  }
  
  return(rhosList)
}

fx <- function(x, c2, c3, c4) {
  return(x*(c2-x)*(c3+x)*(c4-x))
}

main <- function(a, b, c, d) {
  m <- 1.0*b/a
  n <- 1.0*c/a
  l <- 1.0*d/a
  p <- n-1.0*m*m/3
  q <- l-1.0*(n-1.0*m*m/3)*m/3-1.0*m*m*m/27
  delta1 <- q*q+4.0*p*p*p/27
  if(delta1>=0) {
    tmpPow <- delta1^(1.0/2)
    x1 <- sgn((q+tmpPow)/2)*(abs((q+tmpPow)/2)^(1.0/3))
    x2 <- sgn((q-tmpPow)/2)*(abs((q-tmpPow)/2)^(1.0/3))
    delta <- 6.0*x1*x2-3.0*x1*x1-3.0*x2*x2
    a <- (x1+x2)/2.0
    #b <- sqrt(-delta)/2.0
    #print(paste("the result is:", -2*a-(m/3),a-(m/3),"+",b,"i",a-(m/3),"-",b,"i"))
    return(list(x1=(-2*a-(m/3)), x2=(a-(m/3))))
  } else {
    a2Tmp <- a_2(q/2, sqrt(-delta1)/2)
    a1 <- a2Tmp$as
    b1 <- a2Tmp$bs
    delta <- 12*b*b
    x11 <- -2*a1
    x22 <- a1+sqrt(delta)/2
    x33 <- a1-sqrt(delta)/2
    #print(paste("x11=",x11-m/3,"x22=",x22-m/3,"x33=",x33-m/3))
    return(list(x1=(x11-m/3), x2=(x22-m/3), x3=(x33-m/3)))
  }

}


sgn <- function(x) {
  flag <- 0
  if(x>0) {
    flag <- 1
  } else if (x==0) {
    flag <- 0
  } else {
    flag <- -1
  }
  return(flag)
}

a_2 <- function(c, d) {
  if(c==0) {
    m <- pi/2
  } else {
    m <- (pi/2 - atan(d/c)) 
  }
  tmp <- (c*c+d*d)^(1.0/6)
  
  a <- tmp*cos(m/3)
  b <- tmp*sin(m/3)
  
  return(list(as=a, bs=b))
}





readFile <- function(filePath, winLength, stepLength, gapFilePath) {
  
  tol <- 1e-5
  forwPos <- NULL
  proWinLen <- NULL
  
  gapDF <- NULL
  if(!is.null(gapFilePath)) {
    tryCatch({
      gapDF <- read.table(gapFilePath, header=TRUE, colClasses="numeric", comment.char="")
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
  
  
  inFile <- file(filePath, open="r")
  posAndRhoList <- list()
  inBreak <- FALSE
  segCC <- 1
  
  while(length(lines<-readLines(inFile, n=3000, warn = FALSE)) > 0) {
    linesLength <- length(lines)
    for(i in 1:linesLength) {
      posLine <- lines[i]
      
      if(grepl("^Position.*", posLine)) {
        posInfo <- unlist(strsplit(posLine, " |-|:"))
        currPos <- as.numeric(as.character(posInfo[2]))*1000.0
        
        if(is.na(currPos)) {
          stop("Read file failed, please check your file !")
        }
        
        if(is.null(proWinLen)) {
          #proWinLen <- as.integer(as.character(as.numeric(as.character(posInfo[3]))*1000.0 - currPos + 1))
          proWinLen <- as.numeric(as.character(posInfo[3]))*1000.0 - currPos + 1
          
          if(abs(proWinLen-winLength) > tol) {
            stop("Your input 'winLength' is not equal to the window length in the input file !")
          } 
        }
        
        rhoLine <- lines[i+1]
        if(grepl("^Rho.*", rhoLine)) {
          rhoValue <- unlist(strsplit(rhoLine, " "))[1]
          rho <- as.numeric(substr(rhoValue, 5, nchar(rhoValue)))
          inBreak <- FALSE
          if(is.null(forwPos)) {
            posAndRhoList[[segCC]] <- list()
            posAndRhoList[[segCC]][[1]] <- segRho(currPos, currPos+winLength, rho)
          } else {
            if(abs((forwPos+stepLength) - currPos) <= tol) {
              
              tlen <- length(posAndRhoList[[segCC]]) + 1
              posAndRhoList[[segCC]][[tlen]] <- segRho(currPos, currPos+winLength, rho)
              
            } else if((forwPos + stepLength - currPos) > tol) {
              stop("The positions are not correct, please check !")
            } else {
              tempPos <- forwPos+stepLength
              while(tempPos < currPos) {
                # if [tempPos, tempPos+winLength] overlaps with gap, then break and start a new list
                if(insecGap(tempPos, tempPos+winLength, gapDF)) {
                  segCC <- segCC + 1
                  inBreak <- TRUE
                  break
                } else {  # else the value is NULL
                  tlen <- length(posAndRhoList[[segCC]]) + 1
                  posAndRhoList[[segCC]][[tlen]] <- segRho(tempPos, tempPos+winLength, NULL)
                  
                }
                tempPos <- tempPos+stepLength
              }
              
              if(inBreak) {
                posAndRhoList[[segCC]] <- list()  
              }
              
              tlen <- length(posAndRhoList[[segCC]]) + 1
              posAndRhoList[[segCC]][[tlen]] <- segRho(currPos, currPos+winLength, rho)
            }
          }
        }        
        forwPos <- currPos				
      }
    }
  }
  close(inFile)
  
  return(list(parho=posAndRhoList))
}

segRho <- function(leftP, rightP, rho) {
  return(list(lp=leftP, rp=rightP-1, r=rho))
}

