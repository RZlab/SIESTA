fitSigmoid <- function(data, startPars = c("Pl" = 0, "a" = 550, "b" = 10), maxAttempts = 500){
  fctStr <- "(1 - Pl) * 1 / (1+exp(-(a/x-b))) + Pl"
  fitFct <- as.formula(paste("y ~", fctStr))
  xVec <- data[[1]]
  yVec <- data[[2]]
  varyPars <- 0

  attempts <- 0
  repeatLoop <- TRUE

  ## Check if number of non-missing values is sufficient
  ## (NLS can only handle data with at least three non-missing values)
  validValues <- !is.na(yVec)
  if (sum(validValues) <=2){
    m <- NA
    class(m) <- "try-error"
  } else{
    ## Perform fit
    while(repeatLoop & attempts < maxAttempts){
      parTmp <- startPars * (1 + varyPars*(runif(1, -0.2, 0.2)))
      m <- try(nls(formula=fitFct, start=parTmp, data=list(x=xVec, y=yVec), na.action=na.exclude,
                   algorithm="port", lower=c(0.0, 1e-5, 1e-5), upper=c(1.5, 15000, 250)),
               silent=TRUE)
      attempts <- attempts + 1
      varyPars <- 1
      if (class(m)!="try-error") {
        repeatLoop <- FALSE
      }
    }
  }
  return(m)
}
