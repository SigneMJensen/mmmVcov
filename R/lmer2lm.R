function (lmerFit) 
{
  if (inherits(lmerFit, "lmerMod")) {
    yVec <- getME(lmerFit, "y")
    Zmat2 <- getME(lmerFit, "Zt")
    Vinv2 <- Reduce("+", mapply(function(x, y) {
      x * (y^2)
    }, lapply(getME(lmerFit, "Ztlist"), function(x) {
      t(x) %*% x
    }), lapply(VarCorr(lmerFit), function(x) {
      attr(x, "stddev")
    }))) + diag(attr(VarCorr(lmerFit), "sc")^2, length(yVec))
    Xmat2 <- getME(lmerFit, "X")
    naRows <- attr(lmerFit@frame, "na.action")
  }
  if (inherits(lmerFit, "lme")) {
    yVec <- getResponse(lmerFit)
    fullData <- as.data.frame(lmerFit[["data"]])
    if (inherits(lmerFit, "glmmPQL")) {
      ncolFD <- ncol(fullData)
      fullData <- fullData[, -c(ncolFD, ncolFD - 1, ncolFD - 
                                  2)]
    }
    Vinv2 <- bdiag((extract.lme.cov2(lmerFit, fullData))[["V"]])
    Xmat2 <- model.matrix(lmerFit, fullData)
    naRows <- as.vector(lmerFit[["na.action"]])
  }
  Vsqrt2 <- solve(chol(Vinv2))
  totalRows <- length(yVec) + length(naRows)
  XmatNew2 <- matrix(0, totalRows, ncol(Xmat2))
  colnames(XmatNew2) <- colnames(Xmat2)
  respNew2 <- rep(NA, totalRows)
  if (length(naRows) > 0) {
    XmatNew2[-naRows, ] <- as.matrix(t(Vsqrt2) %*% Xmat2)
    respNew2[-naRows] <- as.matrix(t(Vsqrt2) %*% yVec)
  }
  else {
    XmatNew2 <- as.matrix(t(Vsqrt2) %*% Xmat2)
    respNew2 <- as.matrix(t(Vsqrt2) %*% yVec)
  }
  X <- XmatNew2
  fit <- lm(respNew2 ~ X - 1, na.action = na.exclude)
  cnam <- colnames(fit$model$X)
  names(fit$coefficients) <- cnam
  colnames(fit$qr$qr) <- cnam
  names(fit$effects)[1:length(cnam)] <- cnam
  return(fit)
}