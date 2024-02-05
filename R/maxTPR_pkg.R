# Functions provided by Meisner et al. (2017):
# A. Meisner, M. Carone, M. S. Pepe, and K. F. Kerr. Combining biomarkers by maximizing the true positive rate for a fixed false positive rate. 
# UW Biostatistics Working Paper Series (Working Paper 420), 2017. URL https://arxiv.org/abs/1910.02087.
# CRAN repository: https://cran.r-project.org/web/packages/maxTPR/index.html


### Function to maximize TPR

.stpr <- function(parms, Xmat, yvec, h){
  cutoff <- parms[1]
  betavec <- parms[-1]
  XDmat <- as.matrix(Xmat[yvec==1,])
  BxD <- XDmat %*% betavec
  rslt <- -sum(pnorm((BxD-cutoff)/h))/(nrow(XDmat))
  return(rslt)
}

.constrFPR <- function(parms, Xmat, yvec, h){
  cutoff <- parms[1]
  betavec <- parms[-1]
  XnDmat <- as.matrix(Xmat[yvec==0,])
  BxnD <- XnDmat %*% betavec
  rslt <- sum(pnorm((BxnD-cutoff)/h))/(nrow(XnDmat))
  return(rslt)
}

.constrnorm <- function(parms, Xmat, yvec, h){
  betavec <- parms[-1]
  norm(matrix(betavec), type="F")
}

.constrFPRcut <- function(betavec, cutoff, Xmat, yvec, h, tval){
  XnDmat <- as.matrix(Xmat[yvec==0,])
  BxnD <- XnDmat %*% betavec
  rslt <- sum(pnorm((BxnD-cutoff)/h))/(nrow(XnDmat))
  return(rslt-tval)
}

.maxTPR <- function(data, tval, initialval="rGLM", alpha=0.5, approxh=0.5, tolval=1e-4, stepsz=1e-5, multiplier=2){
  if(!is.data.frame(data)){
    stop("data must be a data.frame")
  }
  if(min(sapply(data, is.numeric)) != 1){
    stop("columns of data must be numeric")
  }
  yvals <- unique(data[,1])
  if((length(yvals) != 2) | (min(yvals) != 0) | (max(yvals) != 1)){
    stop("outcome (first column of data) must be 0/1")
  }
  data <- data[complete.cases(data),]
  names(data) <- c("D",paste("V",c(1:(ncol(data)-1)),sep=""))
  varnames <- names(data)
  preds <- varnames[-1] ### assumes that the outcome ("y") is the first column of the dataframe

  ##### Logistic regression
  glmmod <- glm(as.formula(paste(varnames[1], " ~ ", paste(preds, collapse=" + "), sep="")), family="binomial", data=data)
  glmcoef <- glmmod$coef[2:(length(preds)+1)]
  ##### Robust logistic regression
  #robustmod <- aucm::rlogit(as.formula(paste(varnames[1], " ~ ", paste(preds, collapse=" + "), sep="")), dat=data)
  robustmod <- suppressMessages(robustbase::BYlogreg(y = data[,varnames[1]], x0 = data[,preds])) # Modified by Sonia Perez-Fernandez because package `aucm` was removed from the CRAN repository.
  if(robustmod$convergence==TRUE){
    rglmcoef <- robustmod$coef[2:(length(preds)+1)]
  }else{
    rglmcoef <- glmcoef
  }

  normrglm = matrix(rglmcoef/norm(matrix(rglmcoef),type="F"),ncol=1) # have to normalize
  normglm = matrix(glmcoef/norm(matrix(glmcoef),type="F"),ncol=1)
  if(initialval=="rGLM"){
    beta0=normrglm
  }else{
    beta0=normglm
  }
  hval = sd(as.matrix(data[,-1]) %*% matrix(beta0,ncol=1))/(nrow(data)^(approxh))

  guess = quantile(as.matrix(data[data[,1]==0,-1]) %*% matrix(beta0,ncol=1), 1-tval, type=8)
  cutoffguess = uniroot(.constrFPRcut, c(guess-max(0.5,abs(multiplier*guess)), guess+max(0.5,abs(multiplier*guess))), betavec=beta0,
                        Xmat=data[,-1], yvec=data[,1], h=hval, tval=tval)$root
  nD = nrow(data) - sum(data[,1])

  maxlagrange <- Rsolnp::solnp(rbind(cutoffguess,beta0), .stpr, eqfun=.constrnorm, eqB=1, ineqfun=.constrFPR, ineqUB=tval+alpha/nD,
                               ineqLB=0, Xmat=data[,-1], yvec=data[,1], h=hval,
                               control=list("outer.iter"=10^3, "inner.iter"=10^4, "tol"=tolval, "delta"=stepsz))
  maxB <- maxlagrange$pars

  cutoffGLM = as.numeric(quantile(as.matrix(data[data[,1]==0,-1]) %*% matrix(normglm,ncol=1), 1-tval, type=8))
  cutoffrGLM = as.numeric(quantile(as.matrix(data[data[,1]==0,-1]) %*% matrix(normrglm,ncol=1), 1-tval, type=8))
  cutoffsTPRre = as.numeric(quantile(as.matrix(data[data[,1]==0,-1]) %*% matrix(maxB[-1],ncol=1), 1-tval, type=8))

  if(!(maxlagrange$convergence==0)){
    warning("sTPR algorithm failed to converge")
  }
  if(!(robustmod$convergence)){
    warning("rGLM algorithm failed to converge; estimates from GLM used instead")
  }

  return(list(sTPRrslt=c("delta"=maxB[1], "deltaRE"=cutoffsTPRre, "coef"=as.numeric(maxB[-1]), "convergence"=(maxlagrange$convergence==0)),
              rGLMrslt=c("delta"=cutoffrGLM, "coef"=as.numeric(normrglm), "convergence"=robustmod$convergence),
              GLMrslt=c("delta"=cutoffGLM, "coef"=as.numeric(normglm), "convergence"=glmmod$converged), Nobs=nrow(data)))
}
