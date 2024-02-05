
predict.groc <- function(object, FPR = NULL, C = NULL, XL = NULL, XU = NULL, ...){

  obj <- object
  controls <- obj$controls; cases <- obj$cases; X <- c(controls, cases)
  t <- obj$t; roc <- obj$roc; side <- obj$side; param <- obj$param
  if(is.null(FPR) & is.null(C) & (is.null(XL) | is.null(XU))) FPR <- .15
  
  if(!is.null(FPR)){

    index <- which.max(t[t<=FPR]); T <- t[index]; ROC <- roc[index]

    if(side=='right' | side=='left'){
      C <- obj$c[index]
      if(side=='right') ClassSubsets <- c(C, Inf) else ClassSubsets <- c(-Inf, C)
    }else{
      XL <- obj$xl[index]; XU <- obj$xu[index]
      if(XL<min(X)) XL <- -Inf
      if(XU>max(X)) XU <- Inf
      if(side=='both') ClassSubsets <- rbind(c(-Inf, XL), c(XU, Inf)) else ClassSubsets <- c(XL, XU)
    }

  }else{

    if(side %in% c("left", "right") & !is.null(C)){
      if(side == "right") ClassSubsets <- c(C, Inf) else ClassSubsets <- c(-Inf, C)
      if(param){
        T <- pnorm(ClassSubsets[2], mean(controls), sd(controls)) - pnorm(ClassSubsets[1], mean(controls), sd(controls))
        ROC <- pnorm(ClassSubsets[2], mean(cases), sd(cases)) - pnorm(ClassSubsets[1], mean(cases), sd(cases))
      }else{
        T <- sum(controls > ClassSubsets[1] & controls < ClassSubsets[2])/length(controls)
        ROC <- sum(cases > ClassSubsets[1] & cases < ClassSubsets[2])/length(cases)
      }
    }

    if(side %in% c("both", "both2") & !is.null(XL) & !is.null(XU)){
      if(XL<min(X)) XL <- -Inf
      if(XU>max(X)) XU <- Inf
      if(side=='both'){
        ClassSubsets <- rbind(c(-Inf, XL), c(XU, Inf))
        if(param){
          T <- pnorm(ClassSubsets[1,2], mean(controls), sd(controls)) - pnorm(ClassSubsets[1,1], mean(controls), sd(controls)) + pnorm(ClassSubsets[2,2], mean(controls), sd(controls)) - pnorm(ClassSubsets[2,1], mean(controls), sd(controls))
          ROC <- pnorm(ClassSubsets[1,2], mean(cases), sd(cases)) - pnorm(ClassSubsets[1,1], mean(cases), sd(cases)) + pnorm(ClassSubsets[2,2], mean(cases), sd(cases)) - pnorm(ClassSubsets[2,1], mean(cases), sd(cases))
        }else{
          T <- sum((controls >= ClassSubsets[1,1] & controls <= ClassSubsets[1,2]) | (controls >= ClassSubsets[2,1] & controls <= ClassSubsets[2,2]))/length(controls)
          ROC <- sum((cases >= ClassSubsets[1,1] & cases <= ClassSubsets[1,2]) | (cases >= ClassSubsets[2,1] & cases <= ClassSubsets[2,2]))/length(cases)
        }
      }else{
        ClassSubsets <- c(XL, XU)
        if(param){
          T <- pnorm(ClassSubsets[2], mean(controls), sd(controls)) - pnorm(ClassSubsets[1], mean(controls), sd(controls))
          ROC <- pnorm(ClassSubsets[2], mean(cases), sd(cases)) - pnorm(ClassSubsets[1], mean(cases), sd(cases))
        }else{
          T <- sum(controls > ClassSubsets[1] & controls < ClassSubsets[2])/length(controls)
          ROC <- sum(cases > ClassSubsets[1] & cases < ClassSubsets[2])/length(cases)
        }
      }
    }

  }

  results <- list(ClassSubsets=ClassSubsets, Specificity = 1-T, Sensitivity = ROC)

  return(results)

}


predict.hroc <- function(object, FPR = 0.15, ...){

  obj <- object
  X <- obj$X; Y <- obj$Y; D <- obj$D; Sp <- obj$Sp; Se <- obj$Se; type <- obj$type
  indexX <- order(X)
  Xfun <- sort(unique(c(seq(min(X),max(X),length.out=1000),X)))

  C <- ifelse(type=='overfitting', 1-FPR, Y[which.min(ifelse(1-Sp <= FPR, FPR-1+Sp, 1))])
  if(min(ifelse(1-Sp <= FPR, FPR-1+Sp, 1)) == 1){C <- ifelse(type=='overfitting', 1.1, max(Y))}
  h <- approx(X[indexX], Y[indexX], xout=Xfun)$y
  Xcol <- Xfun[h > C]

  subsets <- function(Xcol){
    streaks <- rle(is.element(Xfun, Xcol))
    fin <- Xfun[cumsum(streaks$lengths)[streaks$values]]
    inicio <- Xfun[cumsum(streaks$lengths)[!streaks$values]+1]
    if(streaks$values[1]){
      inicio <- c(Xfun[1],inicio)
    }
    if(!streaks$values[length(streaks$values)]){
      inicio <- inicio[-length(inicio)]
    }
    inicio[inicio==min(X)] <- -Inf; fin[fin==max(X)] <- Inf
    cbind(inicio, fin)
  }

  ClassSubsets <- subsets(Xcol)
  colnames(ClassSubsets) <- c("inf", "sup")

  SP <- as.numeric(Sp[which.min(ifelse(1-Sp <= FPR, FPR-1+Sp, 1))])
  SE <- as.numeric(Se[which.min(ifelse(1-Sp <= FPR, FPR-1+Sp, 1))])

  if(FPR > 1-min(Sp)){ClassSubsets <- c(-Inf, Inf); SP <- 0; SE <- 1}

  results <- list(ClassSubsets=ClassSubsets, Specificity=SP, Sensitivity=SE)

  return(results)

}
