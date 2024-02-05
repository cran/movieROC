multiROC <- function(X, D, ...) {
  UseMethod("multiROC")
}

multiROC.default <- function(X, D, method = c("lrm", "fixedLinear", "dynamicMeisner", "kernelOptimal"),
                             formula.lrm = "D ~ X.1 + I(X.1^2) + X.2 + I(X.2^2) + I(X.1*X.2)", stepModel = TRUE,
                             coefLinear = rep(1,ncol(X)), methodLinear = c("coefLinear", "SuLiu", "PepeThompson", "logistic", "minmax"),
                             alpha = 0.5, approxh = 0.5, multiplier = 2,
                             kernelOptimal.H = c("Hpi", "Hscv", "Hns", "Hlscv", "Hscv.diag", "Hpi.diag", "Hlscv.diag"),
                             eps = sqrt(.Machine$double.eps), verbose = FALSE, ...){

  #X <- na.omit(X)
  X <- as.matrix(X)
  p <- ncol(X)

  levels.names <- levels(as.factor(D))
  controls <- matrix(split(X,D)[[levels.names[1]]], ncol = p)
  cases <- matrix(split(X,D)[[levels.names[2]]], ncol = p)
  D <- ifelse(as.factor(D)==levels.names[1], 0, 1); levels <- levels(as.factor(D))


  for(i in 1:p){
    assign(paste("X.",i, sep=""), X[,i])
  }

  method <- match.arg(method)

  if(method == "lrm"){

    if(is.null(formula.lrm)) stop("formula.lrm should be correct formula to run glm as a character.")
    model <- glm(as.formula(formula.lrm), family=binomial(link='logit'))
    if(stepModel) invisible(capture.output(model <- step(model)))
    Z <- suppressWarnings(predict(model))
    if(length(unique(Z))==1){
      stop("The regression model does not take any variable as a predictor. \n")
    }
    roc <- gROC(Z, D, side='right')

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, formula = formula.lrm, model = model, coefModel = model$coefficients, stepModel = stepModel, Z = Z, t = roc$t, roc = roc$roc, c = roc$c, auc = roc$auc)

  }else if(method == "fixedLinear"){

    methodLinear <- match.arg(methodLinear)

    coefLinear <- if(methodLinear == "SuLiu") .suliu(controls, cases)$coef else
      if(methodLinear == "PepeThompson"){temp <- .pepethompson.multi(X, D); temp$final.beta[order(temp$final.var)]}else
        if(methodLinear == "logistic") .logistic(controls, cases)$coef else
          if(methodLinear == "minmax") .liu.coef(controls, cases)$coef else coefLinear

    if(length(coefLinear)!=p) stop(paste0("coefLinear should be a vector with ", p, " elements indicating the coefficients of each marker in a linear combination coefLinear[1]*X.1 + ... + coefLinear[", p, "]*X.",p))
    Z <- c(tcrossprod(X, matrix(coefLinear, nrow = 1)))
    roc <- gROC(Z, D, side='right')

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, coefLinear = coefLinear, Z = Z, t = roc$t, roc = roc$roc, c = roc$c, auc = roc$auc)

  }else if (method == "dynamicMeisner"){

    m <- nrow(controls)
    t <- seq(0,1,1/m)
    N <- nrow(X)
    t <- t[-1]

    x <- X[,1]; y <- X[,2]
    lx <- (max(x)-min(x))/20; ly <- (max(y)-min(y))/20
    xx <- seq(min(x)-lx,max(x)+lx, length.out=100)
    yy <- seq(min(y)-ly,max(y)+ly, length.out=100)

    if(verbose){
      cat("Progress bar: Estimation of the optimal linear combination for FPR using maxTPR package\n"); flush.console()
      bar <- txtProgressBar(min = 0, max = t[m], style = 3)
    }
    CoefTable <- lapply(t, function(ti){
      if(verbose) setTxtProgressBar(bar, ti)
      invisible(capture.output(coef <- .maxTPR(as.data.frame(cbind(D, X)), tval = ti, alpha = alpha, approxh = approxh, multiplier = multiplier)$sTPRrslt[3:(p+2)]))
      multiroc <- multiROC(X, D, method = "fixedLinear", coefLinear = coef)
      indext <- which(abs(multiroc$t - ti) < .Machine$double.eps)
      c <- multiroc$c[indext]
      Se <- multiroc$roc[indext]
      Z <- multiroc$Z
      lf <- max(min(sort(unique(Z))[-1]-sort(unique(Z))[-length(unique(Z))]), eps)
      #f <- outer(xx,yy,function(x,y) c(tcrossprod(cbind(x,y), matrix(biroc$coefLinear, nrow = 1))))
      list(coef = coef, c = c, t = ti, roc = Se, Z = Z, lf = lf)
    })
    if(verbose) close(bar)

    coef <- sapply(1:m, function(i){CoefTable[[i]]$coef})
    c <- sapply(1:m, function(i){CoefTable[[i]]$c})
    roc <- sapply(1:m, function(i){CoefTable[[i]]$roc})
    Z <- sapply(1:m, function(i){CoefTable[[i]]$Z})
    auc <- as.numeric(sum(roc[-m]*(t[-1] - t[-m])) + t[1]*roc[1])

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, coefLinear = coef, Z = Z, CoefTable = CoefTable, t = t, roc = roc, c = c, auc = auc)

  }else{ # if(method == "kernelOptimal")

    H.method <- match.arg(kernelOptimal.H)

    f_cases <- ks::kde(cases, H = do.call(H.method, list(cases)))
    g_controls <- ks::kde(controls, H = do.call(H.method, list(controls)))

    optimalT <- function(x){
      f.x <- predict(f_cases, x = x); g.x <- predict(g_controls, x = x)
      f.x/(f.x + g.x)
    }
    Z <- optimalT(X)
    roc <- gROC(Z, D, side = "right")

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, H.method = H.method, optimalT = optimalT, Z = Z, t = roc$t, roc = roc$roc, c = roc$c, auc = roc$auc)

  }

  if(method == "fixedLinear") results$methodLinear <- methodLinear

  attr(results, 'class') <- 'multiroc'

  return(results)

}
