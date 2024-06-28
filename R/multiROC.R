multiROC <- function(X, D, ...) {
  UseMethod("multiROC")
}

multiROC.default <- function(X, D, method = c("lrm", "fixedLinear", "fixedQuadratic", "dynamicEmpirical", "dynamicMeisner", "kernelOptimal"),
                             formula.lrm = "D ~ X.1 + I(X.1^2) + X.2 + I(X.2^2) + I(X.1*X.2)", stepModel = TRUE,
                             methodLinear = c("coefLinear", "SuLiu", "PepeThompson", "logistic", "minmax"), coefLinear = rep(1,ncol(X)), 
                             coefQuadratic = c(1,1,0,1,1), K = 201, alpha = 0.5, approxh = 0.5, multiplier = 2,
                             kernelOptimal.H = c("Hbcv", "Hscv", "Hpi", "Hns", "Hlscv", "Hbcv.diag", "Hscv.diag", "Hpi.diag", "Hlscv.diag"),
                             eps = sqrt(.Machine$double.eps), verbose = FALSE, ...){

  if(is.null(ncol(X))) stop("X should be a matrix where each row is a subject.")
  if(nrow(X) != length(D)) stop("X should be a matrix with the number of rows equal to the length of D. Each row of X is a subject.")
    
  pos.NA <- unique(which(is.na(X), arr.ind = TRUE)[,"row"], which(is.na(D)))
  if(length(pos.NA) >= 1){
    warning("Data from", paste(length(pos.NA), "subjects will be removed due to NA's in X or D."))
    X <- X[-pos.NA,]; D <- D[-pos.NA]
  }
  
  X <- as.matrix(X)
  p <- ncol(X)
  if(p == 1) stop("multiROC() function is designed for multivariate markers, not univariate.")

  levels.names <- levels(as.factor(D))
  if(length(levels.names) != 2) stop("D should be a vector with two different levels.")
  controls <- matrix(split(X,D)[[levels.names[1]]], ncol = p)
  cases <- matrix(split(X,D)[[levels.names[2]]], ncol = p)
  D <- ifelse(as.factor(D)==levels.names[1], 0, 1); levels <- levels(as.factor(D))


  for(i in 1:p){
    assign(paste("X.",i, sep=""), X[,i])
  }

  method <- match.arg(method)
  if(p > 2 & method %in% c("fixedQuadratic", "dynamicEmpirical")) stop("fixedQuadratic and dynamicEmpirical methods are only available for bivariate markers.")
  
  if(method == "lrm"){

    if(is.null(formula.lrm)) stop("formula.lrm should be a correct formula to run glm as a character.")
    model <- glm(as.formula(formula.lrm), family=binomial(link='logit'))
    if(stepModel) invisible(capture.output(model <- step(model)))
    Z <- suppressWarnings(predict(model))
    if(length(unique(Z))==1){
      stop("The regression model does not take any variable as a predictor. \n")
    }
    roc <- gROC(Z, D, side = 'right')

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, formula = formula.lrm, model = model, coefModel = model$coefficients, stepModel = stepModel, Z = Z, t = roc$t, roc = roc$roc, c = roc$c, auc = roc$auc)

  }else if(method == "fixedLinear"){

    methodLinear <- match.arg(methodLinear)

    if(p == 2){
      coefLinear <- if(methodLinear == "SuLiu") .suliu(controls, cases)$coef else 
        if(methodLinear == "PepeThompson") .nonpar.combine2.coef(controls, cases)[1:2] else 
          if(methodLinear == "logistic") .logistic(controls, cases)$coef else
            if(methodLinear == "minmax") .liu.coef(controls, cases)$coef else coefLinear
    }else{
      coefLinear <- if(methodLinear == "SuLiu") .suliu(controls, cases)$coef else
        if(methodLinear == "PepeThompson"){temp <- .pepethompson.multi(X, D); temp$final.beta[order(temp$final.var)]}else
          if(methodLinear == "logistic") .logistic(controls, cases)$coef else
            if(methodLinear == "minmax") .liu.coef(controls, cases)$coef else coefLinear
    }

    if(length(coefLinear)!=p) stop(paste0("coefLinear should be a vector with ", p, " elements indicating the coefficients of each marker in a linear combination coefLinear[1]*X.1 + ... + coefLinear[", p, "]*X.",p))
    Z <- c(tcrossprod(as.matrix(X), matrix(coefLinear, nrow = 1)))
    roc <- gROC(Z, D, side = 'right')

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, coefLinear = coefLinear, Z = Z, t = roc$t, roc = roc$roc, c = roc$c, auc = roc$auc)

  }else if(method == "fixedQuadratic"){ # Only for bivariate markers
    
    if(length(coefQuadratic) != 5) stop("coefQuadratic should be a vector with five elements indicating the coefficients of each marker in a quadratic combination coef[1]*X.1 + coef[2]*X.2 + coef[3]*(X.1*X.2) + coef[4]*(X.1)^2 + coef[5]*(X.2)^2.")
    
    x <- X[,1]; y <- X[,2]
    Z <- c(tcrossprod(cbind(x,y,x*y,x^2,y^2), matrix(coefQuadratic, nrow = 1)))
    roc <- gROC(Z, D, side = 'right')
    
    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, coefQuadratic = coefQuadratic, Z = Z, t = roc$t, roc = roc$roc, c = roc$c, auc = roc$auc)
    
  }else if (method %in% c("dynamicMeisner", "dynamicEmpirical")){

    m <- nrow(controls)
    t <- seq(0,1,1/m)
    N <- nrow(X)
    t <- t[-1]

    if(p == 2){
      x <- X[,1]; y <- X[,2]
      lx <- (max(x)-min(x))/20; ly <- (max(y)-min(y))/20
      xx <- seq(min(x)-lx, max(x)+lx, length.out=100)
      yy <- seq(min(y)-ly, max(y)+ly, length.out=100)
    }
    
    if(method == "dynamicMeisner"){
      
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
        output.CoefTable <- list(coef = coef, c = c, t = ti, roc = Se, Z = Z, lf = lf)
        if(p == 2){
          f <- outer(xx, yy, function(x,y) c(tcrossprod(cbind(x,y), matrix(multiroc$coefLinear, nrow = 1))))
          output.CoefTable$f <- f
        }
        output.CoefTable
      })
        if(verbose) close(bar)
    
    }else{  # method == "dynamicEmpirical" (only for bivariate markers)
      
      alpha <- seq(-1, 1, length = K)
      rate <- 1/alpha
      G <- cbind(c(rep(1,2*K), rep(-1,2*K)), rep(c(alpha, rate),2))
      
      if(verbose){
        cat("Progress bar: Estimation of the optimal linear combination for FPR using the empirical estimator\n"); flush.console()
        bar <- txtProgressBar(min = 0, max = t[m], style = 3)
      }
      CoefTable <- lapply(t, function(ti){
        if(verbose) setTxtProgressBar(bar, ti)
        C <- sapply(1:(4*K), function(j){
          Z <- G[j,1] * X[,1] + G[j,2] * X[,2]
          groc <- gROC(X = Z, D = D)
          index.opt <- which.min(groc$c[groc$t <= ti])
          c(tpr = groc$roc[index.opt], c = groc$c[index.opt], g = G[j,])
        })
        index.maxTPR <- which.max(C[1,])
        coef <- C[3:4,index.maxTPR]
        multiroc <- multiROC(X, D, method = "fixedLinear", coefLinear = coef)
        Z <- multiroc$Z
        lf <- max(min(sort(unique(Z))[-1]-sort(unique(Z))[-length(unique(Z))]), eps)
        f <- outer(xx,yy,function(x,y) c(tcrossprod(cbind(x,y), matrix(multiroc$coefLinear, nrow = 1))))
        list(coef = coef, c = C[2,index.maxTPR], t = ti, roc = C[1,index.maxTPR], Z = Z, f = f, lf = lf)
      })
      if(verbose) close(bar)
      
    }

    coef <- sapply(1:m, function(i){CoefTable[[i]]$coef})
    c <- sapply(1:m, function(i){CoefTable[[i]]$c})
    roc <- sapply(1:m, function(i){CoefTable[[i]]$roc})
    Z <- sapply(1:m, function(i){CoefTable[[i]]$Z})
    auc <- as.numeric(sum(roc[-m]*(t[-1] - t[-m])) + t[1]*roc[1])

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, coefLinear = coef, Z = Z, CoefTable = CoefTable, t = t, roc = roc, c = c, auc = auc)

  }else{ # if(method == "kernelOptimal")

    H.method <- match.arg(kernelOptimal.H)
    if(p > 2 & H.method %in% c("Hbcv","Hbcv.diag")){
      H.method <- "Hpi"
      warning("Hpi method is used to estimate the bandwidth matrix H used in the ks::kde() function. Hbcv method is designed for bivariate data.")
      }

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
