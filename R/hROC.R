hROC <- function(X, D, ...) {
  UseMethod("hROC")
}

hROC.default <- function(X, D, type = c("lrm", "h.fun", "overfitting"), formula.lrm = "D ~ pol(X,3)", h.fun = function(x){x}, 
                         plot.h = FALSE, plot.roc = FALSE, new.window = FALSE, main = NULL, xlab = "x", ylab = "h(x)", xaxis = TRUE, ...){

  type <- match.arg(type)
  levels.names <- levels(as.factor(D))
  controls <- split(X,D)[[levels.names[1]]]; cases <- split(X,D)[[levels.names[2]]]
  D <- ifelse(as.factor(D)==levels.names[1], 0, 1); levels <- levels(as.factor(D))

  m <- sum(D==0); n <- sum(D==1)
  indexX <- order(X)

  if(type == 'lrm'){
    model <- rms::lrm(as.formula(formula.lrm))
    Y <- predict(model, data.frame(X = X), type = 'fitted')
    Y.controls <- predict(model, data.frame(X = controls), type = 'fitted')
    Y.cases <- predict(model, data.frame(X = cases), type = 'fitted')
  }else{
    if(type == 'overfitting'){
      h <- function(x, marker0, marker1){
        marker <- c(marker0, marker1)
        xi <- marker[which.min(abs(marker - x))]
        mxi <- sum(marker0 == xi); nxi <- sum(marker1 == xi)
        ifelse(mxi > nxi, 0, 1)
      }
    }else{
      h <- function(x, marker0, marker1) h.fun(x)
    }
    Y <- sapply(X, function(x){h(x, controls, cases)})
    Y.controls <- sapply(controls, function(x){h(x, controls, cases)})
    Y.cases <- sapply(cases, function(x){h(x, controls, cases)})
  }

  if(plot.h){
    if(new.window) dev.new(width = 6, height = 5)
    plot(X[indexX], Y[indexX], 'l', xlab = xlab, ylab = ylab, 
         main = ifelse(is.null(main), paste("Model:", type, ifelse(type=='lrm', formula.lrm, "")), main), 
         xaxt = ifelse(xaxis, "s", "n"), ...)
    # if(xaxis != FALSE){
    #   axis(1, at = seq(min(X), max(X), 0.05))
    #   axis(side = 1, at = seq(min(X), max(X), 0.005), tcl = -0.2, labels = FALSE)
    # }
    readline("Press return for next page....")
  }

  c <- Y
  if(type=='overfitting'){
    Sp <- sapply(c, function(c){sum(Y.controls < c)/m})
    Se <- sapply(c, function(c){sum(Y.cases >= c)/n})
  }else{
    Sp <- sapply(c, function(c){sum(Y.controls <= c)/m})
    Se <- sapply(c, function(c){sum(Y.cases > c)/n})
  }

  FPR <- c(0,1-Sp[order(1-Sp, Se)],1); TPR <- c(0,Se[order(1-Sp, Se)],1); NS <- length(FPR)
  # auc <- abs(sum((Se[order(Y)][-1] + Se[order(Y)][-NS])/2*(Sp[order(Y)][-NS] - Sp[order(Y)][-1])))
  auc <- sum(TPR[-NS]*diff(FPR))

  results <- list(levels = levels.names, X = X, Y = Y, D = D, Sp = Sp, Se = Se, auc = auc, type = type)
  if(type == 'lrm'){
    results$formula <- formula.lrm
    results$model <- model$coefficients
  }

  attr(results, 'class') <- 'hroc'

  if(plot.roc){
    if(new.window) dev.new(width = 5.5, height = 6)
    plot(results)
  }

  return(results)

}
