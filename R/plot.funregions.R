plot.funregions <- function(x, ...) {
  UseMethod("plot.funregions")
}

plot.funregions.hroc <- function(x, FPR = 0.15, FPR2 = NULL, plot.subsets = TRUE, new.window = FALSE, main = NULL, ylim = NULL, ...){

  obj <- x
  X <- obj$X; Y <- obj$Y; Sp <- obj$Sp; Se <- obj$Se; type <- obj$type
  indexX <- order(X)
  rangeX <- max(X)-min(X)
  Xfun <- sort(unique(c(seq(min(X), max(X), length.out=1000), X)))

  C <- ifelse(type == 'overfitting', 1-FPR, Y[which.min(ifelse(1-Sp <= FPR, FPR-1+Sp, 1))])
  if(min(ifelse(1-Sp <= FPR, FPR-1+Sp, 1)) == 1){C <- ifelse(type=='overfitting', 1.1, max(Y))}
  h <- approx(X[indexX], Y[indexX], xout = Xfun)$y
  Xcol <- Xfun[h >= C]

  if(is.null(ylim)) ylim <- c(min(Y)-(max(Y)-min(Y))/20, max(Y))

  if(new.window) dev.new(width = 6, height = 5)

  plot(X[indexX], Y[indexX], type = 'l', xlab = '', ylab = '', yaxt = 'n', frame = FALSE,
       main = ifelse(is.null(main), paste("Model:", type, ifelse(type=='lrm', obj$formula, "")), main),
       ylim = ylim, lwd = 2)
  if(plot.subsets) axis(1, at = Xcol, tcl = 0.6, labels = F, col = 'blue')

  subsets <- function(Xcol){
    streaks <- rle(is.element(Xfun, Xcol))
    fin <- Xfun[cumsum(streaks$lengths)[streaks$values]]
    inicio <- Xfun[cumsum(streaks$lengths)[!streaks$values]+1]
    if(streaks$values[1]){
      inicio <- c(Xfun[1], inicio)
    }
    if(!streaks$values[length(streaks$values)]){
      inicio <- inicio[-length(inicio)]
    }
    cbind(inicio, fin)
  }
  ClassSubsets <- subsets(Xcol)

  if(!is.null(FPR2)){
    C2 <- ifelse(type == 'overfitting', 1-FPR2, Y[which.min(ifelse(1-Sp <= FPR2, FPR2-1+Sp, 1))])
    if(min(ifelse(1-Sp <= FPR2, FPR2-1+Sp, 1)) == 1){C2 <- ifelse(type == 'overfitting', 1.1, max(Y))}
    Xcol <- Xfun[h >= C2]
    if(plot.subsets){
      axis(1, at = Xcol, tcl = 0.3, labels = F, col = 'red')
      abline(h = C2, lty = 2, col = 'lightcoral', lwd = 2)
      text(min(X) + .1*rangeX, C2, paste("FPR=", FPR2, sep = ''), pos = ifelse(FPR2<=0.5, 1, 3), col = 'lightcoral')
    }
    ClassSubsets2 <- subsets(Xcol)
  }

  order.space <- round(log10(rangeX))
  space <- 10^order.space

  if(plot.subsets){
    axis(1, at = Xfun, tcl = 0, labels = F)
    abline(h = C, lty = 2, col = 'lightblue', lwd = 2)
    text(min(X) + .1*rangeX, C, paste("FPR=", FPR, sep = ''), pos = ifelse(FPR <= 0.5, 1, 3), col = 'lightblue')
    axis(side = 1, at = seq(round(min(X),-order.space), round(max(X),-order.space), space/20), tcl = -0.2, labels = FALSE)
  }else{
    axis(side = 1, at = seq(round(min(X),-order.space), round(max(X),-order.space), space/20), tcl = -0.2, labels = FALSE)
  }

  if(is.null(FPR2)){
    return(list(ClassSubsets = ClassSubsets))
  }else{
    return(list(ClassSubsets = ClassSubsets, ClassSubsets2 = ClassSubsets2))
  }

}
