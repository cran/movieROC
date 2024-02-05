### Plot ROC curve and classification regions from min-max method ###

.plot_minmaxROC <- function(X, D, roc.minmax, coef = NULL, FPR = 0.15, col = "blue", type = 's', cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, alpha.points = 1, alpha.contour = 0.25, cex.point = 1.5, lwd.curve = 2, lty.curve = 1, col.curve = 'black',  lf = NULL, xlab = "X1", ylab = "X2"){

  x <- X[,1]; y <- X[,2]
  lx <- (max(x) - min(x))/20;   ly <- (max(y) - min(y))/20

  obj <- roc.minmax
  t <- obj$t; roc <- obj$roc; points <- obj$c; auc <- obj$auc; obj$D <- D

  if(is.null(FPR)) FPR <- roc.minmax$t[which.max(roc.minmax$roc - roc.minmax$t)]
  index.max <- which.max(t[t <= FPR])
  T <- t[index.max]; ROC <- roc[index.max]; C <- points[t==T]

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  new <- FALSE
  par(fig = c(0,0.49,0,1), mar = c(5.1,5.1,4.1,2.1), new = new)

  colTrans <- col
  colTrans <- rgb(red=col2rgb(colTrans)[1], green=col2rgb(colTrans)[2], blue=col2rgb(colTrans)[3], alpha = alpha.contour*255, maxColorValue=255)

  plot(x, y, 'p', col=ifelse(obj$D, adjustcolor('red',alpha.f = alpha.points), adjustcolor('green4', alpha.f = alpha.points)), pch=16, cex=cex,
       xlim=range(x,finite=TRUE)+lx*c(-1,1),  ylim=range(y,finite=TRUE)+ly*c(-1,1), xaxs="i", yaxs="i", xlab = xlab, ylab = ylab, main = "Classification subsets", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)

  if(!is.null(roc.minmax$coef)) coef <- roc.minmax$coef

  Z <- c(tcrossprod(cbind(apply(X,1,max),apply(X,1,min)), matrix(coef, nrow = 1)))
  if(is.null(lf)) lf <- max(min(sort(unique(Z))[-1]-sort(unique(Z))[-length(unique(Z))]), .Machine$double.eps)
  xx <- seq(min(x)-lx,max(x)+lx, length.out=100)
  yy <- seq(min(y)-ly,max(y)+ly, length.out=100)
  f <- outer(xx,yy,function(x,y) c(tcrossprod(cbind(apply(cbind(x,y),1,max),apply(cbind(x,y),1,min)), matrix(coef, nrow = 1))))


  abline(0, 1, col = 'gray', lwd = 2.3, lty = 2)
  .filled.contour(xx, yy, f, levels=c(C,C+lf,max(f)+2*lf), col=c(col,colTrans))

  par(fig = c(0.51,1,0,1), new = TRUE)


  plot(c(0,t,1), c(0,roc,1), type = type, lwd = lwd.curve, lty = lty.curve, col = col.curve, xlim=c(0,1), ylim=c(0,1), xaxt = ifelse(new,"n","s"), yaxt = ifelse(new,"n","s"), xlab="1-Specificity", ylab="Sensitivity", main="ROC curve", cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)

  points(t[index.max],roc[index.max], pch=16, col=col, cex = cex.point)

  abline(0,1, col='gray', lty = 2)
  axis(1, at=seq(0,1,0.01), labels=F, tck=-0.01)
  axis(1, at=seq(0.1,0.9,0.1), labels=F, tck=-0.02)
  axis(2, at=seq(0,1,0.01), labels=F, tck=-0.01)
  axis(2, at=seq(0.1,0.9,0.1), labels=F, tck=-0.02)


}
