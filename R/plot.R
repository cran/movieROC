
plot.groc <- function(x, xlim = c(0,1), ylim = c(0,1), lwd = 3, xlab = "False-Positive Rate", ylab = "True-Positive Rate", main = "ROC curve", cex.lab = 1.25, cex.main = 1.5, type = NULL, new = TRUE, ...){

  obj <- x
  type <- ifelse(is.null(type), ifelse(obj$param, 'l', 's'), type)
  if(new){
    plot(c(0,obj$t,1), c(0,obj$roc,1), type = type, xlim = xlim, ylim = ylim, lwd = lwd, xlab = xlab, ylab = ylab, main = main, cex.lab = cex.lab, cex.main = cex.main, ...)
    axis(1, at = seq(0,1,0.01), labels = F, tck = -0.01)
    axis(1, at = seq(0.1,0.9,0.1), labels = F, tck = -0.02)
    axis(2, at = seq(0,1,0.01), labels = F, tck = -0.01)
    axis(2, at = seq(0.1,0.9,0.1), labels = F, tck = -0.02)
    abline(0, 1, lty = 2, col = "gray")
  }else{
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(new = TRUE)
    plot(c(0,obj$t,1), c(0,obj$roc,1), type = type, xlim = xlim, ylim = ylim, lwd = lwd, xlab = "", ylab = "", axes = FALSE, main = "", cex.lab = cex.lab, cex.main = cex.main, ...)
    abline(0, 1, lty = 2, col = "gray")
  }

}

plot.multiroc <- function(x, ...){
  plot.groc(x, type = 's', ...)
}

plot.hroc <- function(x, type = 'S', xlim = c(0,1), ylim = c(0,1), lwd = 3, xlab = "False-Positive Rate", ylab = "True-Positive Rate", main = "ROC Curve", cex.lab = 1.25, cex.main = 1.5, new = TRUE, ...){

  obj <- x
  if(new){
    plot(c(1,1-obj$Sp[order(obj$Y)],0), c(1,obj$Se[order(obj$Y)],0), type = type, xlim = xlim, ylim = ylim, lwd = lwd, main = main, xlab = xlab, ylab = ylab, cex.lab = cex.lab, cex.main = cex.main, ...)
    axis(1, at = seq(0,1,0.01), labels = F, tck = -0.01)
    axis(1, at = seq(0.1,0.9,0.1), labels = F, tck = -0.02)
    axis(2, at = seq(0,1,0.01), labels = F, tck = -0.01)
    axis(2, at = seq(0.1,0.9,0.1), labels = F, tck = -0.02)
    abline(0, 1, lty = 2, col = "gray")
  }else{
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(new = TRUE)
    plot(c(1,1-obj$Sp[order(obj$Y)],0), c(1,obj$Se[order(obj$Y)],0), type = type, xlim = xlim, ylim = ylim, lwd = lwd, main = "", xlab = "", ylab = "", axes = FALSE, cex.lab = cex.lab, cex.main = cex.main, ...)
    abline(0, 1, lty = 2, col = "gray")
  }

}
