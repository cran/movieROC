plot_densities <- function(obj, h = c(1,1), histogram = FALSE, breaks = 15, col = c('#485C99','#8F3D52'), lwd = 2, xlim = NULL, ylim = NULL, xaxs = "i", yaxs = "i", 
                           xlab = "Marker", ylab = "f(x)", main = "Density functions", legend = FALSE, pos.legend = "topright",
                           cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, cex.legend = 1,
                           eps = sqrt(.Machine$double.eps), new = TRUE, ...){

  x <- obj
  controls <- obj$controls; cases <- obj$cases; side <- obj$side
  mm <- min(controls, cases); MM <- max(controls, cases)
  if(length(h) == 1) h <- h*c(1,1)

  if(side=='right' | side=='left') c <- pmax(pmin(x$c, MM), mm)
  if(side=='both' | side=='both2') xl <- pmax(x$xl, mm); xu <- pmin(x$xu, MM)
  
  if(histogram){
    hist.controls <- hist(controls, breaks=breaks, plot=FALSE)
    hist.cases <- hist(cases, breaks=breaks, plot=FALSE)
    x.dcontrols <- sort(c(hist.controls$breaks, hist.controls$breaks + eps)); y.dcontrols <- c(0, rep(hist.controls$density, each = 2), 0)
    x.dcases <- sort(c(hist.cases$breaks, hist.cases$breaks + eps)); y.dcases <- c(0, rep(hist.cases$density, each=2), 0)
  }else{
    x.dcontrols <- density(controls, adjust=h[1])$x; y.dcontrols <- density(controls, adjust=h[1])$y
    x.dcases <- density(cases, adjust=h[2])$x; y.dcases <- density(cases, adjust=h[2])$y
  }

  mm <- min(c(mm, x.dcontrols[y.dcontrols>0.01], x.dcases[y.dcases>0.01]))
  MM <- max(c(MM, x.dcontrols[y.dcontrols>0.01], x.dcases[y.dcases>0.01]))

  My <- max(c(y.dcases,y.dcontrols))
  if(is.null(ylim)) ylim <- c(-1.2*My, 1.2*My)
  if(is.null(xlim)) xlim <- c(mm, MM)

  if(new){
    plot(x.dcontrols, y.dcontrols, xlim = xlim, ylim = ylim, 'l', col = col[1], lwd = lwd, 
         xlab = xlab, ylab = ylab, main = main, xaxs = xaxs, yaxs = yaxs, frame = FALSE,
         cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main, yaxt = "n")
  }else{
    lines(x.dcontrols, y.dcontrols, col = col[1], lwd = lwd)
  }

  lines(x.dcases, -y.dcases, col = col[2], lwd = lwd)
  abline(h = 0, col = 'gray', lwd = lwd/2, lty = 2)

  ticks.axis1 <- axis(1, cex.axis = cex.axis)
  space0 <- max(diff(ticks.axis1))
  axis(1, at = seq(min(ticks.axis1)-space0, max(ticks.axis1)+space0, space0/4), tcl = -0.4, labels = FALSE)
  axis(1, at = seq(min(ticks.axis1)-space0, max(ticks.axis1)+space0, space0/8), tcl = -0.3, labels = FALSE)

  ticks.axis2 <- as.numeric(axis(2, cex.axis = cex.axis, labels = FALSE))
  axis(2, cex.axis = cex.axis, at = ticks.axis2, labels = abs(ticks.axis2))
  space0 <- max(diff(ticks.axis2))
  axis(2, at = seq(min(ticks.axis2)-space0, max(ticks.axis2)+space0, space0/2), tcl = -0.3, labels = FALSE)

  if(legend) legend(pos.legend, c("Controls", "Cases"), col = col, lwd = lwd, inset = 0.03, bty = "n", cex = cex.legend)

}
