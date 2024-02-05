# plot_densityROC <- function(obj, ...) {
#   UseMethod("plot_densityROC")
# }
# Only for side = "right" or "left"

plot_densityROC <- function(obj, h = c(1,1), C = NULL, build.process = FALSE, completeROC = TRUE, legends = FALSE,
                            rel.tol = 1e-3, par.specify = FALSE, cex.lab = 1.5, cex.axis = 1.25, cex.main = 1.75, lwd = 2, col = c('#485C99','#8F3D52'), col.roc = 'blue', ...){

  x <- obj
  levels.names <- obj$levels
  controls <- obj$controls; cases <- obj$cases; side <- obj$side
  mm <- min(controls, cases); MM <- max(controls, cases)
  if(length(h) == 1) h <- h*c(1,1)

  if(side != "right" & side != "left") stop("This function only works for standard ROC curve (right-sided or left-sided).")

  plot_density <- function(controls, cases, h){
    par(mar = c(4.1, 4.6, 3.1, 1.1))
    density_controls <- density(controls, adjust = h[1])
    density_cases <- density(cases, adjust = h[2])
    plot(density_controls$x, density_controls$y,
         xlim = c(min(c(density_controls$x, density_cases$x)), max(c(density_controls$x, density_cases$x))),
         ylim = c(-0.05-max(density_cases$y), 0.05+max(density_controls$y)),
         'l', col = col[1], lwd = lwd, xlab = "Marker (x)", ylab = "f(x)", main = "Density functions",
         cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    lines(density_cases$x, -density_cases$y, col = col[2], lwd = lwd)
    abline(h = 0, col = 'gray', lwd = 1)
    if(legends) legend('topright', inset = 0.03, c("Controls", "Cases"), col = col, lty = 1)
  }

  density_controls <- density(controls, adjust = h[1])
  density_cases <- density(cases, adjust = h[2])

  f.controls <- approxfun(density_controls$x, density_controls$y, yleft = 0, yright = 0)
  f.cases <- approxfun(density_cases$x, density_cases$y, yleft = 0, yright = 0)

  plot_ROCdensity <- function(controls, cases, side, build.process, p.ROC.C = 0, Se.ROC.C = 0, h){
    density_controls <- density(controls, adjust = h[1])
    density_cases <- density(cases, adjust = h[2])

    par(mar = c(4.1, 4.6, 3.1, 1.1))
    mm <- min(c(density_controls$x,density_cases$x)); MM <-  max(c(density_controls$x,density_cases$x))
    points <- seq(mm, MM, length.out = 200)
    integrate.density.p.ROC <- sapply(points, function(C) integrate(f.controls, ifelse(side=='right',C,mm), ifelse(side=='right',Inf,C),
                                                                    rel.tol = rel.tol, subdivisions = 1000, stop.on.error = FALSE)$value)
    integrate.density.Se.ROC <- sapply(points, function(C) integrate(f.cases, ifelse(side=='right',C,mm), ifelse(side=='right',Inf,C),
                                                                     rel.tol = rel.tol, subdivisions = 1000, stop.on.error = FALSE)$value)
    if(build.process){
      plot(integrate.density.p.ROC, integrate.density.Se.ROC, 'l', lwd = lwd, col = ifelse(completeROC,'gray','white'),
           xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
      if(side == 'right'){
        lines(integrate.density.p.ROC[integrate.density.p.ROC>=p.ROC.C & integrate.density.Se.ROC>=Se.ROC.C],
              integrate.density.Se.ROC[integrate.density.p.ROC>=p.ROC.C & integrate.density.Se.ROC>=Se.ROC.C], lwd = 2*lwd)
      }else{
        lines(integrate.density.p.ROC[integrate.density.p.ROC<=p.ROC.C & integrate.density.Se.ROC<=Se.ROC.C],
              integrate.density.Se.ROC[integrate.density.p.ROC<=p.ROC.C & integrate.density.Se.ROC<=Se.ROC.C], lwd = 2*lwd)
      }
    }else{
      plot(integrate.density.p.ROC, integrate.density.Se.ROC, 'l', lwd = 2*lwd, xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve",
           cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    }
    abline(0, 1, col='gray')
    if(side == 'right'){
      xaxis.C <- integrate.density.p.ROC[integrate.density.p.ROC>=p.ROC.C & integrate.density.Se.ROC>=Se.ROC.C]
      yaxis.C <- points[integrate.density.p.ROC>=p.ROC.C & integrate.density.Se.ROC>=Se.ROC.C]
    }else{
      xaxis.C <- integrate.density.p.ROC[integrate.density.p.ROC<=p.ROC.C & integrate.density.Se.ROC<=Se.ROC.C]
      yaxis.C <- points[integrate.density.p.ROC<=p.ROC.C & integrate.density.Se.ROC<=Se.ROC.C]
    }
    list(xaxis.C = xaxis.C, yaxis.C = yaxis.C)
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if(!par.specify){par(mfrow = c(1, 2))}

  mm <- min(c(density_controls$x, density_cases$x))
  MM <- max(c(density_controls$x, density_cases$x))
  ord.round <- 2 - round(log10(MM-mm))

  if(!is.null(C) && C < mm){
    warning("C is lower than all marker values")
    C <- NULL
  }else if(!is.null(C) && C > MM){
    warning("C is larger than all marker values")
    C <- NULL
  }


  if(is.null(C)){
    plot_density(controls, cases, h = h)
    ROC <- plot_ROCdensity(controls, cases, side = side, h = h, build.process = FALSE)
  }else{
    m <- matrix(c(1,3,3,2,4,4), 3, 2)
    layout(m)
    par(mar = c(4, 4, 2, 1))

    plot_density(controls, cases, h = h)
    if(side == 'right'){
      polygon(c(C, density_controls$x[density_controls$x > C], tail(density_controls$x[density_controls$x > C],1)),
              c(0, density_controls$y[density_controls$x > C], tail(density_controls$y[density_controls$x > C],1)),
              col = col[1], border = NA)
      polygon(c(C, density_cases$x[density_cases$x > C], tail(density_cases$x[density_cases$x > C],1)),
              c(0, -density_cases$y[density_cases$x > C], tail(-density_cases$y[density_cases$x > C],1)),
              col = col[2], border = NA)
    }else{
      polygon(c(C, density_controls$x[density_controls$x < C], tail(density_controls$x[density_controls$x < C],1)),
              c(0, density_controls$y[density_controls$x < C], tail(density_controls$y[density_controls$x < C],1)),
              col = col[1], border = NA)
      polygon(c(C, density_cases$x[density_cases$x < C], tail(density_cases$x[density_cases$x < C],1)),
              c(0, -density_cases$y[density_cases$x < C], tail(-density_cases$y[density_cases$x < C],1)),
              col = col[2], border = NA)
    }

    if(legends) legend('bottomright', inset = 0.03, c("1 - Sp", "Se"), fill = col)
    abline(v = C, col = col.roc, lty = 4)
    text(C, -0.05-max(density_cases$y), adj = 1, round(C, ord.round), cex = 1, col = col.roc)

    p.ROC.C <- integrate(f.controls, ifelse(side=='right',C,mm), ifelse(side=='right',Inf,C), rel.tol = rel.tol, subdivisions = 1000, stop.on.error = FALSE)$value
    Se.ROC.C <- integrate(f.cases, ifelse(side=='right',C,mm), ifelse(side=='right',Inf,C), rel.tol = rel.tol, subdivisions = 1000, stop.on.error = FALSE)$value
    ROC.until.C <- plot_ROCdensity(controls, cases, side=side, build.process, p.ROC.C = p.ROC.C, Se.ROC.C = Se.ROC.C, h = h)
    points(p.ROC.C, Se.ROC.C, cex = 1.25, pch = 19, col = col.roc)
    text(p.ROC.C, Se.ROC.C-0.025, round(C, ord.round), adj = 0, cex = 1, col = col.roc)

    plot(rep(-0.75,length(controls)), controls, 'p', pch = 1, col = col[1], xlim = c(-2.5,2.5), ylim = c(mm,MM), xaxt = 'n', xlab = "", ylab = "Marker values",
         xaxs = 'i', cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    points(rep(0.75,length(cases)), cases, col = col[2])
    boxplot(controls, at = -1.5, col = col[1], add = TRUE, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    boxplot(cases, at = 1.5, col = col[2], add = TRUE, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    abline(h = C, col = col.roc, lty = 4)
    polygon(c(2.3,2.3,2.5,2.5), c(ifelse(side=='right',C, mm),ifelse(side=='right',MM,C),ifelse(side=='right',MM,C),ifelse(side=='right',C,mm)), col = 'gray', border = NA)
    if(legends) legend('topright', inset = 0.03, c("Controls", "Cases"), col = col, pch = c(1,1))

    plot(0.5, controls[1], lwd = 2*lwd, col = 'white', xlim = c(0,1), ylim = c(mm,MM),
         xlab = "1-Specificity", ylab = "Marker intervals", main = "Classification subsets", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    polygon(c(ROC.until.C$xaxis.C, rev(ROC.until.C$xaxis.C)), c(ROC.until.C$yaxis.C, rep(ifelse(side=='right',MM,mm),length(ROC.until.C$yaxis.C))), col = 'gray', border = NA)
    if(legends) legend('topright', inset = 0.03, c("Classified as Controls","Classified as Cases"), fill = c('white', 'gray'))
  }

}
