plot_buildROC <- function(x, ...) {
  UseMethod("plot_buildROC")
}


plot_buildROC.groc <- function(x, FPR = NULL, C = NULL, XL = NULL, XU = NULL, h = c(1,1), histogram = FALSE, breaks = 15,
                               reduce = TRUE, build.process = FALSE, completeROC = FALSE,
                               new.window = FALSE, legends = FALSE, type = 's',
                               cex.point = 1.5, lwd.curve = 2, mar = NULL, lim.density = 0.01, xlim = NULL, ylim = NULL,
                               cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, xlab = "Marker", main.density = "Density functions",
                               col.controlscases = c('#485C99','#8F3D52'), col.threshold = '#FCBA04', col.curve = 'black',
                               eps = sqrt(.Machine$double.eps), ...){

  controls <- x$controls; cases <- x$cases; levels.names <- x$levels; side <- x$side
  t <- x$t; roc <- x$roc
  if(x$param) type <- 'l'
  if(is.null(FPR) & is.null(C) & (is.null(XL) | is.null(XU))) FPR <- .15

  mm <- min(controls, cases); MM <- max(controls,cases)
  ord.round <- 2 - round(log10(MM-mm))

  if(x$param){
    mcontrols <- mean(controls); scontrols <- sd(controls)
    mcases <- mean(cases); scases <- sd(cases)
    t[t==0] <- .Machine$double.eps; t[t==1] <- 1-.Machine$double.eps
    controls <- qnorm(t, mcontrols, scontrols)
    cases <- qnorm(t, mcases, scases)
  }
  if(side == 'right' | side == 'left') c <- pmax(pmin(x$c, MM), mm)
  if(side == 'both' | side == 'both2') xl <- pmax(x$xl, mm); xu <- pmin(x$xu, MM)

  if(histogram){
    hist.controls <- hist(controls, breaks = breaks, plot = FALSE)
    hist.cases <- hist(cases, breaks = breaks, plot = FALSE)
    x.dcontrols <- sort(c(mm, hist.controls$breaks, hist.controls$breaks - eps, MM)); y.dcontrols <- c(0, 0, c(rbind(hist.controls$density, hist.controls$density)), 0, 0)
    x.dcases <- sort(c(mm, hist.cases$breaks, hist.cases$breaks - eps, MM)); y.dcases <- c(0, 0, c(rbind(hist.cases$density, hist.cases$density)), 0, 0)
  }else{
    x.dcontrols <- density(controls, adjust = h[1])$x; y.dcontrols <- density(controls, adjust = h[1])$y
    x.dcases <- density(cases, adjust = h[2])$x; y.dcases <- density(cases, adjust = h[2])$y
  }
  if(x$param){
    x.dcontrols <- c(mm, controls, MM); y.dcontrols <- c(0, dnorm(controls, mcontrols, scontrols), 0)
    x.dcases <- c(mm, cases, MM); y.dcases <- c(0, dnorm(cases, mcases, scases), 0)
  }
  t[1] <- 0; t[length(t)] <- 1

  if(!is.null(FPR) && FPR < 0){
    warning("FPR is lower than 0.")
    FPR <- NULL
  }else if(!is.null(FPR) && FPR > 1){
    warning("FPR is larger than 1.")
    FPR <- NULL
  }

  mm <- min(c(mm, x.dcontrols[y.dcontrols > lim.density], x.dcases[y.dcases > lim.density]))
  MM <- max(c(MM, x.dcontrols[y.dcontrols > lim.density], x.dcases[y.dcases > lim.density]))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if(!is.null(FPR) | (side %in% c("right","left") & !is.null(C)) | (side %in% c("both","both2") & !is.null(XL) & !is.null(XU))){

    if(reduce){

      if(new.window) dev.new(width = 6, height = 3.5, unit = "in")
      m <- matrix(c(1,1,2), 1, 3)
      layout(m)
      if(is.null(mar)) mar <- c(5.1, 5.1, 4.1, 2.1)
      par(mar = mar)

    }else{

      if(new.window) dev.new(width = 6, height = 8, unit = "in")
      m <- matrix(c(1,3,3,2,4,4), 3, 2)
      layout(m)
      if(is.null(mar)) mar <- c(6.1, 6.1, 3.1, 1.1)
      par(mar = mar)

    }

    My <- max(c(y.dcases,y.dcontrols))
    if(is.null(ylim)) ylim <- c(-1.2*My, 1.2*My)
    if(is.null(xlim)) xlim <- c(mm, MM)
    plot(x.dcontrols, y.dcontrols, xlim = xlim, ylim = ylim, 'l', col = col.controlscases[1], lwd = lwd.curve,
         xlab = xlab, ylab = "f(x)", main = main.density, yaxt = "n", xaxs = 'i', cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    lines(x.dcases, -y.dcases, col = col.controlscases[2], lwd = lwd.curve)
    abline(h = 0, col = 'gray', lwd = lwd.curve)

    if(legends) legend('topright', inset = 0.03, c("Controls", "Cases"), col = col.controlscases, lty = 1)

    potax <- trunc(log10(diff(ylim)/2)) - (log10(diff(ylim)/2) < 0)
    axisy <- seq(round(ylim[1],-potax), round(ylim[2],-potax), (10^potax)/2)
    axisy[abs(axisy) < .Machine$double.eps] <- 0
    axis(2, at = axisy[seq(1,length(axisy),2)], abs(axisy[seq(1,length(axisy),2)]), tck = -0.03, cex.axis = cex.axis)
    axis(2, at = axisy, labels = F, tck = -0.02, cex.axis = cex.axis)

    output.predict <- predict(x, FPR = FPR, C = C, XL = XL, XU = XU)
    T <- 1 - output.predict$Specificity; ROC <- output.predict$Sensitivity


    if(side == 'right' | side == 'left'){
      if(side=='right'){
        C <- output.predict$ClassSubsets[1]
        polygon(c(C, C, x.dcontrols[x.dcontrols >= C], max(x.dcontrols), max(x.dcontrols) + eps),
                c(0, y.dcontrols[x.dcontrols >= C][1], y.dcontrols[x.dcontrols >= C], tail(y.dcontrols[x.dcontrols >= C],1), 0),
                col = col.controlscases[1], border = NA)
        polygon(c(C, C, x.dcases[x.dcases >= C], max(x.dcases), max(x.dcases) + eps),
                c(0, -y.dcases[x.dcases >= C][1], -y.dcases[x.dcases >= C], tail(-y.dcases[x.dcases >= C],1), 0),
                col = col.controlscases[2], border = NA)
      }else{
        C <- output.predict$ClassSubsets[2]
        polygon(c(min(x.dcontrols) - eps, min(x.dcontrols), x.dcontrols[x.dcontrols <= C], C, C),
                c(0, y.dcontrols[x.dcontrols <= C][1], y.dcontrols[x.dcontrols <= C], tail(y.dcontrols[x.dcontrols <= C],1), 0),
                col = col.controlscases[1], border = NA)
        polygon(c(min(x.dcases) - eps, min(x.dcases), x.dcases[x.dcases <= C], C, C),
                c(0, -y.dcases[x.dcases <= C][1], -y.dcases[x.dcases <= C], tail(-y.dcases[x.dcases <= C],1), 0),
                col = col.controlscases[2], border = NA)
      }
      if(legends) legend('bottomright', inset = 0.03, c("1 - Sp", "Se"), fill = col.controlscases)
      abline(v = C, col = col.threshold, lty = 4)
      text(C, -0.05-max(y.dcases), adj = 1, round(C, ord.round), cex = cex.point, col = 'gray30')
    }else{
      if(side == 'both'){
        XL <- output.predict$ClassSubsets[1,2]; XU <- output.predict$ClassSubsets[2,1]
        polygon(c(mm, x.dcontrols[x.dcontrols <= XL], XL, XL), c(0, y.dcontrols[x.dcontrols <= XL], max(tail(y.dcontrols[x.dcontrols <= XL],1),0), 0), col = col.controlscases[1], border = NA)
        polygon(c(mm, x.dcases[x.dcases <= XL], XL, XL), c(0, -y.dcases[x.dcases <= XL], min(tail(-y.dcases[x.dcases <= XL],1),0), 0), col = col.controlscases[2], border = NA)
        polygon(c(XU, XU, x.dcontrols[x.dcontrols >= XU], MM), c(0, max(head(y.dcontrols[x.dcontrols >= XU],1),0), y.dcontrols[x.dcontrols >= XU], 0), col = col.controlscases[1], border = NA)
        polygon(c(XU, XU, x.dcases[x.dcases >= XU], MM), c(0, min(head(-y.dcases[x.dcases >= XU],1), 0), -y.dcases[x.dcases >= XU], 0), col = col.controlscases[2], border = NA)
      }else{
        XL <- output.predict$ClassSubsets[1]; XU <- output.predict$ClassSubsets[2]
        polygon(c(XL, XL, x.dcontrols[x.dcontrols <= XU & x.dcontrols >= XL], XU, XU),
                c(0, tail(y.dcontrols[x.dcontrols <= XL],1), y.dcontrols[x.dcontrols <= XU & x.dcontrols >= XL], head(y.dcontrols[x.dcontrols >= XU],1), 0),
                col = col.controlscases[1], border = NA)
        polygon(c(XL, XL, x.dcases[x.dcases <= XU & x.dcases >= XL], XU, XU),
                c(0, tail(-y.dcases[x.dcases <= XL],1), -y.dcases[x.dcases <= XU & x.dcases >= XL], head(-y.dcases[x.dcases >= XU],1), 0),
                col = col.controlscases[2], border = NA)
      }
      if(legends) legend('bottomright', inset = 0.03, c("1 - Sp", "Se"), fill = col.controlscases)
      abline(v = XL, col = col.threshold, lty = 4); abline(v = XU, col = col.threshold, lty = 4)
      text(XL, -0.05-max(y.dcases), adj = 1, round(XL, ord.round), cex = cex.point, col = 'gray30')
      text(XU, -0.05-max(y.dcases), adj = 0, round(XU, ord.round), cex = cex.point, col = 'gray30')
    }

    if(build.process){
      plot(c(0,t,1), c(0,roc,1), type = type, lwd = lwd.curve/2, col = ifelse(completeROC,'gray','white'), xlim = c(0,1), ylim = c(0,1),
           xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
      lines(c(0,t[t<=T & roc<=ROC],T), c(0,roc[t<=T & roc<=ROC],ROC), col = col.curve, lwd = lwd.curve, type = type)
      abline(0, 1, col = 'gray', lty = 2)
    }else{
      plot(c(0,t,1), c(0,roc,1), type = type, col = col.curve, lwd = lwd.curve, xlim = c(0,1), ylim = c(0,1),
           xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    }

    points(T, ROC, cex = cex.point, pch = 19, col = col.threshold)
    if(side == 'right' | side == 'left'){
      text(T, ROC-0.05, round(C, ord.round), adj = 0.5, cex = cex.point, col = 'gray30')
    }else{
      text(T, ROC-0.05, paste0("(", round(XL, ord.round), ",", round(XU, ord.round),")"), adj = 0.5, cex = cex.point, col = 'gray30')
    }

    if(!reduce){
      plot(rep(-0.75,length(controls)), controls, 'p', pch = 1, col = col.controlscases[1], xlim = c(-2.5,2.5), ylim = c(mm,MM),
           xaxt = 'n', xlab = "", ylab = "Marker values", xaxs = 'i', yaxs = 'i', cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
      points(rep(0.75,length(cases)), cases, col = col.controlscases[2])
      boxplot(controls, at = -1.5, col = col.controlscases[1], add = TRUE, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
      boxplot(cases, at = 1.5, col = col.controlscases[2], add = TRUE, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
      if(side == 'right' | side == 'left'){
        abline(h = C, col = col.threshold, lty = 4)
        polygon(c(2.3,2.3,2.5,2.5), c(ifelse(side=='left',mm,C),ifelse(side=='left',C,MM),ifelse(side=='left',C,MM),ifelse(side=='left',mm,C)), col = 'gray30', border = NA)
      }else{
        abline(h = XL, col = col.threshold, lty = 4); abline(h = XU, col = col.threshold, lty = 4)
        if(side=='both'){
          polygon(c(2.3,2.3,2.5,2.5), c(mm, XL, XL, mm), col = 'gray', border = NA)
          polygon(c(2.3,2.3,2.5,2.5), c(XU, MM, MM, XU), col = 'gray', border = NA)
        }else{
          polygon(c(2.3,2.3,2.5,2.5), c(XL, XU, XU, XL), col = 'gray', border = NA)
        }
      }
      if(legends) legend('topright', inset = 0.03, c("Controls", "Cases"), col = col.controlscases, pch = c(1,1))

      colTrans <- 'gray'
      colTrans <- rgb(red = col2rgb(colTrans)[1], green = col2rgb(colTrans)[2], blue = col2rgb(colTrans)[3], alpha = 0.25*255, maxColorValue = 255)

      plot(0.5, controls[1], lwd = lwd.curve, col = 'white', xlim = c(0,1), ylim = c(mm,MM),
           xlab = "1-Specificity", ylab = "Marker intervals", main = "Classification subsets", yaxs = 'i', cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
      if(!build.process) T <- 1
      if(side == 'right' | side == 'left'){
        polygon(c(t[t<=T], T,T,rev(t[t<=T])), c(c[t<=T],C,C, rep(ifelse(side=='left',mm,MM),length(c[t<=T]))), col = colTrans, border = NA)
        segments(c(t[t<=T],T), c(c[t<=T],C), c(t[t<=T],T), rep(ifelse(side=='left',mm,MM),length(c[t<=T])+1), col = 'gray')
      }else{
        if(side == 'both'){
          polygon(c(t[t<=T],T,T, rev(t[t<=T])), c(xu[t<=T],XU,XU, rep(MM,length(xu[t<=T]))), col = colTrans, border = NA)
          polygon(c(t[t<=T],T,T, rev(t[t<=T])), c(xl[t<=T],XL,XL, rep(mm,length(xl[t<=T]))), col = colTrans, border = NA)
          segments(c(t[t<=T],T), c(xu[t<=T],XU), c(t[t<=T],T), rep(MM,length(xu[t<=T])+1), col = 'gray')
          segments(c(t[t<=T],T), c(xl[t<=T],XL), c(t[t<=T],T), rep(mm,length(xl[t<=T])+1), col = 'gray')
        }else{
          polygon(c(t[t<=T], T,T, rev(t[t<=T])), c(xl[t<=T], XL, XU, rev(xu[t<=T])), col = colTrans, border = NA)
          segments(c(t[t<=T],T), c(xl[t<=T],XL), c(t[t<=T],T), c(xu[t<=T],XU), col = 'gray')
        }
      }
      if(legends) legend('topright', inset = 0.03, c("Classified as Controls","Classified as Cases"), fill = c('white', 'gray'))
    }




  }else{

    if(new.window) dev.new(width = 9, height = 5, unit = "in")
    par(mfrow = c(1, 2))
    if(is.null(mar)) mar <- c(5.1, 5.1, 4.1, 2.1)
    par(mar = mar)

    My <- max(c(y.dcases,y.dcontrols))
    if(is.null(ylim)) ylim <- c(-1.2*My, 1.2*My)
    if(is.null(xlim)) xlim <- c(mm, MM)
    plot(x.dcontrols, y.dcontrols, xlim = xlim, ylim = ylim, 'l', col = col.controlscases[1], lwd = lwd.curve,
         xlab = xlab, ylab = "f(x)", main = main.density, yaxt = "n", xaxs = 'i', cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    lines(x.dcases, -y.dcases, col = col.controlscases[2], lwd = lwd.curve)
    abline(h = 0, col = 'gray', lwd = lwd.curve/2)

    potax <- trunc(log10(diff(ylim)/2)) - (log10(diff(ylim)/2) < 0)
    axisy <- seq(round(ylim[1],-potax), round(ylim[2],-potax), (10^potax)/2)
    axisy[abs(axisy) < .Machine$double.eps] <- 0
    axis(2, at = axisy[seq(1,length(axisy),2)], abs(axisy[seq(1,length(axisy),2)]), tck = -0.03, cex.axis = cex.axis)
    axis(2, at = axisy, labels = F, tck = -0.02, cex.axis = cex.axis)

    if(legends) legend('topright', inset = 0.03, c("Controls", "Cases"), col = col.controlscases, lty = 1)

    plot(c(0,t,1), c(0,roc,1), type = type, xlim = c(0,1), ylim = c(0,1), col = col.curve, lwd = lwd.curve,
         xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    abline(0, 1, col = 'gray', lty = 2)

  }

}


plot_buildROC.multiroc <- function(x, FPR = 0.15, display.method = c("PCA", "OV"), displayOV = c(1,2), build.process = FALSE, completeROC = TRUE,
                                   new = FALSE, new.window = FALSE, border = FALSE, cutoff = TRUE, legends = FALSE, type = 's',
                                   col.controlscases = c('#485C99','#8F3D52'), col.threshold = '#FCBA04', col.curve = 'black',
                                   cex.point = 1.5, alpha.points = .75, alpha.contour = 0.25, lwd.curve = 2, lty.curve = 1,
                                   cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, xlab = NULL, ylab = NULL,
                                   lf = NULL, eps = sqrt(.Machine$double.eps), ...){
  
  obj <- x
  X <- obj$X; D <- obj$D; Z <- obj$Z
  t <- obj$t; roc <- obj$roc; points <- obj$c; auc <- obj$auc
  
  if(!is.null(FPR) && FPR < 0){
    warning("FPR is lower than 0.")
    FPR <- NULL
  }else if(!is.null(FPR) && FPR > 1){
    warning("FPR is larger than 1.")
    FPR <- NULL
  }
  
  ## Create the 2D grid and calculate the result of the combination over every pair of points
  
  p <- ncol(X)
  
  display.method <- match.arg(display.method)
  
  if(display.method == "OV"){
    
    if(is.null(xlab)) xlab <- paste0("X.", displayOV[1])
    if(is.null(ylab)) ylab <- paste0("X.", displayOV[2])
    if(is.integer(displayOV) | length(displayOV)!=2) stop("displayOV should be a vector of length two indicating the variables (columns in X) used for displaying the results.")
    if(min(displayOV) < 1 | max(displayOV) > p) stop(paste0("displayOV should be a vector with two integers between 1 and dimension p = ", p, "."))
    
    x <- X[,displayOV[1]]; y <- X[,displayOV[2]]
    lx <- (max(x)-min(x))/20; ly <- (max(y)-min(y))/20
    xx <- seq(min(x)-lx, max(x)+lx, length.out = 100)
    yy <- seq(min(y)-ly, max(y)+ly, length.out = 100)
    
    if(obj$method == c("lrm")){
      
      if(is.null(lf)) lf <- max(min(sort(unique(Z))[-1]-sort(unique(Z))[-length(unique(Z))]), eps) else lf <- lf
      
      if(p == 2){
        f <- outer(xx,yy,function(x,y) suppressWarnings(predict(obj$model, data.frame(X = matrix(cbind(x,y), length(x), 2)))))
      }else{
        f <- outer(xx, yy, function(x,y){
          out <- matrix(0, length(x), p)
          out[,displayOV[1]] <- x; out[,displayOV[2]] <- y
          if(obj$method == "lrm"){
            suppressWarnings(Vectorize(predict(obj$model, data.frame(X = out))))
          }else{
            c(tcrossprod(out, matrix(obj$coefLinear, nrow = 1)))
          }
        })
      }
      
    }else if(obj$method == "fixedLinear"){
      
      if(is.null(lf)) lf <- max(min(sort(unique(Z))[-1]-sort(unique(Z))[-length(unique(Z))]), eps) else lf <- lf
      
      if(p == 2){
        f <- outer(xx, yy, function(x,y) c(tcrossprod(cbind(x,y), matrix(obj$coefLinear, nrow = 1))))
      }else{
        f <- outer(xx, yy, function(x,y){
          out <- matrix(0, length(x), p)
          out[,displayOV[1]] <- x; out[,displayOV[2]] <- y
          if(obj$method == "lrm"){
            suppressWarnings(Vectorize(predict(obj$model, data.frame(X = out))))
          }else{
            c(tcrossprod(out, matrix(obj$coefLinear, nrow = 1)))
          }
        })
      }
      
    }else if(obj$method == "fixedQuadratic"){
      
      if(is.null(lf)) lf <- max(min(sort(unique(Z))[-1]-sort(unique(Z))[-length(unique(Z))]), eps) else lf <- lf
      
      f <- outer(xx, yy, function(x,y) c(tcrossprod(cbind(x,y,x*y,x^2,y^2), matrix(obj$coefQuadratic, nrow = 1))))
      
    }else if(p == 2 & obj$method %in% c("dynamicMeisner", "dynamicEmpirical")){ # p=2 AND obj$method == "dynamicMeisner" OR obj$method == "dynamicEmpirical"
      
      CoefTable <- obj$CoefTable
      m <- sum(D == levels(as.factor(D))[1])
      lf <- ifelse(is.null(lf), mean(sapply(CoefTable, function(x) x$lf)), lf)
      f <- lapply(1:m, function(i){CoefTable[[i]]$f})
      
    }else if(p > 2 & obj$method == "dynamicMeisner"){ # p>2 AND obj$method == "dynamicMeisner"
      
      CoefTable <- obj$CoefTable
      m <- nrow(obj$controls)
      lf <- ifelse(is.null(lf), mean(sapply(CoefTable, function(x) x$lf)), lf)
      f <- lapply(1:m, function(i){
        outer(xx, yy, function(x,y){
          out <- matrix(0, length(x), p)
          out[,displayOV[1]] <- x; out[,displayOV[2]] <- y
          c(tcrossprod(out, matrix(obj$coefLinear[,i], nrow = 1)))
        })
      })
      
    }else if(p == 2 & obj$method == "kernelOptimal"){
      
      if (is.null(lf)) 
        lf <- max(min(sort(unique(Z))[-1] - sort(unique(Z))[-length(unique(Z))]), eps)
      else lf <- lf
      optimalT <- obj$optimalT
      f <- outer(xx, yy, function(x, y) optimalT(cbind(x, y)))
      warning(paste0(round(sum(is.na(f))/length(f) * 100), 
                     "% of values in the grid has a NaN optimalT. These have been replaced by adjacent value."))
      f <- zoo::na.locf(f)
      f <- zoo::na.locf(f, fromLast = TRUE)
      
    }
    
  }else{ # if display.method = "PCA"
    
    if(is.null(xlab)) xlab <- "PC1"
    if(is.null(ylab)) ylab <- "PC2"
    
    PC.output <- prcomp(data.frame(X), scale = FALSE)
    x <- PC.output$x[,1]; y <- PC.output$x[,2]
    lx <- (max(x)-min(x))/20; ly <- (max(y)-min(y))/20
    xx <- seq(min(x)-lx, max(x)+lx, length.out = 100)
    yy <- seq(min(y)-ly, max(y)+ly, length.out = 100)
    
    x.PC <- function(PC.coor){
      diag(p) %*% solve(t(PC.output$rotation)) %*% t(PC.coor)
    }
    
    if(obj$method %in% c("lrm", "fixedLinear")){
      
      f <- outer(xx, yy, function(x,y){
        out.pc <- matrix(0, length(x), p)
        out.pc[,1] <- x; out.pc[,2] <- y
        out <- matrix(x.PC(out.pc), ncol = p, byrow = TRUE)
        if(obj$method == "lrm"){
          suppressWarnings(Vectorize(predict(obj$model, data.frame(X = out))))
        }else{
          c(tcrossprod(out, matrix(obj$coefLinear, nrow = 1)))
        }
      })
      
    }else if(obj$method == "dynamicMeisner"){ # obj$method == "dynamicMeisner"
      
      CoefTable <- obj$CoefTable
      m <- nrow(obj$controls)
      f <- lapply(1:m, function(i){
        outer(xx, yy, function(x,y){
          out.pc <- matrix(0, length(x), p)
          out.pc[,1] <- x; out.pc[,2] <- y
          out <- matrix(x.PC(out.pc), ncol = p, byrow = TRUE)
          c(tcrossprod(out, matrix(obj$coefLinear[,i], nrow = 1)))
        })
      })
      
    }
    
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if(is.null(FPR)){
    
    if(new.window) dev.new(width = 9, height = 5, unit = "in")
    par(fig = c(0,0.49,0,1), mar = c(5.1,5.1,4.1,2.1), new = new)
    
    plot(x, y, 'p', col = ifelse(D, adjustcolor(col.controlscases[2], alpha.f = alpha.points), adjustcolor(col.controlscases[1], alpha.f = alpha.points)),
         pch = 16, cex = cex, xlim = range(x,finite=TRUE)+lx*c(-1,1), ylim = range(y,finite=TRUE)+ly*c(-1,1), xaxs = "i", yaxs = "i", xaxt = ifelse(new,"n","s"), yaxt = ifelse(new,"n","s"),
         xlab = xlab, ylab = ylab, main = "Classification subsets", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    
    par(fig = c(0.51,1,0,1), new = TRUE)
    plot(c(0,t,1), c(0,roc,1), type = type, xlim = c(0,1), ylim = c(0,1), lwd = lwd.curve, lty = lty.curve, col = col.curve,
         cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main, xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve")
    abline(0, 1, col='gray', lty = 2)
    axis(1, at = seq(0,1,0.01), labels = F, tck = -0.01)
    axis(1, at = seq(0.1,0.9,0.1), labels = F, tck = -0.02)
    axis(2, at = seq(0,1,0.01), labels = F, tck = -0.01)
    axis(2, at = seq(0.1,0.9,0.1), labels = F, tck = -0.02)
    
  }else{
    
    index.max <- which.max(t[t <= FPR])
    T <- t[index.max]; ROC <- roc[index.max]; C <- points[t==T]
    
    if(new.window) dev.new(width = 9, height = 5, unit = "in")
    par(fig = c(0,0.49,0,1), mar = c(5.1,5.1,4.1,2.1), new = new)
    
    colTrans <- col.threshold
    colTrans <- rgb(red = col2rgb(colTrans)[1], green = col2rgb(colTrans)[2], blue = col2rgb(colTrans)[3], alpha = alpha.contour*255, maxColorValue = 255)
    
    if(p == 2){
      
      plot(x, y, 'p', col = ifelse(D, adjustcolor(col.controlscases[2],alpha.f = alpha.points), adjustcolor(col.controlscases[1], alpha.f = alpha.points)),
           pch = 16, cex = cex, xlim = range(x,finite=TRUE)+lx*c(-1,1), ylim = range(y,finite=TRUE)+ly*c(-1,1), xaxs = "i", yaxs = "i", xaxt = ifelse(new,"n","s"), yaxt = ifelse(new,"n","s"),
           xlab = xlab, ylab = ylab, main = "Classification subsets", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
      
      if(obj$method %in% c("dynamicMeisner","dynamicEmpirical")){
        .filled.contour(xx, yy, f[[index.max]], levels = c(C,max(f[[index.max]])), col = colTrans)
        if(border){
          abline(C/obj$coefLinear[2,index.max], -obj$coefLinear[1,index.max]/obj$coefLinear[2,index.max], col = col.threshold, lwd = lwd.curve)
        }
      }else{
        .filled.contour(xx, yy, f, levels = c(C,max(f)), col = colTrans)
        if(border){
          if(obj$method != "fixedLinear") contour(xx, yy, f, nlevels = 1, levels = C, col = col.threshold, add = TRUE, drawlabels = FALSE, lwd = lwd.curve)
          if(obj$method == "fixedLinear") abline(C/obj$coefLinear[2], -obj$coefLinear[1]/obj$coefLinear[2], col = col.threshold, lwd = lwd.curve)
        }
      }
      
    }
    
    if(p > 2){
      
      if(obj$method == "dynamicMeisner"){
        
        plot(x, y, 'p', col = ifelse(D, adjustcolor(col.controlscases[2],alpha.f = alpha.points), adjustcolor(col.controlscases[1], alpha.f = alpha.points)),
             pch = ifelse(Z[,index.max] > C, 16, 15), cex = cex, xaxs = "i", yaxs = "i", xaxt = ifelse(new,"n","s"), yaxt = ifelse(new,"n","s"),
             xlim = range(x,finite=TRUE)+lx*c(-1,1), ylim = range(y,finite=TRUE)+ly*c(-1,1),
             xlab = xlab, ylab = ylab, main = "Classification subsets", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
        
        points(x, y, col = ifelse(Z[,index.max] > C, col.threshold, 'gray'), cex = 1.3*cex, pch = ifelse(Z[,index.max] > C, 21, 22))
        if(border) contour(xx, yy, f[[index.max]], nlevels = 1, levels = C, col = col.threshold, add = TRUE, drawlabels = FALSE, lwd = lwd.curve)
        
      }else{
        
        plot(x, y, 'p', col = ifelse(D, adjustcolor(col.controlscases[2],alpha.f = alpha.points), adjustcolor(col.controlscases[1], alpha.f = alpha.points)),
             pch = ifelse(Z > C, 16, 15), cex = cex, xaxs = "i", yaxs = "i", xaxt = ifelse(new,"n","s"), yaxt = ifelse(new,"n","s"),
             xlim = range(x,finite=TRUE)+lx*c(-1,1), ylim = range(y,finite=TRUE)+ly*c(-1,1),
             xlab = xlab, ylab = ylab, main = "Classification subsets", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
        
        points(x, y, col = ifelse(Z > C, col.threshold, 'gray'), cex = 1.3*cex, pch = ifelse(Z > C, 21, 22))
        if(border) contour(xx, yy, f, nlevels = 1, levels = C, col = col.threshold, add = TRUE, drawlabels = FALSE, lwd = lwd.curve)
        
      }
      
    }
    
    if(legends){
      legend('bottomright', inset = 0.03, c("Controls", "Cases"), col = col.controlscases, lty = 1, cex = 0.5)
      legend('topright', inset = 0.03, c("Classified as Controls","Classified as Cases"), fill = c('white', colTrans), cex = 0.5)
    }
    
    par(fig = c(0.51,1,0,1), new = TRUE)
    
    if(build.process){
      plot(c(0,t,1), c(0,roc,1), type = type, col = ifelse(completeROC,'gray','white'), xlim = c(0,1), ylim = c(0,1), xaxt = ifelse(new,"n","s"), yaxt = ifelse(new,"n","s"),
           lwd = lwd.curve, lty = lty.curve, xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
      lines(c(0,t[1:index.max]), c(0,roc[1:index.max]), type = type, lwd = lwd.curve, col = col.curve, lty = lty.curve)
    }else{
      plot(c(0,t,1), c(0,roc,1), type = type, lwd = lwd.curve, lty = lty.curve, col = col.curve, xlim = c(0,1), ylim = c(0,1), xaxt = ifelse(new,"n","s"), yaxt = ifelse(new,"n","s"),
           xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    }
    points(t[index.max], roc[index.max], pch = 16, col = col.threshold, cex = cex.point)
    if(cutoff) text(t[index.max]+0.025, roc[index.max], format(C,1,digits=2), pos = 1, col = 'gray30', cex = cex.point)
    
    abline(0,1, col='gray', lty = 2)
    axis(1, at = seq(0,1,0.01), labels = F, tck = -0.01)
    axis(1, at = seq(0.1,0.9,0.1), labels = F, tck = -0.02)
    axis(2, at = seq(0,1,0.01), labels = F, tck = -0.01)
    axis(2, at = seq(0.1,0.9,0.1), labels = F, tck = -0.02)
    if(legends) legend("bottomright", paste0("AUC=", round(auc,3)), inset = 0.02, bty = "n")
    
  }
  
}
