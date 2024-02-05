plot.regions <- function(x, ...){
  UseMethod("plot.regions")
}

plot.regions.groc <- function(x, FPR = 0.15, plot.roc = TRUE, plot.auc = FALSE, col = c('white','grey'), col.FPR = 'blue', lwd = 2, new.window = TRUE, type.plotroc = 's', xlim = NULL, mar = c(5,6,4,0.25),
                              cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.75, main = NULL, xlab = "", ylab = "False-Positive Rate", main.plotroc = "ROC curve", legend = TRUE, cex.legend = 1, ...){

  obj <- x
  side <- obj$side
  if(is.null(main)) main <- ifelse(obj$param, "Classification subsets [Parametric]", "Classification subsets [Non-Parametric]")

  if(is.null(xlim)){
    inf <- min(c(obj$controls, obj$cases)); sup <- max(c(obj$controls, obj$cases))
  }else{
    inf <- xlim[1]; sup <- xlim[2]
  }
  xlim <- c(inf,sup)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if(new.window){
    par(oma=c(0.2,0.5,0.2,0.5))
    if(plot.roc) layout(rbind(c(1,1,1,2)))
  }
  if(plot.roc) par(mar=mar)

  colTrans <- col
  colTrans[1] <- rgb(red=col2rgb(col[1])[1], green=col2rgb(col[1])[2], blue=col2rgb(col[1])[3], alpha=0.25*255, maxColorValue=255)
  colTrans[2] <- rgb(red=col2rgb(col[2])[1], green=col2rgb(col[2])[2], blue=col2rgb(col[2])[3], alpha=0.25*255, maxColorValue=255)

  if(side=='right' || side=='left'){

    if(length(obj$t) > 151){
      index.t <- sapply(seq(0,1,length.out=151), function(fpr){which.min(abs(obj$t - fpr))})
      obj$c <- obj$c[index.t]
      obj$t <- obj$t[index.t]
      obj$roc <- obj$roc[index.t]
    }

    plot(obj$c, 1-obj$t, xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = c(0,1), yaxt = 'n', xaxs = 'i', col = 'white', cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis)
    polygon(c(rep(inf,length(obj$t)), rev(obj$c)), c(1-obj$t, rev(1-obj$t)), col = ifelse(side=='right',colTrans[1], colTrans[2]), border = NA)
    polygon(c(obj$c, rev(rep(sup,length(obj$t)))), c(1-obj$t, rev(1-obj$t)), col = ifelse(side=='right',colTrans[2], colTrans[1]), border = NA)
    segments(rep(inf,length(obj$t)), 1-obj$t, obj$c, 1-obj$t, col = ifelse(side=='right', col[1], col[2]))
    segments(obj$c, 1-obj$t, rep(sup,length(obj$t)), 1-obj$t, col = ifelse(side=='right', col[2], col[1]))

    ticks.axis1 <- axis(1, cex.axis = cex.axis)
    space0 <- max(diff(ticks.axis1))
    axis(1, at = seq(min(ticks.axis1)-space0, max(ticks.axis1)+space0, space0/4), tcl = -0.4, labels = FALSE)
    axis(1, at = seq(min(ticks.axis1)-space0, max(ticks.axis1)+space0, space0/8), tcl = -0.3, labels = FALSE)

    axis(2, at = seq(1,0,-0.1), labels = F, tck = -0.02, cex.axis = cex.axis)
    axis(2, at = seq(0,1,0.5), labels = seq(1,0,-0.5), tck = -0.04, cex.axis = cex.axis)
    axis(4, at = c(0,1), labels = c("",""), tck = 0, cex.axis = cex.axis)

    if(!is.null(FPR)){
      info <- predict(obj, FPR=FPR)
      info$ClassSubsets <- ifelse(info$ClassSubsets == -Inf, inf, ifelse(info$ClassSubsets == Inf, sup, info$ClassSubsets))
      arrows(info$ClassSubsets[1], info$Specificity, info$ClassSubsets[2], info$Specificity, col = col.FPR, angle = 75, code = 3, length = 0.03)
    }

    if(legend) legend('topleft', obj$levels, pch = 22, col='black', pt.bg = col, title = "Classification:", inset = 0.01, bty = 'n', cex = cex.legend)
  }
  if(side=='both' || side=='both2'){

    if(length(obj$t) > 151){
      index.t <- sapply(seq(0,1,length.out=151), function(fpr){which.min(abs(obj$t - fpr))})
      obj$xl <- obj$xl[index.t]
      obj$xu <- obj$xu[index.t]
      obj$t <- obj$t[index.t]
      obj$roc <- obj$roc[index.t]
    }

    plot(obj$xl, 1-obj$t, xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = c(0,1), yaxt = 'n', xaxs = 'i', col = 'white', cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis)
    polygon(c(obj$xl, rev(obj$xu)), c(1-obj$t, rev(1-obj$t)), col=ifelse(side=='both2',colTrans[2], colTrans[1]), border=NA)
    polygon(c(rep(inf,length(obj$t)), rev(obj$xl)), c(1-obj$t, rev(1-obj$t)), col=ifelse(side=='both2',colTrans[1], colTrans[2]), border=NA)
    polygon(c(obj$xu, rev(rep(sup,length(obj$t)))), c(1-obj$t, rev(1-obj$t)), col=ifelse(side=='both2',colTrans[1], colTrans[2]), border=NA)
    segments(obj$xl, 1-obj$t, obj$xu, 1-obj$t, col=ifelse(side=='both2',col[2], col[1]))
    segments(rep(inf,length(obj$t)), 1-obj$t, obj$xl, 1-obj$t, col=ifelse(side=='both2',col[1], col[2]))
    segments(obj$xu, 1-obj$t, rep(sup,length(obj$t)), 1-obj$t, col=ifelse(side=='both2',col[1], col[2]))

    ticks.axis1 <- axis(1, cex.axis = cex.axis)
    space0 <- max(diff(ticks.axis1))
    axis(1, at =seq(min(ticks.axis1)-space0, max(ticks.axis1)+space0, space0/4), tcl = -0.4, labels = FALSE)
    axis(1, at =seq(min(ticks.axis1)-space0, max(ticks.axis1)+space0, space0/8), tcl = -0.3, labels = FALSE)

    axis(2, at=seq(1,0,-0.1), labels=F, tck=-0.02, cex.axis = cex.axis)
    axis(2, at=seq(0,1,0.5), labels=seq(1,0,-0.5), tck=-0.04, cex.axis = cex.axis)
    axis(4, at=c(0,1), labels=c("",""), tck=0, cex.axis = cex.axis)

    if(!is.null(FPR)){
      info <- predict(obj, FPR=FPR)
      info$ClassSubsets <- ifelse(info$ClassSubsets == -Inf, inf, ifelse(info$ClassSubsets == Inf, sup, info$ClassSubsets))
      if(length(info$ClassSubsets) == 2){info$ClassSubsets <- matrix(info$ClassSubsets, nrow = 1)}
      arrows(info$ClassSubsets[,1], info$Specificity, info$ClassSubsets[,2], info$Specificity, col = col.FPR, angle = 75, code = 3, length = 0.03)
    }

    if(legend) legend('topleft', obj$levels, pch = 22, col='black', pt.bg = col, title = "Classification:", inset = 0.01, bty = 'n', cex = cex.legend)

  }

  if(plot.roc){
    par(mar=c(par("mar")[1],1,par("mar")[3],5))
    plot(c(0,obj$roc,1),c(1,1-obj$t,0), type = type.plotroc,main=" ",xlab=" ", yaxt="n", cex.lab=cex.lab, cex.main=cex.main, cex.axis=cex.axis, xaxt="n", xlim=c(0,1), ylim=c(0,1), ylab=" ", lwd=lwd)
    lines(c(0,1),c(1,0),lty=2)
    text(x=1.15, y=0.5, labels = main.plotroc, srt=-90, xpd=TRUE,font=2, cex = cex.lab)
    axis(1,at=c(0,0.5,1),labels=c(0,0.5,1), cex.axis = cex.axis)
    axis(1,xaxp=c(0,1,40),tcl=-0.2,tcl=-0.2,labels=F, cex.axis = cex.axis)
    axis(3,at=0.5,labels="TPR",tcl=0, cex.axis = cex.axis)
    if(!is.null(FPR)){
      index.FPR <- which.min(abs(obj$t - FPR))
      lines(x = c(0,obj$roc[index.FPR]), y = c(1-obj$t[index.FPR],1-obj$t[index.FPR]), lty = 3, col = col.FPR)
      lines(x = c(obj$roc[index.FPR],obj$roc[index.FPR]), y = c(0,1-obj$t[index.FPR]), lty = 3, col = col.FPR)
      points(obj$roc[index.FPR], 1-obj$t[index.FPR], pch = 16, col = col.FPR, cex = 1.5)
    }
    if(plot.auc) legend('bottomleft', paste("AUC=",round(obj$auc,3),sep=''), cex = 0.75*cex.axis, bty='n', inset=0.01)
  }

}



plot.regions.hroc <- function(x, FPR = 0.15, plot.roc = TRUE, plot.auc = FALSE, col = c('white','grey'), col.FPR = 'blue', lwd = 2, new.window = TRUE, type.plotroc = 's', xlim = NULL, mar = c(5,6,4,0.25),
                              cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.75, main = NULL, xlab = "", ylab = "False-Positive Rate", main.plotroc = "ROC curve", legend = TRUE, cex.legend = 1, verbose = FALSE, ...){

  obj <- x
  X <- obj$X; Y <- obj$Y; Sp <- obj$Sp; Se <- obj$Se; type <- obj$type
  indexX <- order(X)

  if(is.null(xlim)){
    inf <- min(X); sup <- max(X)
  }else{
    inf <- xlim[1]; sup <- xlim[2]
  }
  xlim <- c(inf,sup)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if(new.window){
    par(oma = c(0.2,0.5,0.2,0.5))
    if(plot.roc) layout(rbind(c(1,1,1,2)))
  }
  if(plot.roc) par(mar = mar)

  colTrans <- col
  colTrans[1] <- rgb(red = col2rgb(col[1])[1], green = col2rgb(col[1])[2], blue = col2rgb(col[1])[3], alpha = 0.25*255, maxColorValue = 255)
  colTrans[2] <- rgb(red = col2rgb(col[2])[1], green = col2rgb(col[2])[2], blue = col2rgb(col[2])[3], alpha = 0.25*255, maxColorValue = 255)

  plot(X[indexX], Sp[indexX], 'l', xlim = xlim, ylim = c(0,1), xlab = xlab, ylab = ylab, yaxt = 'n', xaxs = 'i', col = colTrans[1],
       cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, main = ifelse(is.null(main), paste("Classification subsets:", type, ifelse(type=='lrm', obj$formula, "")), main))

  uni.Sp <- unique(Sp)
  if(length(uni.Sp) > 151) uni.Sp <- seq(min(Sp), 1, length.out = 151)
  Segments <- NULL
  if(verbose) bar <- txtProgressBar(min = 0, max = length(uni.Sp), style = 3, title = "Classification subsets for each FPR (progress bar)")
  for(i in 1:length(uni.Sp)){
    subsets <- predict(obj, FPR = 1-uni.Sp[i])$ClassSubsets
    subsets <- ifelse(subsets == -Inf, inf, ifelse(subsets == Inf, sup, subsets))
    Segments <- rbind(Segments, cbind(subsets, rep(uni.Sp[i], nrow(subsets))))
    if(verbose) setTxtProgressBar(bar, i)
  }
  if(verbose) close(bar)
  segments(Segments[,1], Segments[,3], Segments[,2], Segments[,3], col = col[2])

  if(!is.null(FPR)){
    index.FPR <- which.min(abs(Segments[,3] - (1-FPR)))
    index.FPR <- which(Segments[,3] == Segments[index.FPR,3])
    arrows(Segments[index.FPR,1], Segments[index.FPR,3], Segments[index.FPR,2], Segments[index.FPR,3], col = col.FPR, angle = 75, code = 3, length = 0.03)
  }

  ticks.axis1 <- axis(1, cex.axis = cex.axis)
  space0 <- max(diff(ticks.axis1))
  axis(1, at =seq(min(ticks.axis1)-space0, max(ticks.axis1)+space0, space0/4), tcl = -0.4, labels = FALSE)
  axis(1, at =seq(min(ticks.axis1)-space0, max(ticks.axis1)+space0, space0/8), tcl = -0.3, labels = FALSE)

  axis(2, at=seq(1,0,-0.1), labels=F, tck=-0.02, cex.axis = cex.axis)
  axis(2, at=seq(0,1,0.5), labels=seq(1,0,-0.5), tck=-0.04, cex.axis = cex.axis)
  axis(4, at=c(0,1), labels=c("",""), tck=0, cex.axis = cex.axis)
  if(legend) legend('topleft', obj$levels, pch = 22, col = 'black', pt.bg = col, title = "Classification:", inset = 0.01, bty = 'n', cex = cex.legend)

  if(plot.roc){
    par(mar = c(par("mar")[1],1,par("mar")[3],5))
    plot(c(1,Se[order(Y)],0),c(0,Sp[order(Y)],1), type = type.plotroc, main = " ", xlab = " ", yaxt = "n",
         cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, xaxt = "n", xlim = c(0,1), ylim = c(0,1), ylab = " ", lwd = lwd)
    lines(c(0,1), c(1,0), lty = 2)
    text(x = 1.15, y = 0.5, labels = main.plotroc, srt = -90, xpd = TRUE, font = 2, cex = cex.lab)
    axis(1, at = c(0,0.5,1), labels = c(0,0.5,1), cex.axis = cex.axis)
    axis(1, xaxp = c(0,1,40), tcl = -0.2, labels = F, cex.axis = cex.axis)
    axis(3, at = 0.5, labels = "TPR", tcl = 0, cex.axis = cex.axis)
    if(!is.null(FPR)){
      index.FPR <- which.min(abs(1 - Sp[order(Y)] - FPR))
      lines(x = c(0,Se[order(Y)][index.FPR]), y = c(Sp[order(Y)][index.FPR],Sp[order(Y)][index.FPR]), lty = 3, col = col.FPR)
      lines(x = c(Se[order(Y)][index.FPR],Se[order(Y)][index.FPR]), y = c(0,Sp[order(Y)][index.FPR]), lty = 3, col = col.FPR)
      points(Se[order(Y)][index.FPR], Sp[order(Y)][index.FPR], pch = 16, col = col.FPR, cex = 1.5)
    }
    if(plot.auc) legend('bottomleft', paste0("AUC=",round(obj$auc,3)), cex = 0.75*cex.axis, bty = 'n', inset = 0.01)
  }

}
