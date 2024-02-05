movieROC2_densities <- function(obj, h = c(1,1), cut.off = NULL,
                                 completeROC = FALSE, legends = FALSE, videobar = TRUE, file = "animation1.gif", clean = FALSE,
                                 interval = 0.2, ani.width = 500, ani.height = 750, save = TRUE, tpause = 1, verbose = FALSE, ...){

  side <- obj$side
  if(side != "right" & side != "left") stop("This function only works for standard ROC curve (right-sided or left-sided).")

  movie <- function(obj, cut.off, completeROC, videobar, legends){

    X <- sort(c(obj$controls, obj$cases))
    if(is.null(cut.off)){
      range <- max(X) - min(X)
      if(length(unique(X)) < 150 ){
        cut.off <- c(min(X) - range/20, sort(unique(X)), max(X) + range/20)
      }else{
        cut.off <- c(min(X) - range/20, seq(min(X), max(X), length.out=100), max(X) + range/20)
      }
    }

    B <- length(cut.off)

    if(videobar){
      if(verbose){
        cat("\nProgress bar: Construction of GIF with ", B, " thresholds. \n", sep = "")
        bar <- txtProgressBar(min = 0, max = B, style = 3)
      }
    }

    sapply(1:length(cut.off), function(i){
      C <- cut.off[i]
      plot_densityROC(obj, side = side, C = C, h = h, build.process = TRUE, completeROC = completeROC, legends = legends)
      if(videobar){if(verbose) setTxtProgressBar(bar, i)}
    })
    if(videobar){if(verbose) close(bar)}
    if(!save)  Sys.sleep(tpause)
  }


  if(save){
    animation::saveGIF(movie(obj, cut.off = cut.off, completeROC = completeROC, legends = legends, videobar = videobar),
                       movie.name = file, img.name = "Rplot", clean = clean, interval = interval, ani.width = ani.width, ani.height = ani.height, ...)
  }else{
    movie(obj, cut.off = cut.off, completeROC = completeROC, legends = legends, videobar = videobar)
  }

}
