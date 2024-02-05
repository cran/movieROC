\name{movieROC2_densities}
\alias{movieROC2_densities}

\title{
  Create a video with the building procedure of the smooth ROC curve estimate
}
\description{
This function tracks the construction of the standard ROC curve (right- or left-sided, depending on the \code{side} of the object of class \sQuote{groc} included) resulting from the kernel density function estimation for controls and cases. Four graphics are displayed: top-left, the kernel density estimates; top-right, the resulting ROC curve; bottom-left, boxplots and points for controls and cases and classification subset in gray color; bottom-right, classification subsets for every FPR until the current one. It makes use of the \code{\link{plot_densityROC}} function for each screenshot.
}
\usage{
movieROC2_densities(obj, h = c(1, 1), cut.off = NULL, completeROC = FALSE,
  legends = FALSE, videobar = TRUE, file = "animation1.gif", clean = FALSE,
  interval = 0.2, ani.width = 500, ani.height = 750, save = TRUE,
  tpause = 1, verbose = FALSE, ...)
}

\arguments{
  \item{obj}{
  An object of class \sQuote{groc} with \code{side = "right"} or \code{"left"}.}
  \item{h}{
  A vector of length 2 with the bandwidth used to compute kernel density estimation for controls and cases, respectively. See \code{adjust} parameter in \code{density} function from \pkg{stats} package. Default: \code{h = c(1,1)}.}
  \item{cut.off}{
  Vector with marker cutoffs for which the graphics are displayed. Default: if number of unique marker values is lower than 150, these are considered; otherwise, a equally-spaced grid of length 102 in the range of the marker is used.}
  \item{save}{
  If TRUE, video is saved as a GIF by using the \code{saveGIF} function in \CRANpkg{animation} package. Default: TRUE.}
  \item{completeROC}{
  A logical value indicating if the whole ROC curve should be displayed in light gray or not. Default: FALSE.}
  \item{legends}{
  If TRUE, legends with the meaning of colors are displayed. Default: FALSE.}
  \item{videobar}{
  If TRUE, a text progress bar is shown in the R console. Default: TRUE.}
  \item{file}{
  File name of the movie (with the extension). Default: \code{"animation1.gif"}.}
  \item{clean}{
  Whether to delete the individual image frames in \CRANpkg{animation} package. Default: FALSE.}
  \item{interval}{
  A positive number to set the time interval of the animation (unit in seconds) in \CRANpkg{animation} package. Default: 0.2.}
  \item{ani.width, ani.height}{
   Width and height of image frames (unit in px) in \CRANpkg{animation} package.}
  \item{tpause}{
  If \code{save = FALSE}, time interval to suspend execution for, in seconds. Default: 1.}
  \item{verbose}{
  If TRUE, a progress bar is displayed for computationally intensive methods. Default: FALSE.}
  \item{\dots}{
  Other parameters to be passed to the \code{saveGIF} function in \CRANpkg{animation} package.}
}
\value{
A video with the building procedure of the smooth ROC curve estimate with the selected graphical parameters
}
\examples{
data(HCC)

# Standard ROC curve for gene 20202438
roc_cg20202438 <- gROC(X = HCC$cg20202438, D = HCC$tumor)
\donttest{
### Do not change the \code{save} parameter to save the video in a GIF
### with the name "video_cg20202438_smooth":
movieROC2_densities(roc_cg20202438, file = "video_cg20202438_smooth.gif", save = FALSE)}
}
