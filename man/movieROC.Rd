\name{movieROC}
\alias{movieROC}
\alias{movieROC.groc}
\alias{movieROC.biroc}
\alias{movieROC.multiroc}

\title{
  Create a video with the building procedure of the ROC curve
}
\description{
  Save a video as a GIF illustrating the construction of the ROC curve. Particularly, for each
element in vector \code{fpr} (optional input parameter), the function executed is

- for \sQuote{groc} and \sQuote{biroc} objects:

  \code{plot.buildROC(obj, FPR = fpr[i], build.process = TRUE,...)};

- for \sQuote{multiroc} objects:

  \code{plot.buildROC(obj, FPR = fpr[i], build.process = TRUE, display.method, displayOV,...)}
}
\usage{
movieROC(x, ...)
\S3method{movieROC}{groc}(x, fpr = NULL, h = c(1,1), histogram = FALSE, breaks = 15,
  reduce = TRUE, completeROC = FALSE, videobar = TRUE, file = "animation1.gif",
  save = TRUE, legends = FALSE, speedcorrection = FALSE, tpause = 1, interval = 0.2,
  ani.width, ani.height, xlab = "Marker", main.density = "Density functions",
  cex.lab = 2.5, cex.axis = 1.75, cex.main = 2.25 + as.numeric(reduce),
  xlim = NULL, ylim = NULL, cex.point = 1.5, lwd.curve = 2, mar = NULL,
  lim.density = 0.01, col.controlscases = c("#485C99", "#8F3D52"),
  col.curve = "black", col.threshold = "#FCBA04", verbose = FALSE, ...)
\S3method{movieROC}{biroc}(x, fpr = NULL, border = TRUE, completeROC = FALSE,
  videobar = TRUE, file = "animation1.gif", save = TRUE, legends = FALSE,
  tpause = 1, interval = 0.2, ani.width, ani.height, xlab = "X1", ylab = "X2",
  cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,
  alpha.points = 1, alpha.contour = 0.25, lwd.curve = 2, lty.curve = 1,
  lf = NULL, col.controlscases = c('#485C99','#8F3D52'),
  col.curve = 'black', col.threshold = '#FCBA04', verbose = FALSE, ...)
\S3method{movieROC}{multiroc}(x, fpr = NULL, display.method = c("PCA", "OV"),
  displayOV = c(1,2), border = TRUE, completeROC = FALSE, videobar = TRUE,
  file = "animation1.gif", save = TRUE, legends = FALSE, tpause = 1,
  interval = 0.2, ani.width, ani.height, xlab = NULL, ylab = NULL,
  cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, alpha.points = 1,
  alpha.contour = 0.25, lwd.curve = 2, lty.curve = 1, lf = NULL,
  col.controlscases = c('#485C99','#8F3D52'), col.curve = 'black',
  col.threshold = '#FCBA04', verbose = FALSE, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  An ROC curve object from the \pkg{movieROC} package. Possible classes are: \sQuote{groc} (output of \code{gROC} function), \sQuote{biroc} (output of \code{biROC} function), and \sQuote{multiroc} (output of \code{multiROC} function).}
  \item{fpr}{
  A vector of false-positive rates for which the building process of the ROC curve is displayed. Default: if length of \code{x$t} is lower than 150, such vector is taken as \code{fpr}; otherwise, an equally-space vector of length 100 covering the range of the marker is considered.}
  \item{h}{
  A vector of length 2 with the bandwidth used to compute kernel density estimation for controls and cases, respectively. See \code{adjust} parameter in \code{density} function from \pkg{stats} package. Default: \code{h = c(1,1)}.}
  \item{histogram}{
  If TRUE, histograms are displayed on the left instead of kernel density estimates. Default: FALSE. Only available for a \sQuote{groc} object.}
  \item{breaks}{
  If \code{histogram = TRUE}, number of breaks used for the histograms. Default: 15.}
  \item{reduce}{
  If FALSE, two extra graphics are displayed at the bottom (see Details for more information). Default: TRUE. Only available for a \sQuote{groc} object.}
  \item{completeROC}{
  A logical value indicating if the whole ROC curve should be displayed in light gray or not. Default: FALSE.}
  \item{videobar}{
  If TRUE, a text progress bar is shown in the R console. Default: TRUE.}
  \item{verbose}{
  If TRUE, a progress bar is displayed for computationally intensive methods. Default: FALSE.}
  \item{file}{
  File name of the movie (with the extension). Default: \code{"animation1.gif"}.}
  \item{save}{
  If TRUE, video is saved as a GIF by using the \code{saveGIF} function in \CRANpkg{animation} package. Default: TRUE.}
  \item{legends}{
  If TRUE, a legend with colors meaning is displayed on the left graphic and the AUC is shown on the right graphic. Default: FALSE.}
  \item{speedcorrection}{
  If TRUE, only some FPR in \code{fpr} vector are considered to make the video faster. Default: FALSE. Only available for a \sQuote{groc} object.}
  \item{tpause}{
  If \code{save = FALSE}, time interval to suspend execution for, in seconds. Default: 1.}
  \item{interval}{
  If \code{save = TRUE}, a positive number to set the time interval of the animation (unit in seconds) in \CRANpkg{animation} package. Default: 0.2.}
  \item{ani.width, ani.height}{
  If \code{save = TRUE}, width and height of image frames (unit in px) in \CRANpkg{animation} package.}
  \item{xlab, ylab}{
  Label for x- and y-axis on the left plot. \code{ylab} not available for object of class \sQuote{groc}.}
  \item{main.density}{
  Title for the left plot. Only available for a \sQuote{groc} object.}
  \item{cex.lab, cex.axis, cex.main}{
  The magnification to be used for labels, axis annotation and main titles, respectively, relative to the current setting of \code{cex}.}
  \item{cex}{
  A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. Default: 0.8. Not available for object of class \sQuote{groc} (\code{cex=1} in this case).}
  \item{cex.point}{
   The magnification to be used for the particular point on the ROC curve and its text, relative to the current setting of \code{cex}. Default: \code{cex.point = 1.5}. Only available for a \sQuote{groc} object.}
  \item{xlim, ylim}{
  Range for x- and y-axis on the left plot. Only available for a \sQuote{groc} object.}
  \item{lim.density}{
  If \code{xlim} is \code{NULL}, lower limit for the density estimate to compute the range for the x-axis on the left plot. Default: 0.01. Only available for a \sQuote{groc} object.}
  \item{lty.curve, lwd.curve, col.curve}{
  The line type, width and color for ROC curve. Default: \code{lty.curve = 1} (solid), \code{lwd.curve = 2} and \code{col.curve = "black"}.}
  \item{mar}{
  A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. Only available for a \sQuote{groc} object.}
  \item{col.controlscases}{
  Vector of length 2 with the color used to control and case group, respectively. Default: \code{c('#485C99','#8F3D52')}.}
  \item{col.threshold}{
  Color for the chosen point. Default: \code{"FCBA04"}.}
  \item{border}{
  If TRUE, a border for the classification subsets is drawn. Default: TRUE. Not available for object of class \sQuote{groc}.}
    \item{alpha.points, alpha.contour}{
  Number in [0,1] modifying the opacity alpha of the color for the points and classification region. The \code{adjustcolor} function in \pkg{grDevices} package is used. Default: \code{alpha.points = 1, alpha.contour = 0.25}. Not available for object of class \sQuote{groc}.}
  \item{lf}{
  Epsilon value for steps. Not used.}
  \item{display.method}{
  Method to display the marker values from a \sQuote{multiroc} object on the left plot. Methods available: \code{"OV"} (projected over two selected components of the marker indicated in \code{displayOV}), or \code{"PCA"} (projected over the two principal components from a Principal Component Analysis). Default: \code{"PCA"}. Only available for object of class \sQuote{multiroc}.}
  \item{displayOV}{
  If \code{display.method = "OV"}, the two components of the marker used to project the marker values on the left plot. Default: \code{c(1,2)} (two first components). Only available for object of class \sQuote{multiroc}.}
  \item{\dots}{
  Other parameters to be passed to the \code{saveGIF} function in \CRANpkg{animation} package.}
}
\value{
A video with the building procedure of the ROC curve estimate with the selected graphical parameters
}
\examples{
data(HCC)

# Standard ROC curve for gene 20202438
roc_cg20202438 <- gROC(X = HCC$cg20202438, D = HCC$tumor, side = "right")
\donttest{
### Do not change the \code{save} parameter to save the video in a GIF
### with the name "video_cg20202438":
movieROC(roc_cg20202438, file = "video_cg20202438.gif", save = FALSE)}

# Multivariate ROC curve for genes 0202438, 18384097, and 03515901
multiroc_PT <- multiROC(X = cbind(HCC$cg20202438, HCC$cg18384097, HCC$cg03515901),
    D = HCC$tumor, method = "fixedLinear", methodLinear = "PepeThompson")
\donttest{
# Two first components from PCA:
### Do not change the \code{save} parameter to save the video in a GIF
### with the name "video_multiroc_pca":
movieROC(multiroc_PT, file = "video_multiroc_pca.gif", save = FALSE)
}
\donttest{
### Do not change the \code{save} parameter to save the video in a GIF
### with the name "video_multiroc_orig":
movieROC(multiroc_PT, display.method = "OV", displayOV = c(1,3),
    file = "video_multiroc_orig.gif", cex = 1.2, xlab = "Gene 20202438",
    ylab = "Gene 03515901", lwd.curve = 4, save = FALSE)
  }
}
