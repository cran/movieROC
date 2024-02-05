\name{plot.buildROC}
\alias{plot.buildROC}
\alias{plot.buildROC.groc}
\alias{plot.buildROC.biroc}
\alias{plot.buildROC.multiroc}

\title{
  Plot the building procedure of the ROC curve
}

\description{
  This function tracks the ROC curve. It plots two graphics in the same figure:
  right, the empirical ROC curve estimate; left, classification subset for a 
  particular FPR or threshold(s) for the marker (threshold only for object of class \sQuote{groc}).
}
\usage{
\S3method{plot}{buildROC}(x, ...)
\S3method{plot.buildROC}{groc}(x, FPR = NULL, C = NULL, XL = NULL, XU = NULL, h = c(1,1), 
    histogram = FALSE, breaks = 15, reduce = TRUE, build.process = FALSE, 
    completeROC = FALSE, new.window = FALSE, legends = FALSE, type = 's', 
    cex.point = 1.5, lwd.curve = 2, mar = NULL, lim.density = 0.01, xlim = NULL, 
    ylim = NULL, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, xlab = "Marker", 
    main.density = "Density functions", col.controlscases = c('#485C99','#8F3D52'), 
    col.threshold = '#FCBA04', col.curve = 'black', eps = sqrt(.Machine$double.eps), 
    ...)
\S3method{plot.buildROC}{biroc}(x, FPR = 0.15, build.process = FALSE, completeROC = TRUE, 
    new = FALSE, new.window = FALSE, border = TRUE, cutoff = TRUE, 
    legends = FALSE, type = 's', col.controlscases = c('#485C99','#8F3D52'), 
    col.threshold = '#FCBA04', col.curve = 'black', cex.point = 1.5, 
    alpha.points = .75, alpha.contour = 0.25, lwd.curve = 2, lty.curve = 1, 
    cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, xlab = "X1", ylab = "X2", 
    lf = NULL, eps = sqrt(.Machine$double.eps), ...)
\S3method{plot.buildROC}{multiroc}(x, FPR = 0.15, display.method = c("PCA", "OV"), 
    displayOV = c(1,2), build.process = FALSE, completeROC = TRUE, 
    new = FALSE, new.window = FALSE, border = FALSE, cutoff = TRUE, legends = FALSE, 
    type = 's', col.controlscases = c('#485C99','#8F3D52'), 
    col.threshold = '#FCBA04', col.curve = 'black', cex.point = 1.5, 
    alpha.points = .75, alpha.contour = 0.25, lwd.curve = 2, lty.curve = 1, 
    cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, xlab = NULL, ylab = NULL, 
    lf = NULL, eps = sqrt(.Machine$double.eps), ...)
}
\arguments{
  \item{x}{
  An ROC curve object from movieROC package. Possible classes are: \sQuote{groc} (output of \code{gROC} function), \sQuote{biroc} (output of \code{biROC} function), and \sQuote{multiroc} (output of \code{multiROC} function).}
  \item{FPR}{
  False-positive rate for which the left plot is computed. Default: 0.15.}
  \item{C}{
  Marker cutoff for which the left plot is computed. Only available for a \sQuote{groc} object.}
  \item{XL, XU}{
  Marker cutoffs for which the left plot is computed. Only available for object of class \sQuote{groc} with \code{side = "both"} or \code{"both2"}.}
  \item{h}{
  A vector of length 2 with the bandwidth used to compute kernel density estimation for controls and cases, respectively. See \code{adjust} parameter in \code{density} function from \pkg{stats} package. Default: \code{h = c(1,1)}.}
  \item{histogram}{
  If TRUE, histograms are displayed on the left instead of kernel density estimates. Default: FALSE. Only available for a \sQuote{groc} object.}
  \item{breaks}{
  If \code{histogram = TRUE}, number of breaks used for the histograms. Default: 15.}
  \item{reduce}{
  If FALSE, two extra graphics are displayed at the bottom (see Details for more information). Default: TRUE. Only available for a \sQuote{groc} object.}
  \item{build.process}{
  If FALSE, the whole ROC curve is displayed; otherwise, if \code{completeROC = TRUE}, the portion of the ROC curve until the fixed FPR is highlighted in black and the rest is shown in gray. Default: FALSE.}
  \item{completeROC}{
  If \code{build.process = TRUE}, a logical value indicating if the whole ROC curve should be displayed in light gray or not. Default: FALSE for \sQuote{groc} object; TRUE for \sQuote{biroc} and \sQuote{multiroc} objects.}
  \item{new.window, new}{
  If TRUE, graphics are displayed in a new window. Default: FALSE. \code{new} not available for object of class \sQuote{groc}.}
  \item{legends}{
  If TRUE, a legend with colors meaning is displayed on the left graphic and the AUC is shown on the right graphic. Default: FALSE.}
  \item{type, lty.curve, lwd.curve, col.curve}{
  The line type, width and color for ROC curve. Default: \code{type = "s"} (stair steps), \code{lty.curve = 1} (solid), \code{lwd.curve = 2} and \code{col.curve = "black"}.}
  \item{cex.lab, cex.axis, cex.main}{
  The magnification to be used for labels, axis annotation and main titles, respectively, relative to the current setting of \code{cex}. Default: \code{cex.lab = 1.5, cex.axis = 1.5, cex.main = 2}.}
  \item{cex}{
  A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. Default: 0.8. Not available for object of class \sQuote{groc} (\code{cex=1} in this case).}
  \item{cex.point}{
   The magnification to be used for the particular point on the ROC curve and its text, relative to the current setting of \code{cex}. Default: \code{cex.point = 1.5}.}
  \item{xlim, ylim}{
  Range for x- and y-axis on the left plot. Only available for a \sQuote{groc} object.}
  \item{lim.density}{
  If \code{xlim} is \code{NULL}, lower limit for the density estimate to compute the range for the x-axis on the left plot. Default: 0.01. Only available for a \sQuote{groc} object.}
  \item{xlab, ylab}{
  Label for x- and y-axis on the left plot. \code{ylab} not available for object of class \sQuote{groc}.}
  \item{main.density}{
  Title for the left plot. Only available for a \sQuote{groc} object.}
  \item{col.controlscases}{
  Vector of length 2 with the color used to control and case group, respectively. Default: \code{c('#485C99','#8F3D52')}.}
  \item{col.threshold}{
  Color for the chosen point. Default: \code{"FCBA04"}.}
  \item{mar}{
  A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. Only available for a \sQuote{groc} object.}
  \item{eps}{
  Epsilon value to consider. Default: \code{sqrt(.Machine$double.eps)}.}
  \item{border}{
  If TRUE, a border for the classification subsets is drawn. Default: TRUE for \sQuote{biroc} object and FALSE for \sQuote{multiroc} object. Not available for object of class \sQuote{groc}.}
  \item{cutoff}{
  If TRUE, the cutoff value for the resulting univariate marker is displayed on the right plot (ROC curve). Default: TRUE. Not available for object of class \sQuote{groc}.}
  \item{alpha.points, alpha.contour}{
  Number in [0,1] modifying the opacity alpha of the color for the points and classification region. The \code{adjustcolor} function in \pkg{grDevices} package is used. Default: \code{alpha.points = 0.75, alpha.contour = 0.25}. Not available for object of class \sQuote{groc}.}
  \item{lf}{
  Epsilon value for steps. Not used. Not available for object of class \sQuote{groc}.}
  \item{display.method}{
  Method to display the marker values from a \sQuote{multiroc} object on the left plot. Methods available: \code{"OV"} (projected over two selected components of the marker indicated in \code{displayOV}), or \code{"PCA"} (projected over the two principal components from a Principal Component Analysis). Default: \code{"PCA"}. Only available for object of class \sQuote{multiroc}.}
  \item{displayOV}{
  If \code{display.method = "OV"}, the two components of the marker used to project the marker values on the left plot. Default: \code{c(1,2)} (two first components). Only available for object of class \sQuote{multiroc}.}
  \item{\dots}{
  Other parameters to be passed. Not used.}
}

\details{
- For univariate ROC curves:

It plots two graphics in the same figure:
  left, density function estimates for the marker in both populations with the areas 
  corresponding to FPR and TPR colored (blue and red, respectively) for a particular 
  FPR or threshold(s) for the marker; right, the empirical ROC curve estimate. 
  
For an object of class \sQuote{groc}, if the input parameter \code{reduce = FALSE}, 
two more graphics are displayed at the bottom: 
left, boxplots for the marker distribution in both populations;  
right, classification subsets for every false-positive rate (grey color by default).

- For multivariate ROC curves:

Plot two graphics in the same figure: right, the ROC curve highlighting the point and the
threshold for the resulting univariate marker; left, scatterplot with the marker values for controls (blue color) and cases (red color), and the classification subset (in gold color) reporting the false-positive rate selected by the user (if \code{FPR} is not \code{NULL}).

If \code{build.process = FALSE}, the whole ROC curve is displayed; otherwise, if \code{completeROC = TRUE}, 
the portion of the ROC curve until the fixed FPR is highlighted in black and the rest is
shown in gray, while if \code{completeROC = FALSE}, only the first portion of the curve is illustrated.
}
\value{
A plot of the building procedure of the ROC curve with the selected graphical parameters
}
\examples{
data(HCC)

# Standard ROC curve for gene 20202438
roc_cg20202438 <- gROC(X = HCC$cg20202438, D = HCC$tumor, side = "right")
plot.buildROC(roc_cg20202438)
plot.buildROC(roc_cg20202438, C = .77, build.process = TRUE, reduce = FALSE)

# Multivariate ROC curve for genes 0202438, 18384097, and 03515901
multiroc_PT <- multiROC(X = cbind(HCC$cg20202438, HCC$cg18384097, HCC$cg03515901), 
    D = HCC$tumor, method = "fixedLinear", methodLinear = "PepeThompson")
plot.buildROC(multiroc_PT, cex = 1.2, lwd.curve = 4) # Two first components from PCA
plot.buildROC(multiroc_PT, display.method = "OV", displayOV = c(1,3), 
    cex = 1.2, xlab = "Gene 20202438", ylab = "Gene 03515901", lwd.curve = 4)
}
