\name{plot.regions}
\alias{plot.regions}
\alias{plot.regions.groc}
\alias{plot.regions.hroc}
\title{
  Plot the classification regions underlying a ROC curve
}
\description{
  This function plots the classification regions for univariate markers.
  It works for objects of class \sQuote{groc} and \sQuote{hroc}.
  Two graphics are displayed in the same figure: left, classification subsets 
  for every false-positive rate (grey color by default); right, 90º rotated ROC curve. 
}
\usage{
\S3method{plot}{regions}(x, ...)
\S3method{plot.regions}{groc}(x, FPR = 0.15, plot.roc = TRUE, plot.auc = FALSE, 
  col = c("white", "grey"), col.FPR = "blue", lwd = 2, new.window = TRUE, 
  type.plotroc = "s", xlim = NULL, mar = c(5, 6, 4, 0.25), 
  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.75, main = NULL, 
  xlab = "", ylab = "False-Positive Rate", main.plotroc = "ROC curve", 
  legend = TRUE, cex.legend = 1, ...)
\S3method{plot.regions}{hroc}(x, FPR = 0.15, plot.roc = TRUE, plot.auc = FALSE, 
  col = c('white','grey'), col.FPR = 'blue', lwd = 2, new.window = TRUE, 
  type.plotroc = 's', xlim = NULL, mar = c(5,6,4,0.25), 
  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.75, main = NULL, 
  xlab = "", ylab = "False-Positive Rate", main.plotroc = "ROC curve", 
  legend = TRUE, cex.legend = 1, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  An ROC curve object from a univariate marker. Possible classes are: \sQuote{groc} (output of \code{gROC} function) and \sQuote{hroc} (output of \code{hROC} function),}
  \item{FPR}{
  False-positive rate to be highlighted, both on the ROC curve (right plot) and classification subset (left plot). Default: 0.15.}
  \item{plot.roc}{
  If TRUE, the ROC curve is displayed. Default: TRUE.
  }
  \item{plot.auc}{
  If TRUE, the AUC is displayed on the plot of the bottomleft corner of the ROC curve plot. Default: FALSE.}
  \item{col}{
  Vector of length 2 with colors used for outside and inside of classification subsets, respectively. Default: \code{c("white", "grey")}.}
  \item{col.FPR}{
  Color used to highlight the \code{FPR} chose. Default: \code{"blue"}.}
  \item{type.plotroc, lwd}{
  Line type and width for the ROC curve. Default: \code{type.plotroc = "s", lwd = 2}.
  }
  \item{new.window}{
  If TRUE, graphics are displayed in a new window. Default: TRUE.}
  \item{xlim}{
  Range for x-axis on the left plot. Default: range of the marker values.}
  \item{mar}{
  A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot.}
  \item{cex.lab, cex.axis, cex.main}{
  The magnification to be used for labels, axis annotation and main titles, respectively, relative to the current setting of \code{cex}. Default: \code{cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.75}.}
  \item{xlab, ylab}{
  Label for x- and y-axis on the left plot.}
  \item{main, main.plotroc}{
  Title for the left and the right plot, respectively.}
  \item{legend}{
  If TRUE, a legend with the meaning of colors is displayed on the left plot. Default: TRUE.}
  \item{cex.legend}{
  The magnification to be used for legend, relative to the current setting of \code{cex}. Default: 1.}
  \item{verbose}{
  If TRUE, a progress bar is displayed. Default: FALSE.}
  \item{\dots}{
  Other parameters to be passed. Not used.}
}
\value{
A plot of the classification regions underlying a ROC curve with the selected graphical parameters
}
\examples{
data(HCC)
# 1. Standard ROC curve
roc_cg18384097 <- gROC(X = HCC$cg18384097, D = HCC$tumor)
plot.regions(roc_cg18384097, plot.auc = TRUE)
# 2. gROC curve
groc_cg18384097 <- gROC(X = HCC$cg18384097, D = HCC$tumor, side = "both")
plot.regions(groc_cg18384097, plot.auc = TRUE)
# 3. hROC curve with a restricted cubic splines transformation
hroc_cg18384097 <- hROC(X = HCC$cg18384097, D = HCC$tumor, 
    formula.lrm = "D ~ rcs(X,8)") 
plot.regions(hroc_cg18384097, plot.auc = TRUE)
}
