\name{plot}
\alias{plot}
\alias{plot.groc}
\alias{plot.hroc}
\alias{plot.multiroc}

\title{
  Plot an ROC curve
}
\description{
This is one of the core functions of the movieROC package. It displays the empirical ROC curve estimate from an object of class \sQuote{groc}, \sQuote{hroc}, or \sQuote{multiroc}.
}
\usage{
\S3method{plot}{groc}(x, xlim = c(0, 1), ylim = c(0, 1), lwd = 3, 
  xlab = "False-Positive Rate", ylab = "True-Positive Rate", main = "ROC curve", 
  cex.lab = 1.25, cex.main = 1.5, type = NULL, new = TRUE, ...)
\S3method{plot}{hroc}(x, type = 'S', xlim = c(0,1), ylim = c(0,1), 
  lwd = 3, xlab = "False-Positive Rate", ylab = "True-Positive Rate", 
  main = "ROC Curve", cex.lab = 1.25, cex.main = 1.5, new = TRUE, ...)
\S3method{plot}{multiroc}(x, ...)
}

\arguments{
  \item{x}{
  An ROC curve object from movieROC package. Possible classes are: \sQuote{groc} (output of \code{\link{gROC}()} function), \sQuote{hroc} (output of \code{\link{hROC}()} function), and \sQuote{multiroc} (output of \code{\link{multiROC}()} function).}
  \item{xlim, ylim}{
  Range for x- and y-axis. Default: unit interval.}
  \item{lwd}{
  Line width of the ROC curve. Default: 3.}
  \item{xlab, ylab}{
  Label for x- and y-axis.}
  \item{main}{
  Title for the plot.}
  \item{cex.lab, cex.main}{
  The magnification to be used for labels and main title, respectively, relative to the current setting of \code{cex}. Default: \code{cex.lab = 1.25, cex.main = 1.5}.}
  \item{type}{
  What type of plot should be drawn (see help from \code{plot} function in \pkg{base} package). Default: \code{type = "s"} (stair steps), except for \code{x} object from \code{gROC_param} function, in which case \code{type = "l"} (lines).}
  \item{new}{
  If TRUE, a new plot is displayed; otherwise, the ROC curve is plotted over the existing graphic. Default: TRUE.}
  \item{\dots}{
  Other graphical parameters to be passed.}
}
\value{
A plot of the ROC curve with the selected graphical parameters
}
\examples{
data(HCC)
# ROC curve estimates for gene 03515901 and response tumor
rroc <- gROC(X = HCC[,"cg03515901"], D = HCC$tumor) # Right-sided
lroc <- gROC(X = HCC[,"cg03515901"], D = HCC$tumor, side = "left") # Left-sided
hroc <- hROC(X = HCC[,"cg03515901"], D = HCC$tumor) # Transformed by a cubic polinomial

plot(rroc, lty = 2, frame = FALSE)
plot(lroc, new = FALSE)
plot(hroc, new = FALSE, col = "blue")
legend("topleft", legend = c("Right-sided", "Left-sided", "Transformed marker"),
    col = c("black", "black", "blue"), lty = c(1,2,1), lwd = 2, bty = "n")

# ROC curve estimate for genes 20202438 and 18384097 to simultaneously identify tumor
# by a logistic regression model with quadratic formula
biroc <- multiROC(X = cbind(HCC$cg20202438, HCC$cg18384097), D = HCC$tumor)
plot(biroc)
legend("bottomright", paste("AUC = ", format(biroc$auc, digits = 3)))

# ROC curve estimate for genes 20202438, 18384097 and 03515901 to simultaneously 
# identify tumor by a linear combinations with fixed parameters by Pepe and Thompson (2000)
multiroc <- multiROC(X = cbind(HCC$cg20202438, HCC$cg18384097, HCC$cg03515901), 
    D = HCC$tumor, method = "fixedLinear", methodLinear = "PepeThompson")
plot(multiroc)
legend("bottomright", paste("AUC = ", format(multiroc$auc, digits = 3)))
}
