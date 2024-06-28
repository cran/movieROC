\name{plot_densities}
\alias{plot_densities}

\title{
  Plot density function estimates for controls and cases
}
\description{
This function plots the kernel density function estimates for controls and cases from an object of class \sQuote{groc}. It offers the possibility of plotting the histograms instead of the smooth estimates.
}
\usage{
plot_densities(obj, h = c(1, 1), histogram = FALSE, breaks = 15, 
  col = c("#485C99", "#8F3D52"), lwd = 2, xlim = NULL, ylim = NULL, 
  xaxs = "i", yaxs = "i", xlab = "Marker", ylab = "f(x)", 
  main = "Density functions", legend = FALSE, pos.legend = "topright", 
  cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, cex.legend = 1, 
  eps = sqrt(.Machine$double.eps), new = TRUE, ...)
}

\arguments{
  \item{obj}{
  An object of class \sQuote{groc}.}
  \item{h}{
  A vector of length 2 with the bandwidth used to compute kernel density estimation for controls and cases, respectively. See \code{adjust} parameter in \code{density()} function from \pkg{stats} package. Default: \code{h = c(1,1)}.}
  \item{histogram}{
  If TRUE, histograms are displayed instead of kernel density estimates. Default: FALSE.}
  \item{breaks}{
  If \code{histogram = TRUE}, number of breaks used for the histograms. Default: 15.}
  \item{col}{
  A vector of length 2 with color used for controls and cases, respectively. Default: \code{c("#485C99", "#8F3D52")}.}
  \item{lwd}{
  Line width for the density function or histogram. Default: 2.}
  \item{xlim, ylim}{
  Range for x- and y-axis.}
  \item{xaxs, yaxs}{
  The style of axis interval calculation to be used for the x- and the y-axis, respectively. Default: \code{"i"} (internal; just finds an axis with pretty labels that fits within the original data range). For more information, see help for \code{par}.}
  \item{xlab, ylab}{
  Label for x- and y-axis.}
  \item{main}{
  Title for the plot.}
  \item{cex.lab, cex.axis, cex.main}{
  The magnification to be used for labels, axis annotation and main titles, respectively, relative to the current setting of \code{cex}. Default: \code{cex.lab = 1.5, cex.axis = 1.5, cex.main = 2}.}
  \item{legend}{
  If TRUE, a legend with the meaning of colors is displayed on the left plot. Default: FALSE.}
  \item{pos.legend, cex.legend}{
  The position and magnification to be used for legend, relative to the current setting of \code{cex}. Default: \code{pos.legend = "topright", cex.legend = 1}.}
  \item{eps}{
  Epsilon value to consider. Default: \code{sqrt(.Machine$double.eps)}.}
  \item{new}{
  If TRUE, a new plot is displayed; otherwise, density estimates are displayed over the current plot. Default: TRUE.}
  \item{\dots}{
  Other parameters to be passed. Not used.}
}
\value{
A plot of the kernel density function estimates for both populations with the selected graphical parameters
}
\examples{
data(HCC)
roc_cg20202438 <- gROC(X = HCC$cg20202438, D = HCC$tumor)
plot_densities(roc_cg20202438, main = "Density functions and histogram")
plot_densities(roc_cg20202438, histogram = TRUE, new = FALSE)
}
