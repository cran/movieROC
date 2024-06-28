\name{plot_funregions}
\alias{plot_funregions}
\alias{plot_funregions.hroc}
\alias{plot_funregions.groc}
\title{
  Plot the transformation function used for the marker
}
\description{
  This function plots the transformation function used for the marker for an object of class \sQuote{hroc} or class \sQuote{groc} (this one is only allowed for objects with self-contained classification subsets).
  Over this graphic, the classification region for a particular FPR chosen by the user is displayed in blue color.
  If the user specifies a second FPR (\code{FPR2}), the classification region is displayed on the same graphic in red color.
}
\usage{
\S3method{plot_funregions}{hroc}(x, FPR = 0.15, FPR2 = NULL, 
  plot.subsets = TRUE, new.window = FALSE, main = NULL, ylim = NULL, ...)
\S3method{plot_funregions}{groc}(x, FPR = 0.15, FPR2 = NULL, 
  plot.subsets = TRUE, new.window = FALSE, main = NULL, ylim = NULL, ...)
}

\arguments{
  \item{x}{
  An object of class \sQuote{hroc} or \sQuote{groc}.}
  \item{FPR}{
  False-positive rate used to display the classification region in blue. Default: 0.15.}
  \item{FPR2}{
  Other false-positive rate used to display the classification region in red. Default: none.}
  \item{plot.subsets}{
  If TRUE, the classification subsets are displayed. Otherwise, only the transformation function. Default: TRUE.}
  \item{new.window}{
  If TRUE, graphics are displayed in a new window. Default: FALSE.}
  \item{main}{
  Title for the plot.}
  \item{ylim}{
  Range for the y-axis.}
  \item{\dots}{
  Other parameters to be passed. Not used.}
}
\value{
A plot of the transformation function used for the marker with the selected graphical parameters
}
\examples{
data(HCC)

hroc_cg18384097 <- hROC(X = HCC$cg18384097, D = HCC$tumor, 
    formula.lrm = "D ~ rcs(X,8)") 
plot_funregions(hroc_cg18384097)
plot_funregions(hroc_cg18384097, FPR = .1, FPR2 = .5)
}