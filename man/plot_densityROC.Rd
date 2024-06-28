\name{plot_densityROC}
\alias{plot_densityROC}

\title{
Plot standard smooth ROC curve estimate
}
\description{
This function estimates the standard ROC curve (right- or left-sided, depending on the \code{side} of the object of class \sQuote{groc} included) resulting from the kernel density function estimation for controls and cases. Two graphics are displayed: left, the kernel density estimates; right, the resulting ROC curve. Two extra graphics may be shown at the bottom if the used introduces a value for the input parameter \code{C}: left, boxplots and points for controls and cases and corresponding classification subset in gray color; right, classification subsets for every FPR until the one corresponding to the chosen \code{C}.
}
\usage{
plot_densityROC(obj, h = c(1, 1), C = NULL, build.process = FALSE, 
  completeROC = TRUE, legends = FALSE, rel.tol = 0.001, 
  par.specify = FALSE, cex.lab = 1.5, cex.axis = 1.25, cex.main = 1.75, 
  lwd = 2, col = c("#485C99", "#8F3D52"), col.roc = "blue", ...)
}

\arguments{
  \item{obj}{
  An object of class \sQuote{groc} with \code{side = "right"} or \code{"left"}.}
  \item{h}{
  A vector of length 2 with the bandwidth used to compute kernel density estimation for controls and cases, respectively. See \code{adjust} parameter in \code{density()} function from \pkg{stats} package. Default: \code{h = c(1,1)}.}
  \item{C}{
  Marker cutoff for which the graphics are displayed. Default: none.}
  \item{build.process}{
  If FALSE, the whole ROC curve is displayed; otherwise, if \code{completeROC = TRUE}, the portion of the ROC curve until the fixed FPR (resulting from the chosen \code{C}) is highlighted in black and the rest is shown in gray. Default: FALSE.}
  \item{completeROC}{
  If \code{build.process = TRUE}, a logical value indicating if the whole ROC curve should be displayed in light gray or not. Default: TRUE.}
  \item{legends}{
   If TRUE, legends with the meaning of colors are displayed. Default: FALSE.}
  \item{rel.tol}{
  Relative accuracy requested for the \code{integrate()} function from \pkg{stats} package. Default: 0.001.}
  \item{par.specify}{
  If FALSE, graphics are organized in one row and two columns. Default: FALSE.}
  \item{cex.lab, cex.axis, cex.main}{
  The magnification to be used for labels, axis annotation and main titles, respectively, relative to the current setting of \code{cex}. Default: \code{cex.lab = 1.5, cex.axis = 1.25, cex.main = 1.75}.}
  \item{lwd, col.roc}{
  Line width and color for the ROC curve. Default: \code{lwd = 2, col.roc = "blue"}.}
  \item{col}{
  A vector of length 2 with color used for controls and cases, respectively. Default: \code{c("#485C99", "#8F3D52")}.}
  \item{\dots}{
  Other parameters to be passed. Not used.}
}
\value{
A plot of the standard smooth ROC curve estimate with the selected graphical parameters
}
\examples{
data(HCC)
roc_cg20202438 <- gROC(X = HCC$cg20202438, D = HCC$tumor)
plot_densityROC(roc_cg20202438)
plot_densityROC(roc_cg20202438, h = c(2,2))
}

