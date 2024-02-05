\name{gROC}
\alias{gROC}
\alias{gROC.default}
\title{
 Build a ROC curve for a univariate marker
}
\description{
  This is one of the main functions of the movieROC package.
  It builds a univariate ROC curve (standard or general) and returns a
  \sQuote{groc} object, a list of class \sQuote{groc}.
  This object can be \code{print}ed, \code{plot}ted, or \code{predict}ed
  for a particular point. It may be also passed to other functions:
  \code{\link{plot.regions}}, \code{\link{plot.buildROC}}, \code{\link{movieROC}},
  \code{\link{plot_densities}}, and \code{\link{plot_densityROC}}.
}
\usage{
gROC(X, D, ...)
\S3method{gROC}{default}(X, D, side = c("right", "left", "both", "both2"),
    restric = FALSE, optim = TRUE, t0 = NULL, t0max = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{X}{
    Vector of marker values.}
  \item{D}{
    Vector of response values. Two levels; if more, the two first ones are used.}
  \item{side}{
    Type of ROC curve. One of \code{"right"} (\eqn{\mathcal{R}_r(\cdot)}),
    \code{"left"} (\eqn{\mathcal{R}_l(\cdot)}), \code{"both"} (\eqn{\mathcal{R}_g(\cdot)})
    or \code{"both2"} (\eqn{\mathcal{R}_{g'}(\cdot)}). Default: \code{"right"}.}
  \item{restric}{
    If TRUE, the gROC curve with restriction \emph{(C)} is computed. Default: FALSE.}
  \item{optim}{
    If TRUE (and \code{restric = TRUE}), the computation of the optimal gROC curve
    under restriction \emph{(C)} is performed by using Floyd's algorithm. Default: TRUE.}
  \item{t0}{
    An integer number between 1 and \eqn{m+1} (where \eqn{m} is the negative sample size).
    If \code{restric = TRUE}, the restricted gROC curve is computed departing
    from (\code{t0}\eqn{-1})/\eqn{m}. Default: the one reporting the Youden index.}
  \item{t0max}{
    If TRUE, the computation of the gROC curve under restriction \emph{(C)} is performed
    departing from every possible \code{t0} and the one reporting the maximum AUC is selected.}
  \item{verbose}{
  If TRUE, a progress bar is displayed for computationally intensive methods. Default: FALSE.}
  \item{\dots}{
    Other parameters to be passed. Not used.}
}
\details{
  This function's main job is to estimate an ROC curve for a univariate marker
  under one of these considerations: larger values of the marker are associated
  with a higher probability of being positive (resulting in the \emph{right-sided}
  ROC curve, \eqn{\mathcal{R}_r (\cdot)}), the opposite (\emph{left-sided} ROC curve,
  \eqn{\mathcal{R}_l (\cdot)}), when both smaller and larger values of the marker are
  associated with having more probability of being positive (\emph{gROC curve},
  \eqn{\mathcal{R}_g(\cdot)}), the opposite (\emph{opposite gROC curve}, \eqn{\mathcal{R}_{g'}(\cdot)}).
}
\value{
  A list of class \sQuote{groc} with the following fields:
  \item{controls, cases}{Marker values of negative and positive subjects, respectively.}
  \item{levels}{Levels of response values.}
  \item{side}{Type of ROC curve.}
  \item{t}{Vector of false-positive rates.}
  \item{roc}{Vector of values of the ROC curve for \code{t}.}
  \item{c}{Vector of marker thresholds resulting in (\code{t}, \code{roc})
    if \code{side = "right" | "left"}.}
  \item{xl, xu}{Vectors of marker thresholds resulting in (\code{t}, \code{roc})
    if \code{side = "both" | "both2"}.}
  \item{auc}{Area under the curve estimate.}
  \item{aucfree}{Area under the curve estimate without restrictions.}
  \item{aucs}{Area under the curve with restriction \emph{(C)} departing from every
    false-positive rate, \eqn{FPR \in \{ 0, 1/m, \dots, 1 \}}.}
}
\examples{
data(HCC)

# ROC curve estimates for gene 03515901 and response tumor
gROC(X = HCC[,"cg03515901"], D = HCC$tumor) # Standard right-sided ROC curve
gROC(X = HCC[,"cg03515901"], D = HCC$tumor, side = "left") # Left-sided ROC curve
gROC(X = HCC[,"cg03515901"], D = HCC$tumor, side = "both") # gROC curve without restrictions

\donttest{
### Warning: Next line of code is time consuming. gROC curve with restriction (C)
gROC(X = HCC[,"cg03515901"], D = HCC$tumor, side = "both", restric = TRUE)
}
}
