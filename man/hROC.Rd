\name{hROC}
\alias{hROC}
\alias{hROC.default}
\title{
 Build a ROC curve for a transformation of a univariate marker
}
\description{
  This is one of the main functions of the movieROC package. 
  It builds a univariate ROC curve for a transformed marker \eqn{h(X)} 
  and returns a \sQuote{hroc} object, a list of class \sQuote{hroc}. 
  This object can be \code{print}ed, \code{plot}ted, or \code{predict}ed
  for a particular point. It may be also passed to  
  \code{\link{plot.funregions}} and \code{\link{plot.regions}} functions.
}
\usage{
hROC(X, D, ...)
\S3method{hROC}{default}(X, D, type = c("lrm", "h.fun", "overfitting"), 
    formula.lrm = "D ~ pol(X,3)", h.fun = function(x) {x}, 
    plot.h = FALSE, plot.roc = FALSE, new.window = FALSE, 
    main = NULL, xlab = "x", ylab = "h(x)", xaxis = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    Vector of marker values.}
  \item{D}{
    Vector of response values. Two levels; if more, the two first ones are used.}
  \item{type}{
    Type of transformation considered. One of \code{"lrm"} (a binary logistic 
    regression is computed by using \code{lrm} function in \code{rms} package), 
    \code{"h.fun"} (the transformation indicated in the input parameter \code{h.fun} 
    is considered) or \code{"overfitting"} (the overfitting transformation, 
    \eqn{\ell_{of}(\cdot)} is taken). Default: \code{"lrm"}.}
  \item{formula.lrm}{
    If \code{type = "lrm"}, the transformation employed in the right-hand side 
    of the logistic regression model (in terms of \code{X} and \code{D}). 
    Default: \code{'D ~ pol(X, 3)'}.}
  \item{h.fun}{
      If \code{type = "h.fun"}, the transformation employed (as a function in R). 
      Default: \code{function(x){x}}.}
  \item{plot.h}{
    If TRUE, the transformation employed is illustrated.}
  \item{plot.roc}{
  If TRUE, the resulting ROC curve is illustrated.}
  \item{new.window}{
  If TRUE, two previous graphics are plotted separately in different windows.
}
  \item{main}{
  A main title for the plot used if \code{plot.h = TRUE}. 
}
  \item{xlab, ylab}{
  A label for the x and y axis of the plot used if \code{plot.h = TRUE}. 
}
  \item{xaxis}{
  Graphical parameter used if \code{plot.h = TRUE}. If FALSE, plotting of the axis is supressed.
}
  \item{\dots}{
    Other parameters to be passed. Not used.}
}
\value{
  A list of class \sQuote{hroc} with the following fields: 
  \item{levels}{Levels of response values.}
  \item{X, Y}{Original and transformed marker values, respectively.}
  \item{Sp, Se}{Vector of true-negtive and true-positive rates, respectively.}
  \item{auc}{Area under the curve estimate.}
  \item{model}{If \code{type = "lrm"}, the coefficients of the logistic regression model fitted by \code{formula}.}
}
\examples{
data(HCC)

# ROC curve for gene 18384097 to identify tumor by considering 4  different transformations:
X <- HCC$cg18384097; D <- HCC$tumor
## 1. Ordinary cubic polynomial formula for binary logistic regression
hROC(X, D)
## 2. Linear tail-restricted cubic splines for binary logistic regression
hROC(X, D, formula.lrm = "D ~ rcs(X,8)")
## 3. Overfitting transformation for this particular sample
hROC(X, D, type = "overfitting")
## 4. Optimal transformation in terms of likelihood ratio 
##    by kernel density estimation with bandwidth h=3
EstDensTransf_FUN <- function(X, D, h = 1){
  D <- as.factor(D)
  controls <- X[D == levels(D)[1]]; dens_controls <- density(controls, adjust = h)
  cases <- X[D == levels(D)[2]]; dens_cases <- density(cases, adjust = h)
  dens_controls_FUN <- approxfun(dens_controls$x, dens_controls$y, rule = 0)
  dens_cases_FUN <- approxfun(dens_cases$x, dens_cases$y, rule = 0)
  function(x) dens_cases_FUN(x)/(dens_controls_FUN(x) + dens_cases_FUN(x))
}
hROC(X, D, type = "h.fun", h.fun = function(x) EstDensTransf_FUN(X, D, h = 3)(x))
}
