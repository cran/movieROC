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
  This object can be \code{\link{print}}ed, \code{\link{plot}}ted, or \code{\link{predict}}ed
  for a particular point. It may be also passed to  
  \code{\link{plot_funregions}()} and \code{\link{plot_regions}()} functions.
}
\usage{
hROC(X, D, ...)
\S3method{hROC}{default}(X, D, type = c("lrm", "overfitting", "kernel", "h.fun"), 
    formula.lrm = "D ~ pol(X,3)", h.fun = function(x) {x},  kernel.h = 1,
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
    Type of transformation considered. One of 
    \code{"lrm"} (a binary logistic regression is computed by using \code{lrm()} function in \CRANpkg{rms} package), 
    \code{"kernel"} (the transformation included in Martínez-Camblor et al. (2021) estimated by the kernel density approach),
    \code{"overfitting"} (the overfitting transformation, \eqn{\hat{h}_{of}(\cdot)} is taken), or
    \code{"h.fun"} (the transformation indicated in the input parameter \code{h.fun} 
    is considered). Default: \code{"lrm"}.}
  \item{formula.lrm}{
    If \code{type = "lrm"}, the transformation employed in the right-hand side 
    of the logistic regression model (in terms of \code{X} and \code{D}). 
    Default: \code{'D ~ pol(X, 3)'}.}
  \item{kernel.h}{
    If \code{type = "kernel"}, the bandwidth used for the kernel density estimation by using the \code{density()} function in \pkg{stats} package. Default: 1.
    }
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
\details{
A theoretical and practical discussion about the \code{type} of transformation considered and its basis may be found in Martínez-Camblor et al. (2019) and Martínez-Camblor et al. (2021).

The overfitting function estimate is defined as follows:

\deqn{\hat{h}_{of}(x) = \sum_{i=1}^{n_1} I(x = y_i) + \sum_{i=1}^{n_2} \dfrac{\#(\xi = z_i)}{\#(\xi = z_i) + \#(\chi = z_i)} I(x = z_i)}

where \eqn{I(A)} denotes the indicator function (which takes the value 1 if \eqn{A} is true and 0 otherwise), \eqn{\#(B)} is the cardinal of the subset \eqn{B}, \eqn{\{y_1, \dots, y_{n_1}\} \subseteq \left\{ \xi_1, \dots, \xi_n \right\}} are the positive sample values without ties and \eqn{\{z_1, \dots, z_{n_2}\} \subseteq \left\{ \xi_1, \dots, \xi_n \right\}} are the positive sample values with ties with any negative sample value. Classification based on this transformation is the optimal one in the AUC sense, but the resulting decision rules cannot be extended to any other sample.
}
\value{
  A list of class \sQuote{hroc} with the following fields: 
  \item{levels}{Levels of response values.}
  \item{X, Y}{Original and transformed marker values, respectively.}
  \item{Sp, Se}{Vector of true-negtive and true-positive rates, respectively.}
  \item{auc}{Area under the curve estimate.}
  \item{model}{If \code{type = "lrm"}, the coefficients of the logistic regression model fitted by \code{formula}.}
}
\section{Dependencies}{
If \code{type = "lrm"}, the \code{lrm()} function in the \CRANpkg{rms} package is used. This library is also loaded to consider special transformation functions such as \code{pol()} and \code{rcs()}.
}
\references{
    P. Martínez-Camblor, S. Pérez-Fernández, and S. Díaz-Coto (2019) ``Improving the biomarker diagnostic capacity via functional transformations''. \emph{Journal of Applied Statistics}, \bold{46}(9): 1550–1566. DOI: \doi{10.1080/02664763.2018.1554628}.
    
    P. Martínez-Camblor, S. Pérez-Fernández, and S. Díaz-Coto (2021) ``Optimal classification scores based on multivariate marker transformations''. \emph{AStA Advances in Statistical Analysis}, \bold{105}(4): 581–599. DOI: \doi{10.1007/s10182-020-00388-z}.
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
##    by kernel density estimation with bandwidth 3
hROC(X, D, type = "kernel", kernel.h = 3)
}
