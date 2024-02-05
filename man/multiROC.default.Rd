\name{multiROC}
\alias{multiROC}
\alias{multiROC.default}
\title{
 Build a ROC curve for a multivariate marker with dimension higher than two
}
\description{
  This is one of the main functions of the movieROC package.
  It builds a multivariate ROC curve by considering one of these methods:
  i) fitting a binary logistic regression model with a particular combination
  (fixed by the user) of the two components on the right-hand side,
  ii) linear combinations with fixed parameters, or
  iii) linear combinations with dynamic parameters.
  It returns a \sQuote{multiroc} object, a list of class \sQuote{multiroc}.
  This object can be \code{print}ed or \code{plot}ted. It may be also passed to
  \code{\link{plot.buildROC}} and \code{\link{movieROC}} function.
}
\usage{
multiROC(X, D, ...)
\S3method{multiROC}{default}(X, D, method = c("lrm", "fixedLinear",
        "dynamicMeisner", "kernelOptimal"),
    formula.lrm = "D ~ X.1 + I(X.1^2) + X.2 + I(X.2^2) + I(X.1*X.2)",
    stepModel = TRUE, coefLinear = rep(1, ncol(X)),
    methodLinear = c("coefLinear", "SuLiu", "PepeThompson", "logistic", "minmax"),
    alpha = 0.5, approxh = 0.5, multiplier = 2, kernelOptimal.H =
        c("Hpi", "Hscv", "Hns", "Hlscv", "Hscv.diag", "Hpi.diag", "Hlscv.diag"),
    eps = sqrt(.Machine$double.eps), verbose = FALSE, ...)
}
\arguments{
  \item{X}{
  Matrix (dimension \eqn{n \times p}) of marker values where \eqn{n} is the sample size and \eqn{p} is the dimension of the multivariate marker.}
  \item{D}{
  Vector of response values. Two levels; if more, the two first ones are used.}
  \item{method}{
  Method used to build the classification regions. One of \code{"lrm"} (fitting a binary logistic regression model by the input parameter \code{formula}), \code{"fixedLinear"} (linear frontiers with fixed parameters given in \code{coefLinear}), \code{"dynamicMeisner"} (linear frontiers with dynamic parameters reported by Meisner et al. (2021) method), or \code{"kernelOptimal"} (estimating optimal transformation based on bivariate kernel density estimation by Martínez-Camblor et al. (2021) using the \code{kde} function in the \CRANpkg{ks} package). Default: \code{"lrm"}.}
  \item{formula.lrm}{
  If \code{method = "lrm"}, the transformation employed in the right-hand side of the logistic regression model (in terms of \code{X.1}, \code{X.2} \code{dots}, \code{X.p}, and \code{D}). Default: quadratic formula for the two first components \code{X.1} and \code{X.2}.}
  \item{stepModel}{
  If TRUE and \code{method = "lrm"}, a model selection is performed based on the AIC (Akaike information criterion) in a stepwise algorithm (see \code{step} function for more information). Default: TRUE.}
  \item{coefLinear}{
  If \code{method = "fixedLinear"} and \code{methodLinear = "coefLinear"}, a vector of length \eqn{p} with the coefficients \eqn{\beta_i} (\eqn{i \in  \{1, \dots, p\}}) used to \eqn{\mathcal{L}_{\boldsymbol{\beta}}(\boldsymbol{X}) = \beta_1 X_1 + \dots + \beta_p X_p}. Default: \eqn{(1,\dots,1)}.}
  \item{methodLinear}{
  If \code{method = "fixedLinear"}, method used to build the classification regions. One of \code{"coefLinear"} (particular fixed coefficients in \code{coefLinear}), \code{"SuLiu"} (Su and Liu, 1993), \code{"PepeThompson"} (Pepe and Thompson, 2000), \code{"logistic"} (logistic regression model), \code{"minmax"} (Liu et al., 2011). Default: \code{"coefLinear"}.}
  \item{alpha, approxh, multiplier}{
  If \code{method = "dynamicMeisner"}, input parameters used in the \code{maxTPR} function of the \CRANpkg{maxTPR} package. Default: \code{alpha = 0.5}, \code{approxh = 0.5} and \code{multiplier = 2}.}
  \item{kernelOptimal.H}{
  If \code{method = "kernelOptimal"}, the bandwidth matrix \code{H} used in the \code{kde} function of the \CRANpkg{ks} package. Default: \code{"Hpi"} (plug-in bandwidth selector).}
  \item{eps}{
  Epsilon value to consider. Default: \code{sqrt(.Machine$double.eps)}.}
  \item{verbose}{
  If TRUE, a progress bar is displayed for computationally intensive methods. Default: FALSE.}
  \item{\dots}{
  Other parameters to be passed. Not used.}
}
\value{
  A list of class \sQuote{multiroc} with the following fields:
  \item{controls, cases}{Marker values of negative and positive subjects, respectively.}
  \item{levels}{Levels of response values.}
  \item{t}{Vector of false-positive rates.}
  \item{roc}{Vector of values of the ROC curve for \code{t}.}
  \item{auc}{Area under the curve estimate.}
  \item{Z}{If \code{method = "lrm"}, \code{"fixedLinear"}, or \code{"kernelOptimal"}, resulting
univariate marker values.}
  \item{c}{If \code{method = "lrm"}, \code{"fixedLinear"}, or \code{"kernelOptimal"}, vector of final
marker thresholds resulting in (\code{t}, \code{roc}).}
  \item{CoefTable}{If \code{method = "dynamicMeisner"} or \code{"dynamicEmpirical"}, a list of length equal to length of vector \code{t}. Each element of the list keeps the linear coefficients (\code{coef}), threshold for such linear combination (\code{c}), the corresponding point in the ROC curve (\code{t}, \code{roc}), the resulting univariate marker values (\code{Z}) and a matrix of dimension 100 \eqn{\times} 100 with the marker values over a grid of (\eqn{X_1}, \eqn{X_2}) bivariate values (\code{f}).}
}
\references{
  J. Q. Su and J. S. Liu. (1993) ``Linear combinations of multiple diagnostic markers''. \emph{Journal of the American Statistical Association}, \bold{88}(424): 1350–1355. DOI: \doi{10.1080/01621459.1993.10476417}.

  M. S. Pepe and M. L. Thompson (2000) ``Combining diagnostic test results to increase accuracy''. \emph{Biostatistics}, \bold{1}
(2):123–140. DOI: \doi{10.1093/biostatistics/1.2.123}.

  C. Liu, A. Liu, and S. Halabi (2011) ``A min–max combination of biomarkers to improve diagnostic accuracy''. \emph{Statistics in Medicine}, \bold{30}(16): 2005–2014. DOI: \doi{10.1002/sim.4238}.

    P. Martínez-Camblor, S. Pérez-Fernández, and S. Díaz-Coto (2021) ``Optimal classification scores based on multivariate marker transformations''. \emph{AStA Advances in Statistical Analysis}, \bold{105}(4): 581–599. DOI: \doi{10.1007/s10182-020-00388-z}.

    A. Meisner, M. Carone, M. S. Pepe, and K. F. Kerr (2021) ``Combining biomarkers by maximizing the true positive rate for a fixed false positive rate''. \emph{Biometrical Journal}, 63(6): 1223–1240. DOI: \doi{10.1002/bimj.202000210}.
}

\examples{
data(HCC)

# ROC curve for genes 20202438, 18384097, and 03515901 to identify tumor
# by 4 different methods:
X <- cbind(HCC$cg20202438, HCC$cg18384097, HCC$cg03515901); D <- HCC$tumor
## 1. Linear combinations with fixed parameters by Pepe and Thompson (2000)
multiROC(X, D, method = "fixedLinear", methodLinear = "PepeThompson")
## 2.Linear combinations with dynamic parameters by Meisner et al. (2021)
\donttest{
### Time consuming
multiROC(X, D, method = "dynamicMeisner")}
## 3. Logistic regression model with quadratic formula by default
multiROC(X, D)
## 4. Optimal transformation with multivariate KDE by Martínez-Camblor et al. (2021)
multiROC(X, D, method = "kernelOptimal")
}

