\name{gROC_param}
\alias{gROC_param}
\title{
 Build a binormal ROC curve for a univariate marker
}
\description{
  This function builds a univariate ROC curve (standard or general) 
  assuming the binormal scenario with parameters being the sample estimates. 
  It returns a \sQuote{groc} object, a list of class \sQuote{groc}.  
}
\usage{
gROC_param(X, D, side = c("right", "left", "both", "both2"), N = NULL, ...)
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
  \item{N}{
    Number indicating the length of the vector of FPR considered to build the ROC curve: 
    \eqn{t \in \{ 0, 1/N, 2/N, \dots, 1 \}}. Default: \code{1000}.}
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
  \item{a, b}{Estimates for parameters \eqn{a} and \eqn{b} considered for the ROC curve estimation: 
    \eqn{\hat{a} = \left[ \overline{\xi_n} - \overline{\chi_m} \right]/\hat{s}_\xi} and 
    \eqn{\hat{b} = \hat{s}_\chi / \hat{s}_\xi}.}
  \item{p0}{Estimate of the "central value", \eqn{\mu^*}, about to which the thresholds \eqn{x^L} 
    and \eqn{x^U} are symmetrical.}
}
\references{
  P. Martínez-Camblor and J. C. Pardo-Fernández (2019) ``Parametric estimates for the receiver operating characteristic curve generalization for non-monotone relationships''. \emph{Statistical Methods in Medical Research}, \bold{28}(7): 2032–2048. DOI: \doi{10.1177/0962280217747009}.
}
\examples{
data(HCC)

# ROC curve estimates for gene 03515901 and response tumor assuming the binormal scenario
gROC_param(X = HCC[,"cg03515901"], D = HCC$tumor) # Standard right-sided ROC curve
gROC_param(X = HCC[,"cg03515901"], D = HCC$tumor, side = "left") # Left-sided ROC curve 
gROC_param(X = HCC[,"cg03515901"], D = HCC$tumor, side = "both") # gROC curve 
}
