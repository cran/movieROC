\name{predict}
\alias{predict}
\alias{predict.groc}
\alias{predict.hroc}

\title{
  Predict the classification regions for a particular specificity
}
\description{
This function prints the classification subsets corresponding to a particular false-positive rate \code{FPR} or to cutoff value(s) \code{C} or \code{XL, XU} introduced by the user.}
\usage{
\S3method{predict}{groc}(object, FPR = NULL, C = NULL, XL = NULL, XU = NULL, ...)
\S3method{predict}{hroc}(object, FPR = 0.15, ...)
}

\arguments{
  \item{object}{
  An object of class \sQuote{groc} or \sQuote{hroc}.}
  \item{FPR}{
  False-positive rate used to predict the classification region. Default: 0.15 if no cutoff value is provided by the next input parameters.}
  \item{C}{
  Cutoff value used to predict the classification region for \sQuote{groc} object with \code{side = "right"} or \code{"left"}. If \code{FPR} is provided, \code{C} is not used. Default: none.}
  \item{XL, XU}{
  Cutoff values used to predict the classification region for \sQuote{groc} object with \code{side = "both"} or \code{"both2"}. If \code{FPR} is provided, \code{C} is not used. Default: none.}
  \item{\dots}{
  Other parameters to be passed. Not used.}
}

\value{
  A list of length 3 with the following fields:
  \item{ClassSubsets}{A matrix with the classification region. Number of rows indicate the number of intervals whose union defines the classification region.}
  \item{Specificity}{Resulting specificity value.}
  \item{Sensitivity}{Resulting sensitivity value.}
}

\examples{
data(HCC)

roc <- gROC(X = HCC$cg18384097, D = HCC$tumor) # Right-sided ROC curve
predict(roc, FPR = 0.5)
groc <- gROC(X = HCC$cg18384097, D = HCC$tumor, side = "both") # gROC curve
predict(groc, FPR = 0.5)
hroc_cg18384097 <- hROC(X = HCC$cg18384097, D = HCC$tumor, 
    formula.lrm = "D ~ rcs(X,8)") 
predict(hroc_cg18384097, FPR = 0.5)
}
