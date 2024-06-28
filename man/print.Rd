\name{print}
\alias{print.groc}
\alias{print.hroc}
\alias{print.multiroc}

\title{
  Print an ROC curve object
}
\description{
  This function prints information about an ROC curve.
}

\usage{
\S3method{print}{groc}(x, ...)
\S3method{print}{hroc}(x, ...)
\S3method{print}{multiroc}(x, ...)
}
\arguments{
  \item{x}{
An ROC curve object from the \pkg{movieROC} package. Possible classes are: 
\sQuote{groc} (output of \code{\link{gROC}()} function), \sQuote{hroc} (output of \code{\link{hROC}()} function), and \sQuote{multiroc} (output of \code{\link{multiROC}()} function).}
  \item{\dots}{
  Other parameters to be passed. Not used.
  }
}
\value{
A character vector with the details of an ROC curve object.
}
\examples{
## See examples for gROC(), hROC(), and multiROC() function
}
