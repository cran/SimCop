% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NewBEVSplineCopula.R
\name{NewBEVSplineCopula}
\alias{NewBEVSplineCopula}
\title{Creates a flexible extreme value copula}
\usage{
NewBEVSplineCopula(spl)
}
\arguments{
\item{spl}{a spline function.}
}
\value{
A function that evaluates the bivariate EV copula (whose dependence function is given by the spline) at a given \eqn{2}-dimensional vector in the unit square.  The environment of the function also contains a function called \code{pdfCopula} that evaluates the probability density function of the bivariate asymmetric mixed model EV copula via automatic differentation.
}
\description{
Creates a bivariate extreme value copula from a spline estimate of its dependence function.
}
\seealso{
\code{\link{SplineFitDepFct}}
}
\author{
Berwin A. Turlach \email{berwin.turlach@gmail.com}
}
\keyword{distribution}
