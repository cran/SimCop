#' SimCop: A package to simulate random variates from an arbitrary multivariate copula
#'
#' @description
#' R code to support Tajvidi and Turlach (2017).  The main functions implemented for the SimCop package are:
#' \itemize{
#' \item\code{New*Copula}, various functions that create objects of \code{\link{class}} \sQuote{SimCop}. These functions return a copula function with various helpful information stored in the environment of the function.
#'
#'     Details of the implementation are subject to change and should not be relied on.
#'
#'     Only a  \code{\link{print}} method is implemented for this class so far.
#'
#' \item\code{GetApprox}, a function that calculates approximations to a copula and returns an object of \code{\link{class}} \sQuote{CopApprox}.
#'
#'     For bivariate copulae a method for \code{\link{plot}} is implemented for this class.
#'
#' \item\code{GenerateRV}, a generic function that generates random variates from an object, together with a method for objects of class \sQuote{CopApprox}.
#' }
#'
#' @references Tajvidi, N. and Turlach, B.A. (2017). A general approach to generate random variates for multivariate copulae, \emph{Australian & New Zealand Journal of Statistics} \bold{60}(1): 140--155. \doi{10.1111/anzs.12209}.
#'
#' @name SimCop
"_PACKAGE"
