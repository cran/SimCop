#' Creates a Gumpel copula
#'
#' Creates an instance of the Gumbel copula with parameter r.  This family is also known as the Gumbel--Hougaard copula or the logistic model.
#'
#' The following parameterisation of the copula is used:
#' \deqn{C(u_1,\dots,u_d) = \exp\left(- \left\{ \sum_{j=1}^d \bar u_j^r \right\}^{1/r}\right)}{C(u_1,\dots,u_d) = exp(- { \sum_j v_j^r }^{1/r} )}
#' where \eqn{\bar u_j = -\log(u_j)}{v_j = -log(u_j)}, \eqn{j=1,\dots,d}.
#'
#' @param r real, the parameter of the Gumbel copula, defaults to 2, must be larger or equal to one.
#'
#' @return A function that evaluates the Gumbel copula (with parameter r) at a given \eqn{d}-dimensional vector in the unit cube.  The environment of the function also contains a function called \code{pdfCopula} that evaluates the probability density function of the Gumbel copula via automatic differentation.
#'
#' @seealso \code{\link{NewMEVAsyLogisticCopula}}
#'
#' @author  Berwin A. Turlach \email{berwin.turlach@@gmail.com}
#'
#' @keywords distribution
#'
#' @export
NewMEVGumbelCopula <- function(r = 2){

  if(r < 1)
    stop("Parameter r has to be larger or equal to one.")

  param <- list(r = r)
  type <- c("Gumbel", "Multivariate")

  pdfCopula <- function(x){
    if( any(x<0) || any(x>1) )
      return(0)

    dim <- length(x)
    v <- NewCube(0, dim=dim)
    v$val[1] <- sum((-log(x))^r)
    ii <- 2^(0:(dim-1))+1
    v$val[ii] <- r*(-log(x))^(r-1)*(-1/x)
    res <- CrossPow(v, 1/r)
    res$val <- -res$val
    res <- CrossExp(res)

    if( !isTRUE(all.equal(res$val[1], cdfCopula(x), check.attributes=FALSE)) )
      stop("Inconsistency between implementation of cdf and pdf of this copula.")
    res$val[2^dim]
  }

  res <- cdfCopula <- function(x){
    x <- pmin(pmax(x,0),1)
    exp(-sum((-log(x))^r)^(1/r))
  }
  class(res) <- "SimCop"
  res
}

