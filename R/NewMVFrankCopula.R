#' Creates a Frank copula
#'
#' Creates an instance of the Frank copula with parameter \eqn{\alpha}.
#'
#' The following parameterisation of the copula is used:
#' \deqn{C(u_1,\dots,u_d) = -\log(1+\exp(s) * t)/\alpha}{C(u_1,\dots,u_d) = -log( 1 + exp(s) * t) / alpha}
#' where \eqn{s = \sum_{j=1}^d \log\left(\frac{\exp(-\alpha u_j) -1 }{t}\right)}{s = \sum_j log( (exp(-alpha u_j) -1) / t )} and \eqn{t=\exp(-\alpha)-1}{t = exp(-alpha) - 1}.
#'
#' @param alpha real, the parameter of the Frank copula, defaults to 2; must be positive.
#'
#' @return A function that evaluates the Frank copula (with parameter \eqn{\alpha}{alpha}) at a given \eqn{d}-dimensional vector in the unit cube.  The environment of the function also contains a function called \code{pdfCopula} that evaluates the probability density function of the Frank copula via automatic differentation.
#'
#' @author  Berwin A. Turlach <berwin.turlach@gmail.com>
#'
#' @keywords distribution
#'
#' @export
NewMVFrankCopula <- function(alpha = 1){
  if(alpha == 0)
    stop("Parameter alpha has to be positive.")

  param <- list(alpha = alpha)
  type <- c("Frank", "Multivariate")

  tmp <- exp(-alpha) - 1

  pdfCopula <- function(x){
    if( any(x<0) || any(x>1) )
      return(0)

    dim <- length(x)
    v <- NewCube(0, dim=dim)
    v$val[1] <- sum(log((exp(-alpha*x) - 1)/tmp))
    ii <- 2^(0:(dim-1))+1
    tt <- exp(-alpha*x)
    v$val[ii] <- -alpha*tt/(tt-1)
    res <- CrossExp(v)
    res$val <- res$val*tmp
    res$val[1] <- res$val[1] + 1
    res <- CrossLog(res)
    res$val <- -res$val/alpha

    if( !isTRUE(all.equal(res$val[1], cdfCopula(x), check.attributes=FALSE)) )
      stop("Inconsistency between implementation of cdf and pdf of this copula.")
    res$val[2^dim]
  }

  res <- cdfCopula <- function(x){
    sm <- sum(log((exp(-alpha*x) - 1)/tmp))
    -log(1 + exp(sm) * tmp )/alpha
  }
  class(res) <- "SimCop"
  res
}

