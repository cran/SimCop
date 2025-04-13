#' Creates a Clayton copula
#'
#' Creates an instance of the Clayton copula with parameter \eqn{\theta}{theta}.
#'
#' The following parameterisation of the copula is used:
#' \deqn{C(u_1,\dots,u_d) = \left(\left\{ \sum_{j=1}^d u_j^{-\theta} \right\} - (d-1) \right)^{-1/\theta}}{C(u_1,\dots,u_d) = ( { \sum_j u_j^{-theta} } - d + 1)^{-1/theta}}
#'
#' @param theta real, the parameter of the Clayton copula, defaults to 1; must be positive.
#'
#' @return A function that evaluates the Clayton copula (with parameter \eqn{\alpha}{alpha}) at a given \eqn{d}-dimensional vector in the unit cube.  The environment of the function also contains a function called \code{pdfCopula} that evaluates the probability density function of the Clayton copula via automatic differentation.
#'
#' @author  Berwin A. Turlach \email{berwin.turlach@@gmail.com}
#'
#' @keywords distribution
#'
#' @export
NewMVClaytonCopula <- function(theta = 1){
  if(theta <= 0)
    stop("Parameter theta has to be postive.")

  param <- list(theta = theta)
  type <- c("Clayton", "Multivariate")

  pdfCopula <- function(x){
    if( any(x<0) || any(x>1) )
      return(0)

    dim <- length(x)
    v <- NewCube(0, dim=dim)
    v$val[1] <- sum(x^(-theta)) - dim + 1
    ii <- 2^(0:(dim-1))+1
    v$val[ii] <- (-theta)*x^(-theta-1)
    res <- CrossPow(v, -1/theta)

    if( !isTRUE(all.equal(res$val[1], cdfCopula(x), check.attributes=FALSE)) )
      stop("Inconsistency between implementation of cdf and pdf of this copula.")
    res$val[2^dim]
  }

  res <- cdfCopula <- function(x){
    d <- length(x) - 1
    (sum(x^(-theta)) - d)^(-1/theta)
  }
  class(res) <- "SimCop"
  res
}
