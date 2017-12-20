#' Creates a bivariate mixed model extreme value copula
#'
#' Creates an instance of the bivariate asymmetric mixed model extreme value copula with parameter \eqn{\theta}.
#'
#' The dependence function for this bivariate EV copula is \deqn{A(w) = \theta w^2 - \theta w + 1}
#' Necessary and sufficient conditions for the dependence function to be valid are
#' \itemize{
#'   \item \eqn{0 \le \theta \le 1}
#' }
#'
#' @param theta real.
#'
#' @return A function that evaluates the bivariate asymmetric mixed model EV copula (with parameter \eqn{\theta}) at a given \eqn{2}-dimensional vector in the unit square.  The environment of the function also contains a function called \code{pdfCopula} that evaluates the probability density function of the bivariate asymmetric mixed model EV copula via automatic differentation.
#'
#' @seealso \code{\link{NewBEVAsyMixedModelCopula}}
#'
#' @author  Berwin A. Turlach <berwin.turlach@gmail.com>
#'
#' @keywords distribution
#'
#' @export
NewBEVMixedModelCopula <- function(theta){

  if( theta < 0 || theta > 1)
    stop("'theta' must be between 0 and 1 (inclusive).")

  param <- list(theta = theta)
  type <- c("mixed model extreme value", "Bivariate")

  pdfCopula <- function(x){
    if( length(x) != 2 )
      stop("x has the wrong length, should be 2.")
    if( any(x<0) || any(x>1) )
      return(0)

    slx <- NewCube(0, dim=2)
    slx$val[1] <- sum(log(x))
    slx$val[2:3] <- 1/x

    lx1 <- NewCube(0, dim=2)
    lx1$val[1] <- log(x[1])
    lx1$val[2] <- 1/x[1]
    lx1 <- CrossDivide(lx1, slx)

    res <- NewCube(0, dim=2)
    res$val <- theta*lx1$val
    res$val[1] <- res$val[1] - theta
    res <- CrossMult(res, lx1)
    res$val[1] <- res$val[1] + 1

    res <- CrossExp(CrossMult(slx,res))
    if( !isTRUE(all.equal(res$val[1], cdfCopula(x), check.attributes=FALSE)) ){
      stop("Inconsistency between implementation of cdf and pdf of this copula.")
    }
    res$val[4]
  }

  res <- cdfCopula <- function(x){
    if( length(x) != 2 )
      stop("x has the wrong length, should be 2.")

    if( any(x==0) ) return(0)
    if( all(x==1) ) return(1)

    tmp <- log(x)
    sumtmp <- sum(tmp)

    tmp1 <- tmp[1]/sumtmp
    tmp1 <- (theta*tmp1 - theta)*tmp1 + 1

    exp(sumtmp*tmp1)
  }
  class(res) <- "SimCop"
  res
}
