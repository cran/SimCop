#' Creates a bivariate asymmetric logistic model extreme value copula
#'
#' Creates an instance of the bivariate asymmetric logistic model extreme value copula with parameters \eqn{r}, \eqn{\theta} and \eqn{\phi}.
#'
#' The dependence function for this bivariate EV copula is \deqn{A(w) = (\theta (1-w)^r + (\phi w)^r)^{1/r} + (\theta - \phi) w + 1 - \theta}
#' Necessary and sufficient conditions for the dependence function to be valid are
#' \itemize{
#'   \item \eqn{r \ge 1}
#'   \item \eqn{0 \le \theta \le 0}
#'   \item \eqn{0 \le \phi \le 1}
#' }
#' For \eqn{\theta = \phi = 1} this model reduces to the symmetric logistic model.
#'
#'
#' @param r real.
#' @param theta real.
#' @param phi  real.
#'
#' @return A function that evaluates the bivariate asymmetric logistic model EV copula (with parameters \eqn{r}, \eqn{\theta} and \eqn{phi}) at a given \eqn{2}-dimensional vector in the unit square.  The environment of the function also contains a function called \code{pdfCopula} that evaluates the probability density function of the bivariate asymmetric mixed model EV copula via automatic differentation.
#'
#' @seealso \code{\link{NewBEVLogisticCopula}}, \code{\link{NewMEVAsyLogisticCopula}}
#'
#' @author  Berwin A. Turlach <berwin.turlach@gmail.com>
#'
#' @keywords distribution
#'
#' @export
NewBEVAsyLogisticCopula <- function(r, theta, phi){

  if( r <  1 )
    stop("'r' must be equal or larger than 1.")
  if( theta <  0 )
    stop("'theta' must be equal or larger than 0.")
  if( phi < 0 || phi >  1 )
    stop("'phi' must be between 0 and 1 (inclusive).")

  param <- list(r = r)
  type <- c("logistic model extreme value", "Bivariate")

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
    res$val <- -lx1$val
    res$val[1] <- res$val[1] + 1
    res$val <- theta * res$val
    tmp <- NewCube(0, dim=2)
    tmp$val <- phi * lx1$val
    res <- CrossPow(CrossSum(CrossPow(res, r), CrossPow(tmp, r)), 1/r)
    res$val <- res$val + (theta-phi)*lx1$val
    res$val[1] <- res$val[1] + 1 - theta

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
    tmp1 <- ((theta*(1-tmp1))^r + (phi*tmp1)^r)^(1/r) + (theta-phi)*tmp1 + 1 - theta

    exp(sumtmp*tmp1)
  }
  class(res) <- "SimCop"
  res
}
