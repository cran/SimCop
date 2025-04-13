#' Creates a flexible extreme value copula
#'
#' Creates a bivariate extreme value copula from a spline estimate of its dependence function.
#'
#' @param spl a spline function.
#'
#' @return A function that evaluates the bivariate EV copula (whose dependence function is given by the spline) at a given \eqn{2}-dimensional vector in the unit square.  The environment of the function also contains a function called \code{pdfCopula} that evaluates the probability density function of the bivariate asymmetric mixed model EV copula via automatic differentation.
#'
#' @author  Berwin A. Turlach \email{berwin.turlach@@gmail.com}
#'
#' @seealso \code{\link{SplineFitDepFct}}
#'
#' @keywords distribution
#'
#' @export
NewBEVSplineCopula <- function(spl){

  param <- list(spl = "Spline function")
  type <- c("extreme value", "Bivariate")

  z <- get("z", envir=environment(spl))
  lzxm1 <- length(z$x) - 1
  ii <- 1:lzxm1
  iip1 <- ii + 1

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

    x0 <- lx1$val[1]
    if(is.na(x0))
      return(NaN)
    jj <- min(which(x0 >= z$x[ii] & x0 < z$x[iip1]), lzxm1)

    lx1$val[1] <- lx1$val[1] - z$x[jj]
    res <- NewCube(0, dim=2)
    res$val <- z$d[jj] * lx1$val
    res$val[1] <- res$val[1] + z$c[jj]
    res <- CrossMult(res, lx1)
    res$val[1] <- res$val[1] + z$b[jj]
    res <- CrossMult(res, lx1)
    res$val[1] <- res$val[1] + z$y[jj]

    res <- CrossExp(CrossMult(slx, res))
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

    tmp1 <- spl(tmp[1]/sumtmp)

    exp(sumtmp*tmp1)
  }
  class(res) <- "SimCop"
  res
}
