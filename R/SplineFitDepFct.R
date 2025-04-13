#' Fit a dependence function by spline smoothing
#'
#' Given estimates for the dependence function of a bivariate extreme value copula at specified points, this function fits a natural cubic smoothing spline, that is constrained to fulfill all the conditions of a dependence function, to the given estimates via quadratic programming.
#'
#' \code{integ.value} should be between 0 and 2.  If a value is specified, then an additional constraint is added to the quadratic program to ensure that the integeral (over 0 to 1) of the second derivative of the spline is larger or equal to \code{integ.value}.  Choosing values close to 2 may lead to quadratic programms on which \code{\link[quadprog]{solve.QP}} reports inconsistent constraints.
#'
#' @param x,y vectors giving the coordinates of the points to be approximated. Alternatively a single plotting structure can be specified: see \code{\link{xy.coords}}.
#' @param alpha real, the smoothing parameter for the smoothing splines.
#' @param integ.value real, non-negative value that should be less than two; see \sQuote{Details}
#'
#' @return A function, created by \code{\link{splinefun}}, that evaluates the natural cubic spline that was fitted to the data.
#'
#' @examples
#' ## Data from Hall and Tajvidi (2004, ANZJS)
#' EstDF2 <- NonparEstDepFct(MaxTemp, convex = FALSE)
#'
#' ## Plot modified Pickands Function and area in which
#' ## dependence function must lie
#' plot(EstDF2, ylim = c(0.5,1), xlab = "w", ylab = "A(w)", type="l", lty="longdash")
#' polygon(c(0, 0.5, 1, 0), c(1, 0.5, 1, 1))
#'
#' ## Fit spline to Pickands function and add to plot
#' splfit <- SplineFitDepFct(EstDF2)
#' curve(splfit, n = 301, add = TRUE, lty = "dashed")
#'
#' @references Hall, P. and Tajvidi, N. (2000). Distribution and dependence-function estimation for bivariate extreme-value distributions.  \emph{Bernoulli} \bold{6}(5): 835--844. \doi{10.2307/3318758}.
#'
#' Hall, P. and Tajvidi, N. (2004). Prediction regions for bivariate extreme events. \emph{Australian & New Zealand Journal of Statistics} \bold{46}(1): 99--112. \doi{10.1111/j.1467-842X.2004.00316.x}.
#'
#' @author Nader Tajvidi \email{Nader.Tajvidi@@matstat.lu.se}
#' @author Berwin A Turlach \email{Berwin.Turlach@@gmail.com}
#'
#' @seealso \code{\link{NonparEstDepFct}}, \code{\link{NewBEVSplineCopula}}
#'
#' @keywords smooth
#'
#' @export
SplineFitDepFct <- function(x, y = NULL, alpha = 0.01, integ.value)
{
  tmp <- grDevices::xy.coords(x, y, setLab=FALSE)
  w <- tmp$x
  y <- tmp$y


  n.x <- length(w)
  n.y <- length(y)
  if(n.x != n.y)
    stop("observations don't have the same length\n")

  n <- length(w)

  hvec <- w[2:n] - w[1:(n - 1)]
  Rmat <- matrix(0, n - 2, n - 2)
  diag(Rmat) <- (hvec[1:(n - 2)] + hvec[2:(n - 1)])/3
  diag(Rmat[-1,  ]) <- hvec[2:(n - 2)]/6
  diag(Rmat[, -1]) <- hvec[2:(n - 2)]/6
  Qmat <- matrix(0, n, n - 2)
  diag(Qmat) <- 1/hvec[1:(n - 2)]
  diag(Qmat[-1,  ]) <- -1/hvec[1:(n - 2)] - 1/hvec[2:(n - 1)]
  diag(Qmat[ - (1:2),  ]) <- 1/hvec[2:(n - 1)]

  Amat <- rbind(Qmat,  - t(Rmat))
  meq <- ncol(Amat)
  bvec <- rep(0, meq)

  meq <- meq + 2
  bvec <- c(1, 1, bvec)
  Amat <- cbind(c(1, rep(0, 2 * n - 3)),
                c(rep(0, n - 1), 1, rep(0, n - 2)),
                Amat)

  dvec <- c(y, rep(0, n - 2))
  Dmat <- matrix(0, n, n)
  diag(Dmat) <- 1
  Dmat <- rbind(cbind(Dmat,
                      matrix(0, n, n - 2)),
                cbind(matrix(0,n - 2, n),
                      alpha * Rmat))
  nadd <- n - 2
  Tmat <- matrix(0, n - 2, nadd)
  diag(Tmat) <- 1
  Amat <- cbind(Amat, rbind(matrix(0, n, nadd), Tmat))
  bvec <- c(bvec,rep(0, nadd))

  if(!missing(integ.value)) {
    cat("\n   Calculating with constrained integral of second derivative larger or equal to",
        integ.value, ".\n")
    delt <- 1/n
    bvec <- c(bvec, integ.value/delt)
    Amat <- cbind(Amat, c(rep(0, n), rep(1, n - 2)))
  }

  h <- w[n] - w[n - 1]
  newcol <- c(rep(0, n - 2), -1/h, 1/h, rep(0, n - 3), h/6)
  Amat <- cbind(Amat, newcol)
  bvec <- c(bvec, 0)
  Amat <- cbind(Amat,  - newcol)
  bvec <- c(bvec, -1)
  h <- w[2] - w[1]
  newcol <- c(-1/h, 1/h, rep(0, n - 2),  - h/6, rep(0, n - 3))
  Amat <- cbind(Amat,  - newcol)
  bvec <- c(bvec, 0)
  Amat <- cbind(Amat, newcol)
  bvec <- c(bvec, -1)

  res2 <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)

  gvi <- res2$solution[1:n]
  gami <- c(0, res2$solution[(n+1):length(res2$solution)], 0)

  ac <- gvi
  bc <- c( (gvi[2:n]-gvi[1:(n-1)])/hvec - hvec/6*(2*gami[1:(n-1)]+gami[2:n]),
          (gvi[n]-gvi[n-1])/hvec[n-1] + hvec[n-1]/6*gami[n-1]  )
  cc <- gami/2
  dc <- c( (gami[2:n]-gami[1:(n-1)])/(6*hvec), 0)

  tt <- stats::splinefun(w, ac, method = "natural")

  zz <- get("z", envir=environment(tt))
  if( !isTRUE(all.equal(zz$y, ac, check.attributes=FALSE)) ){
    stop("Inconsistency between values of splines.")
  }
  if( !isTRUE(all.equal(zz$b, bc, check.attributes=FALSE)) ){
    stop("Inconsistency between first derivatives of splines.")
  }
  if( !isTRUE(all.equal(zz$c, cc, check.attributes=FALSE)) ){
    stop("Inconsistency between second derivatives of splines.")
  }
  if( !isTRUE(all.equal(zz$d, dc, check.attributes=FALSE)) ){
    stop("Inconsistency between third derivatives of splines.")
  }

  tt
}

