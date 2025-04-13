#' Nonparametric estimator of bivariate dependence function
#'
#' Function to calculate nonparametric estimates of the dependence functions of bivariate extreme value copula.
#'
#' If \code{transf.to.frechet} is \code{TRUE}, the default, then a generalised extreme value (GEV) distribution is fitted to each margin and the fitted parameters are used to transform the data to have standard \enc{Fréchet}{Frechet} margins.  The parameterisation of the cumulative distribution of the GEV that is used is, if \eqn{\gamma \neq 0}{\gamma <> 0}:
#' \deqn{G(z) = \exp\left(-\left[1+\gamma\left(\frac{z-\mu}{\sigma}\right)\right]^{-1/\gamma}\right)}{G(z) = exp( -{1 + \gamma [(z - \mu)/\sigma] }^(-1/\gamma) )}
#' and for \eqn{\gamma = 0}:
#' \deqn{G(z) = \exp(-\exp(-z))}{G(z) = exp( -exp(-z) )}
#' If \eqn{\gamma < 0}, then the support of the GEV is the interval \eqn{(-\infty, \mu - \sigma/\gamma]}, while it is \eqn{[\mu - \sigma/\gamma, \infty)} if \eqn{\gamma > 0}.  For \eqn{\gamma = 0}, the support is the real line.
#'
#' If \code{verbose} is \code{TRUE}, not the default, and \code{transf.to.frechet} is \code{TRUE}, the estimates for the fitted GEV distribution are printed out using \code{\link{cat}}.
#'
#' @param x,y vectors giving the observations of the extreme values. Alternatively a single plotting structure can be specified: see \code{\link{xy.coords}}.
#' @param w.length number of grid points (using an equidistant grid from 0 to 1) on which the dependence function is estimated.
#' @param transf.to.frechet logical, controls whether \code{x} and \code{y} are first transformed to have standard \enc{Fréchet}{Frechet} margins: see \sQuote{Details}; defaults to \code{TRUE}.
#' @param convex.hull logical, controls whether the convex hull of the modified Pickands estimator is returned; defaults to \code{TRUE}.
#' @param verbose logical, controls whether progress messages are given; defaults to \code{FALSE}.
#'
#' @return A list with two named components.  The component \code{x} contains a vector with the grid points at which the dependence function was estimated.  The component \code{y} contains the estimated dependence functions.
#'
#' @examples
#' ## Data from Hall and Tajvidi (2004, ANZJS)
#' EstDF1 <- NonparEstDepFct(MaxTemp)
#'
#' ## Plot modified Pickands Function and area in which
#' ## dependence function must lie
#' plot(EstDF1, ylim = c(0.5,1), xlab = "w", ylab = "A(w)", type="l", lty="longdash")
#' polygon(c(0, 0.5, 1, 0), c(1, 0.5, 1, 1))
#'
#' @seealso \code{\link{SplineFitDepFct}}
#'
#' @references Hall, P. and Tajvidi, N. (2000). Distribution and dependence-function estimation for bivariate extreme-value distributions.  \emph{Bernoulli} \bold{6}(5): 835--844. \doi{10.2307/3318758}.
#'
#' Hall, P. and Tajvidi, N. (2004). Prediction regions for bivariate extreme events. \emph{Australian & New Zealand Journal of Statistics} \bold{46}(1): 99--112. \doi{10.1111/j.1467-842X.2004.00316.x}.
#'
#' @keywords nonparametric
#'
#' @author Nader Tajvidi \email{Nader.Tajvidi@@matstat.lu.se}
#'
#' @export
NonparEstDepFct <- function(x , y = NULL, w.length = 101,
                            transf.to.frechet = TRUE,
                            convex.hull = TRUE, verbose = FALSE)
{
  tmp <- grDevices::xy.coords(x, y, setLab=FALSE)
  obs.x <- tmp$x
  obs.y <- tmp$y


  if(transf.to.frechet) {
    if(verbose){
      cat("\n\t Transforming to Frechet marginals ... \n")
      cat("\t\tEstimates of the parameters in GEV for the first margin are:\n")
    }
    x <- transformToFrechet(obs.x,verbose=verbose)
    if(verbose){
      cat("\t\tEstimates of the parameters in GEV for the second margin are:\n")
    }
    y <- transformToFrechet(obs.y,verbose=verbose)
    non.valid.index <- (is.na(x)) | (is.na(y))
    if(any(non.valid.index)) {
      cat("some non valid observations will be removed\n")
      x <- x[!non.valid.index]
      y <- y[!non.valid.index]
    }
  }
  else {
    if(verbose)
      cat("We use the original observations without transforming to Frechet marginals ... \n")
    x <- obs.x
    y <- obs.y
  }
  w <- seq(0, 1, length = w.length)
  a.hat.pickands <- numeric(w.length)
  x.for.pick.orig <- 1/x
  y.for.pick.orig <- 1/y
  x.for.pick <- x.for.pick.orig/mean(x.for.pick.orig)
  y.for.pick <- y.for.pick.orig/mean(y.for.pick.orig)
  n.for.pick <- length(x.for.pick)
  for(i in seq(w.length)) {
    a.hat.pickands[i] <- 1/mean(pmin(x.for.pick/(1 - w[i]),
                                     y.for.pick/w[i]))
  }

  if(convex.hull){
    if(verbose)
      cat("Taking the convex hull of the estimate ...\n")
    convexPoints <- grDevices::chull(w,a.hat.pickands)
    return(invisible(list(x = w[sort(convexPoints)  ],
                          y = a.hat.pickands[sort(convexPoints)])))
  }else{
    invisible(list(x = w, y = a.hat.pickands))
  }
}
