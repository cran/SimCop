#' Plot a bivariate copula or its density
#'
#' Plot a bivariate copula or its density function.  Essentially hands over immediately to \code{\link[rgl]{persp3d}}.
#'
#' @param x an object of \code{\link{class}} \sQuote{SimCop}.
#' @param \dots additional arguments passed to \code{\link[rgl]{persp3d}}.
#'
#' @keywords hplot
#'
#' @author Berwin A Turlach \email{berwin.turlach@@gmail.com}
#'
#' @method plot SimCop
#' @export
plot.SimCop <- function(x, ...)
  persp3d(x, ...)

#' @rdname plot.SimCop
#'
#' @details If \code{type} is \dQuote{\code{pdf}}, then \code{x} must have a function with name \dQuote{\code{pdfCopula}} in its environment which evaluates the probability density function of the copula stored in \code{x}.
#'
#' For plotting the copula, the \dQuote{x}-grid and \dQuote{y}-grid are produced by using, respectively,  \code{nx} and \code{ny} equispaced points between 0 and 1 (inclusive).  To avoid evaluating the density function at the boundary, if \code{type} is \dQuote{\code{pdf}}, then grids are initially generated in the same manner, after which the first grid point is removed and all remaining grid points are shifted towards the origin by half the distance between two neighbouring grid points.
#'
#' The function to be plotted is evaluated at all possible combination of the grid points.  As the function values of the density function can be large, in particular for extreme value copulae, the values are truncated using either \code{cut} (absolute value) or \code{qcut} (corresponding quantile of all values as determined by the \code{\link{quantile}} function).
#'
#' Finally, the colour scheme passed by argument \code{col} is used to assign a colour to every point at which the function was evaluated by linearly mapping the function values onto the provided colours.
#'
#' @param type specifies whether the copula (\dQuote{\code{cdf}}) or its probability density function (\dQuote{\code{pdf}}) is plotted; see \sQuote{Details}. Can be abbreviated.
#' @param nx,ny length of grids that define the gridlines on which the plotted function is evaluated.
#' @param col colour scheme to use for the surface; see \sQuote{Details}.
#' @param qcut,cut used if \code{type} is \dQuote{\code{pdf}}; see \sQuote{Details}.
#' @param xlab,ylab default labels for the two axes.
#'
#' @return Returns a list with four components \link[=invisible]{invisibly}.  Components \code{x} and \code{y} contain the location of grid lines at which the function was evaluate.  Component \code{z} contains the function evaluations (possibly truncated) and component \code{col} contains the col used.
#'
#' @examples
#' Cop1 <- NewMVFrankCopula(2)
#' plot(Cop1)
#' plot(Cop1, type = "p", qcut = 1)
#'
#' @importFrom rgl persp3d
#' @importFrom grDevices heat.colors
#'
#' @method persp3d SimCop
#' @export
persp3d.SimCop <- function(x, ...,
                           type = c("cdf", "pdf"),
                           nx = 129, ny=129,
                           col = heat.colors(100),
                           qcut = 0.95, cut,
                           xlab=expression(u[1]),
                           ylab=expression(u[2])){
  type <- match.arg(type)
  fct <- x
  x <- seq(from=0, to=1, length=nx)
  y <- seq(from=0, to=1, length=ny)
  if(type == "pdf"){
    if( exists("pdfCopula", envir=environment(fct)) ){
      fct <- get("pdfCopula", envir=environment(fct))
      if(!missing(cut) && !missing(qcut) )
        stop("Please specify 'qcut' or 'cut', but not both.")

      x <- x[-1] - (x[2]-x[1])/2
      y <- y[-1] - (y[2]-y[1])/2
    }else{
      stop("'x' does not have a pdf in its environment.")
    }
  }

  z <- matrix(0, nrow=length(x), ncol=length(y))
  for(i in seq_along(x))
    for(j in seq_along(y))
      z[i,j] <- fct(c(x[i],y[j]))

  if( type == "cdf"){
    cut <- 1
  }else{
    if(missing(cut))
      cut <- stats::quantile(z, qcut, na.rm=TRUE)
  }

  z[z>cut] <- cut
  ncol <- length(col)
  col1 <- col[trunc((ncol-1)*z/max(z, na.rm=TRUE))+1]

  persp3d(x,y,z, ..., xlab=xlab, ylab=ylab, col=col1)
  invisible(list(x=x, y=y, z=z, col=col1))
}
