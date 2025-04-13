#' Plot the histogram density approximation to a copula
#'
#' Plots the histogram density approximation to a copula as determined by \code{\link{GetApprox}}.  Currently works only for bivariate copulae.
#'
#' If \code{type} is \dQuote{\code{original}} then plots are produced of the kind shown in Tajvidi and Turlach (2017).  In this case the dot arguments are passed to \code{\link{plot}} when the initial plot is created.  Thus, they can be used to set the main title, axes labels and so forth.  If the approximation is of type Approximation II, then argument \code{col} is used to colour the squares while they are plotted and filled via \code{\link{polygon}}.  For this, the ranks of the probability masses of the squares are mapped (linearly) onto the provided colours.
#'
#' If \code{type} is \dQuote{\code{rgl}} then plots are used using \code{\link[rgl]{rgl}}.  The code used is based on the \sQuote{\code{hist3d}} demo of the \code{rgl} package.  The argument \code{col} sets the background colour of the plot.  Arguments \code{topcol} and \code{sidecol} are used to set the colour of the top and sides of the cuboids, and the edges of the cuboids are drawn using \code{linecol}.  Argument \code{alpha} sets the transparency.   Finally, as the heights of the cuboids can be large, in particular for extreme value copulae, their heights are truncated using either \code{cut} (absolute value) or \code{qcut} (corresponding quantile of all heights as determined by the \code{\link{quantile}} function).  After truncation, the heights are rescaled to be between 0 and 1, thus the unit on the \dQuote{z}-axis is meaningless.
#'
#' @param x an object of \code{\link{class}} \sQuote{CopApprox}.
#' @param type specifies the type of plot to produce.  Possible values are \dQuote{\code{rgl}} and \dQuote{\code{original}}; see \sQuote{Details}. Can be abbreviated.
#' @param col colour(s) to be used; see \sQuote{Details}.
#' @param qcut,cut used if \code{type} is \dQuote{\code{rgl}}; see \sQuote{Details}.
#' @param alpha,topcol,sidecol,linecol used if \code{type} is \dQuote{\code{rgl}}; see \sQuote{Details}.
#' @param \dots used if \code{type} is \dQuote{\code{original}}; see \sQuote{Details}.
#'
#' @return \code{NULL} is returned \link[=invisible]{invisibly}.
#'
#' @references Tajvidi, N. and Turlach, B.A. (2017). A general approach to generate random variates for multivariate copulae, \emph{Australian & New Zealand Journal of Statistics} \bold{60}(1): 140--155. \doi{10.1111/anzs.12209}.
#'
#' @examples
#' Cop <- NewMEVGumbelCopula(4)
#' CopApprox1 <- GetApprox(Cop, dim=2)
#' plot(CopApprox1)
#' plot(CopApprox1, type = "o")
#' CopApprox2 <- GetApprox(Cop, dim=2, type=2)
#' plot(CopApprox2)
#' plot(CopApprox2, type = "o", xlab = expression(u[1]), ylab = expression(u[2]))
#' plot(CopApprox2, type = "o", col = heat.colors(100))
#'
#' @author  Berwin A. Turlach \email{berwin.turlach@@gmail.com}
#'
#' @keywords hplot
#'
#' @importFrom grDevices grey.colors
#'
#' @export
plot.CopApprox <- function(x, type = c("rgl", "original"),
                           col = if(type == "rgl" ) "#cccccc"
                                 else grey.colors(100, start=0, end=0.8),
                           qcut = 0.95, cut,
                           alpha = 1,
                           topcol = "#ff0000",
                           sidecol = "#aaaaaa",
                           linecol = "#000000",
                           ...){

  type <- match.arg(type)

  if( type == "rgl" ){
    tmp <- x$res
    if( x$type == 1 ){
      mss <- 1/NROW(tmp)
      areas <- apply(tmp[,3:4] - tmp[,1:2], 1, prod)
    }else{
      mss <- x$probs
      h <- x$h
      areas <- h^2
    }
    heights <- mss/areas
    if(missing(cut)){
      cut <- stats::quantile(heights, qcut)
    }else if(!missing(qcut))
      stop("Please specify 'qcut' or 'cut', but not both.")

    heights[heights > cut] <- cut
    heights <- heights/cut

    rgl::open3d()
    rgl::bg3d(col="#cccccc")

    save <- rgl::par3d(skipRedraw=TRUE)
    on.exit(rgl::par3d(save))

    if(x$type == 1){
      for(i in 1:NROW(tmp)){
        SimCop_binplot_3d(tmp[i,c(1,3)], tmp[i,c(2,4)], heights[i],
                          alpha = alpha,
                          topcol = topcol, sidecol = sidecol, linecol = linecol)
      }
    }else{
      for(i in 1:NROW(tmp)){
        SimCop_binplot_3d(tmp[i,1]+c(0,h), tmp[i,2]+c(0,h), heights[i],
                          alpha = alpha,
                          topcol = topcol, sidecol = sidecol, linecol = linecol)
      }
    }
    rgl::axes3d(c('x','y'))
    rgl::axis3d('z', labels=FALSE, tick=FALSE, nticks=0)
    rgl::title3d(xlab = 'x', ylab= 'y', zlab='z')
  }else{
    oldpin <- graphics::par("pin")
    newpin <- rep(min(oldpin),2)
    graphics::par(pin=newpin)
    on.exit(graphics::par(pin=oldpin))

    dot.args <- list(...)
    names.da <- names(dot.args)
    if ( !("xlab" %in% names.da) )  dot.args$xlab = ""
    if ( !("ylab" %in% names.da) )  dot.args$ylab = ""
    dot.args$x <- dot.args$y <- c(0, 1)
    dot.args$type <- "n"
    do.call(graphics::plot, dot.args)

    if(x$type == 1){
      d <- NCOL(x$res)/2
      if(d != 2 ) stop("Can only plot 2-dim approximations")

      for(i in 1:NROW(x$res)){
        tt <- x$res[i,]
        graphics::polygon(cbind(c(rep(tt[c(1,3)], each=2), tt[1]),
                                c(tt[c(2,4)], tt[c(4,2)], tt[2])))
      }
    }else{
      d <- NCOL(x$res)
      if(d != 2 ) stop("Can only plot 2-dim approximations")

      h <- x$h
      xh <- c(0, 0, h, h, 0)
      yh <- c(0, h, h, 0, 0)
      num.c <- length(col)
      cols <- rank(x$probs)
      cols <- num.c - trunc(cols/(max(cols)+1)*num.c)
      for(i in 1:NROW(x$res)){
        tt <- x$res[i,]
        graphics::polygon(cbind(tt[1]+xh, tt[2]+yh), col=col[cols[i]],
                          border=NA)
      }
    }
  }
  invisible()
}

SimCop_binplot_3d<-function(x, y, z, alpha = alpha,
                            topcol = topcol, sidecol = sidecol,
                            linecol = linecol)
{
#  save <- rgl::par3d(skipRedraw=TRUE)
#  on.exit(rgl::par3d(save))

  x1 <- c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
  z1 <- c(rep(0,4),rep(c(0,0,z,z),4))
  y1 <- c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
  x2 <- c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2))
  z2 <- c(rep(c(0,z),4),rep(0,8),rep(z,8) )
  y2 <- c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2) )
  rgl::rgl.quads(x1,y1,z1,col=rep(sidecol,each=4),alpha=alpha)
  rgl::rgl.quads(c(x[1],x[2],x[2],x[1]),c(y[1],y[1],y[2],y[2]),rep(z,4),
                 col=rep(topcol, each=4), alpha = 1)
  rgl::rgl.lines(x2, y2, z2, col= linecol)
}
