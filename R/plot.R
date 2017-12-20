#' Plot the histogram density approximation to a copula
#'
#' Plots the histogram density approximation to a copula as determined by \code{\link{GetApprox}}.  Currently works only for bivariate copulae.
#'
#' @param x an object of \code{\link{class}} \sQuote{CopApprox}.
#' @param \dots not used.
#'
#' @references Tajvidi, N. and Turlach, B.A. (2017). A general approach to generate random variates for multivariate copulae, \emph{Australian & New Zealand Journal of Statistics}. Doi:10.1111/anzs.12209.
#'
#' @examples
#' Cop <- NewMEVGumbelCopula(4)
#' CopApprox1 <- GetApprox(Cop, dim=2)
#' plot(CopApprox1)
#' CopApprox2 <- GetApprox(Cop, dim=2, type=2)
#' plot(CopApprox2)
#'
#' @author  Berwin A. Turlach <berwin.turlach@gmail.com>
#'
#' @keywords hplot
#' 
#' @export
plot.CopApprox <- function(x, ...){
  oldpin <- graphics::par("pin")
  newpin <- rep(min(oldpin),2)
  graphics::par(pin=newpin)
  on.exit(graphics::par(pin=oldpin))

  if(x$type == 1){
    d <- NCOL(x$res)/2
    if(d != 2 ) stop("Can only plot 2-dim approximations")

    graphics::plot(c(0,1), c(0,1), type="n", xlab="", ylab="")
    for(i in 1:NROW(x$res)){
      tt <- x$res[i,]
      graphics::polygon(cbind(c(rep(tt[c(1,3)], each=2), tt[1]),
                              c(tt[c(2,4)], tt[c(4,2)], tt[2])))
    }
  }else{
    d <- NCOL(x$res)
    if(d != 2 ) stop("Can only plot 2-dim approximations")
    gc <- grDevices::grey.colors(100, start=0, end=0.8)

    graphics::plot(c(0,1), c(0,1), type="n", xlab="", ylab="")
    h <- x$h
    xh <- c(0, 0, h, h, 0)
    yh <- c(0, h, h, 0, 0)
    cols <- rank(x$probs)
    cols <- 100-trunc(cols/max(cols)*99)
    for(i in 1:NROW(x$res)){
      tt <- x$res[i,]
      graphics::polygon(cbind(tt[1]+xh, tt[2]+yh), col=gc[cols[i]], border=NA)
    }
  }
  invisible()
}
