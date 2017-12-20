#' Calculations on cubes
#'
#' Implements several operations of objects of \code{\link{class}} \sQuote{ADCube}.
#'
#' Assume \code{u} and \code{w} contain the evaluations of functions \eqn{f()} and \eqn{g()}, together with evaluations of all their cross-derivatives, at points \eqn{x} and \eqn{y}, respectively.
#'
#'  \code{CrossSum} returns the evaluation of \eqn{f()+g()}, together with all evaluated cross-derivates.
#'
#'  \code{CrossMult} returns the evaluation of \eqn{f()g()}, together with all evaluated cross-derivates.
#'
#'  \code{CrossDivide} returns the evaluation of \eqn{f()/g()}, together with all evaluated cross-derivates.
#'
#'  \code{CrossSquare} returns the evaluation of \eqn{(f())^2}, together with all evaluated cross-derivates.
#'
#'  \code{CrossExp} returns the evaluation of \eqn{\exp(f())}{exp(f())}, together with all evaluated cross-derivates.
#'
#'  \code{CrossLog} returns the evaluation of \eqn{\log(f())}{log(f())}, together with all evaluated cross-derivates.
#'
#'  \code{CrossPow} returns the evaluation of \eqn{(f())^r}, together with all evaluated cross-derivates.
#'
#'  \code{CrossSqrt} returns the evaluation of \eqn{\sqrt{f()}}{sqrt(f())}, together with all evaluated cross-derivates.
#'
#'  \code{CrossRaisePow} returns the evaluation of \eqn{f()^{g()}}, together with all evaluated cross-derivates.
#'
#' @param u an object of \code{\link{class}} \sQuote{ADCube}.
#' @param w an object of \code{\link{class}} \sQuote{ADCube}.
#'
#' @return An object of \code{\link{class}} \sQuote{ADCube} is returned.  See \sQuote{Details}.
#'
#' @author Berwin A. Turlach <berwin.turlach@gmail.com>
#'
#' @references Griewank, A., Lehmann, L., Leovey, H. and Zilberman, M. (2014). Automatic evaluations of cross-derivatives, \emph{Mathematics of Computation} \strong{83}(285): 251-274.
#'
#' @seealso \code{\link{NewCube}}
#'
#' @keywords internal
#'
CrossSum <- function(u,w){

  if(u$dim != w $dim)
    stop("u and w are not compatible")
  v <- NewCube(0, dim=u$dim)

  v$val <- u$val+w$val
  v
}

#' @rdname CrossSum
#' @inheritParams CrossSum
CrossMult <- function(u,w){

  if(u$dim != w $dim)
    stop("u and w are not compatible")
  v <- NewCube(0, dim=u$dim)
  h <- 2^u$dim

  .crossmult(h, u, w, v, 1, 1, 1)
  v
}

#' @rdname CrossSum
#' @inheritParams CrossSum
CrossDivide <- function(u,w){

  if(u$dim != w $dim)
    stop("u and w are not compatible")
  v <- NewCube(0, dim=u$dim)
  h <- 2^u$dim

  .crossdivide(h, u, w, v, 1, 1, 1)
  v
}

#' @rdname CrossSum
#' @inheritParams CrossSum
CrossSquare <- function(u){

  v <- NewCube(0, dim=u$dim)
  v$val[1] <- u$val[1]*u$val[1]
  h <- 2^v$dim

  i <- 1
  while(i < h){
    .crossmult(i, u, u, v, 1, i+1, i+1)
    i <- i*2
  }
  v$val[-1] <- 2*v$val[-1]

  v
}

#' @rdname CrossSum
#' @inheritParams CrossSum
CrossExp <- function(u){

  v <- NewCube(0, dim=u$dim)

  v$val[1] <- exp(u$val[1])
  h <- 2^v$dim

  i <- 1
  while(i < h){
    .crossmult(i, v, u, v, 1, i+1, i+1)

    i <- i*2
  }
  v
}

#' @rdname CrossSum
#' @inheritParams CrossSum
CrossLog <- function(u){

  v <- NewCube(0, dim=u$dim)

  u0 <- u$val[1]
  u$val[1] <- 0

  v$val[1] <- log(u0)
  v$val[-1] <- u$val[-1] <- u$val[-1]/u0

  h <- 2^v$dim
  i <- 1
  while(i<h){
    .crossdeconv(i, u, v, v, 1, i+1, i+1)
    i <- 2*i
  }

  u$val[1] <- u0
  u$val[-1] <- u$val[-1]*u0
  v
}

#' @rdname CrossSum
#' @inheritParams CrossSum
#' @param r a positive real number.
CrossPow <- function(u, r){

  v <- NewCube(0, dim=u$dim)

  u0 <- u$val[1]
  u$val[1] <- 0

  v$val[1] <- u0^r
  u$val[-1] <- u$val[-1]/u0

  h <- 2^v$dim
  i <- 1
  while(i<h){
    .crossmult(i, v, u, v, 1, i+1, i+1)
    j <- (i+1):(2*i)
    v$val[j] <- v$val[j] * r
    .crossdeconv(i, u, v, v, 1, i+1, i+1)
    i <- 2*i
  }

  u$val[1] <- u0
  u$val[-1] <- u$val[-1]*u0
  v
}

#' @rdname CrossSum
#' @inheritParams CrossSum
CrossSqrt <- function(u){

  v <- NewCube(0, dim=u$dim)

  v$val[-1] <- 0.5*u$val[-1]/u$val[1]

  h <- 2^v$dim
  i <- 1
  while(i<h){
    .crossdeconv(i, v, v, v, 1, i+1, i+1)
    i <- 2*i
  }
  v$val[1] <- sqrt(u$val[1])
  v$val[-1] <- v$val[-1] * v$val[1]
  v
}

#' @rdname CrossSum
#' @inheritParams CrossSum
CrossRaisePow <- function(u,w){

  if(u$dim != w $dim)
    stop("u and w are not compatible")

  tmp1 <- CrossLog(u)
  tmp2 <- CrossMult(w, tmp1)
  CrossExp(tmp2)
}
