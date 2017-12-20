#' Cross-derivatives via automatic differentiation
#'
#' Basic building blocks for evaluating functionals \eqn{f:R^d \to R}{f:R^d -> R} and all their cross-derivatives at a given point \eqn{x \in R^d}{x in R^d}.
#'
#' If the optional argument \code{j} is specfied, then the function \eqn{f(x)=x_j} and all its cross-derivatives (all of which but one will be zero, the derivative with respect to the \eqn{j}th component being 1) are evaluated with \eqn{x_j} being set to the value of \code{x}.
#'
#' If the optional argument \code{j} is not used, then the function \eqn{f(x) =c} and all its cross-derivatives (all of which will be zero) are evaluated with \eqn{c} beting set to the value of \code{x}.
#'
#' From these primitive function evaluations, more complicated functions can be constructed using the operations documented in \code{\link{CrossSum}}.
#'
#' @seealso \code{\link{CrossSum}}
#'
#' @param x (scalar) value at which the function is evaluated.
#' @param j optional input. See \sQuote{Details}.
#' @param dim dimension \eqn{d} of the input vector, defaults to two.
#'
#' @return \code{NewCube} returns an object of \code{\link{class}} \sQuote{ADCube} according to its inputs. See \sQuote{Details}.
#'
#' @author Berwin A. Turlach <berwin.turlach@gmail.com>
#'
#' @references Griewank, A., Lehmann, L., Leovey, H. and Zilberman, M. (2014). Automatic evaluations of cross-derivatives, \emph{Mathematics of Computation} \strong{83}(285): 251-274.
#'
#' @keywords internal
#'
NewCube <- function(x, j, dim=2){

  object <- new.env(parent=.GlobalEnv)

  object$val <- rep(0,2^dim)
  object$val[1] <- x
  object$dim <- dim
  if(!missing(j)){
    if(j <= dim){
      object$val[2^(j-1)+1] <- 1
    }else{
      stop("j is larger than dim")
    }
  }
  class(object) <- "ADCube"

  object
}
