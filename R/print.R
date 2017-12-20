#' Print basic information on a copula
#'
#' Prints basic information on a copula created with the methods in this package.
#'
#' @param x an object of \code{\link{class}} \sQuote{SimCop}.
#' @param \dots not used.
#'
#' @author  Berwin A. Turlach <berwin.turlach@gmail.com>
#'
#' @keywords print
#'
#' @export
print.SimCop <- function(x, ...){
  param <- get("param", envir=environment(x))
  type <- get("type", envir=environment(x))
  if(length(param) == 0){
    cat(paste0(type[2], " copula of class ", type[1], ", with no parameters:"), "\n")
  }else{
    cat(paste0(type[2], " copula of class ", type[1], ", with parameters:"), "\n")
    for(i in 1:length(param)){
      cat(names(param)[i], "\n")
      print(param[[i]])
    }
  }
  invisible(x)
}

#' Print method for cubes
#'
#' \code{\link{print}} method for \sQuote{ADCube} objects.
#'
#' Currently, the component \sQuote{\code{val}} of the \sQuote{ADCube} object is printed using \code{base::print}.
#'
#'
#' @param x an object of \code{\link{class}} \sQuote{ADCube}.
#' @param \dots arguments passed on to \code{base::print}.
#'
#' @keywords internal
#'
#' @export
print.ADCube <- function(x, ...)
  print(x$val, ...)

