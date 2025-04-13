#' GenerateRV generic
#'
#' A generic to sample random variates from an object.
#'
#' @param obj object from which to sample.
#' @param n   number of items to sample.
#' @param \dots further arguments for methods.
#'
#' @author  Berwin A. Turlach \email{berwin.turlach@@gmail.com}
#'
#' @keywords distribution
#'
#' @export
GenerateRV <- function(obj, n, ...) UseMethod("GenerateRV")
