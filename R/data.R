#' Extreme temperatures at two West Australian meteorological stations
#'
#' A dataset on maximum annual values of average daily temperature measurements at two meteorological stations---Leonora (latitude 28.53S, longitude 121.19E) and Menzies (latitude 29.42S, longitude 121.02E)---in Western Australia, for the period 1898--1993.
#'
#' @format A data frame with 96 rows and 2 variables:
#' \describe{
#'   \item{Leonora}{annual maximal temperature at Leonora, in degrees Celsius}
#'   \item{Menzies}{annual maximal temperature at Menzies, in degrees Celsius}
#' }
#'
#' @examples
#' plot(Menzies ~ Leonora, MaxTemp,
#'      xlab = expression("Temperature at Leonora ("*degree*"C)"),
#'      ylab = expression("Temperature at Menzies ("*degree*"C)"))
#'
#' @references Hall, P. and Tajvidi, N. (2004). Prediction regions for bivariate extreme events. \emph{Australian & New Zealand Journal of Statistics} \bold{46}(1): 99--112. \doi{10.1111/j.1467-842X.2004.00316.x}.
#'
#' @keywords datasets
#'
"MaxTemp"
