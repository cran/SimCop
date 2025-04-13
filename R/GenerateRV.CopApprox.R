#' Generates random variates from a copula approximation
#'
#' Method to sample random variates from an object of \code{\link{class}} \sQuote{CopApprox}.
#'
#' If argument \code{MH} is \code{FALSE}, the default, random variates are directly sampled from the approximation, as described in Tajvidi and Turlach (2017).
#'
#' If \code{MH} is \code{TRUE}, \code{GenerateRV} uses additionally a Metropolis-Hastings refinement.  It first samples from the approximation, but uses those samples then as proposals in a Metropolis-Hastings algorithm.  The latter needs the probability density function of the copula.  This density function has either to be passed to the argument \code{PDF}, or the copula (stored in argument \code{obj}) belonging to the approximation must have the density function (with name \sQuote{\code{pdfCopula}}) stored in its environment. In the latter case, the argument \code{PDF} can be \code{NULL} (its default).
#'
#' @inheritParams GenerateRV
#' @param MH logical, should a Metropolis-Hastings algorithm be used to refine the sample?  Default is \code{FALSE}.
#' @param trace   logical, indicating whether the function should be verbose.
#' @param PDF  probability density function corresponding to copula used to create \code{obj}, only used if \code{MH} is \code{TRUE}; see \sQuote{Details}.
#' @param burnin   the number of burn-in iterations of the MH sampler, only used if \code{MH} is \code{TRUE}, defaults to 500.
#' @param thinning     the thining parameter, only used if \code{MH} is \code{TRUE}, defaults to 5.
#' @param \dots  not used.
#'
#' @return A matrix of dimension \code{n} times \code{dim}, where \code{dim} is the dimension for which the copula approximation was determined.
#'
#' If \code{MH} was \code{TRUE} the return value has an attribute called \sQuote{\code{AcceptanceRate}}, indicating the fraction of samples that were accepted in the Metropolis-Hastings step.  This fraction is based on all \code{burnin + (n-1)*thinning + 1} samples that are  initially generated from the approximation.
#'
#' @seealso \code{\link{GetApprox}}
#'
#' @references Tajvidi, N. and Turlach, B.A. (2017). A general approach to generate random variates for multivariate copulae, \emph{Australian & New Zealand Journal of Statistics} \bold{60}(1): 140--155. \doi{10.1111/anzs.12209}.
#'
#' @examples
#' cop <- NewBEVAsyMixedModelCopula(theta=1, phi=-0.25)
#' approx1 <- GetApprox(cop)
#' approx2 <- GetApprox(cop, type = 1)
#' sample1 <- GenerateRV(approx1, 100)
#' plot(sample1)
#' sample2 <- GenerateRV(approx2, 100)
#' plot(sample2)
#' sample1 <- GenerateRV(approx1, 50, MH = TRUE, trace = TRUE)
#' plot(sample1)
#' sample2 <- GenerateRV(approx2, 50, MH = TRUE)
#' plot(sample2)
#'
#' @keywords distribution
#'
#' @export
GenerateRV.CopApprox <- function(obj, n,
                            MH = FALSE,
                            trace = FALSE,
                            PDF = NULL,
                            burnin = 500, thinning = 5, ...){

  if(MH){
    if( is.null(PDF) ){
      if( exists("pdfCopula", envir=environment(obj$Copula), inherits=FALSE) ){
        PDF <- get("pdfCopula", envir=environment(obj$Copula))
      }else{
        stop("Cannot use Metropolis-Hastings.\n'PDF' is not provided and 'obj' does not contain function to calculate pdf.\nPlease use 'MH = FALSE'.")
      }
    }

    n0 <- n
    n <- burnin + (n-1)*thinning + 1
  }

  if(trace) cat("\tSampling from proposal distribution.\n")
  if(obj$type==1){
    d <- NCOL(obj$res)/2
    ind <- base::sample.int(NROW(obj$res), n, replace=TRUE)

    tmp <- obj$res[ind,]
    il <- 1:d
    ir <- il+d
    tmp1 <- tmp[,ir]-tmp[,il]
    res <- tmp[,il] + stats::runif(d*n) * tmp1

    logqval <- -rowSums(log(tmp1))
  }else{
    d <- NCOL(obj$res)
    ind <- base::sample.int(NROW(obj$res), n, replace=TRUE, prob=obj$probs)

    res <- obj$res[ind,] + stats::runif(d*n) * obj$h

    logqval <- log(obj$probs[ind])
  }

  if(MH){
    if(trace) cat("\tEvaluating density at sampled points.\n")
    logdf <- log(apply(res, 1, PDF))

    if(trace) cat("\tPerforming M-H steps.\n")
    tmp <- log(stats::runif(n))
    AcceptRate <- n-1
    for(i in 2:n){
      lalphaxy <- logdf[i]-logdf[i-1]+logqval[i-1]-logqval[i]
      if( tmp[i] > lalphaxy ){
        res[i,] <- res[i-1,]
        logqval[i] <- logqval[i-1]
        logdf[i] <- logdf[i-1]
        AcceptRate <- AcceptRate -1
      }
    }

    ii <- burnin + 1 + (0:(n0-1))*thinning
    res <- res[ii,]
    attr(res, "AcceptanceRate") <- AcceptRate/(n-1)
  }

  dimnames(res) <- dimnames(obj$res)
  res
}
